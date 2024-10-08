/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/TrackCheck.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/sqr.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/TrackCheckParameters.h"
#include "ufo/filters/TrackCheckUtils.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

TrackCheck::TrackObservation::TrackObservation(float latitude, float longitude,
                                               const util::DateTime &time, float pressure,
                                               float speed, int64_t timeOffset)
  : obsLocationTime_(latitude, longitude, time),  pressure_(pressure),
    speed_(speed), timeOffset_(timeOffset),
    rejectedInPreviousSweep_(false), rejectedBeforePreviousSweep_(false),
    numNeighborsVisitedInPreviousSweep_{NO_PREVIOUS_SWEEP, NO_PREVIOUS_SWEEP}
{}

void TrackCheck::TrackObservation::checkAgainstBuddy(
    const TrackObservation &buddyObs,
    const TrackCheckParameters &options,
    float maxSpeed,
    CheckResults & results) const {
  results = CheckResults();
  const float temporalDistance =
    std::abs(static_cast<float>(buddyObs.timeOffset() - this->timeOffset()));

  const float spatialDistance = TrackCheckUtils::distance(this->obsLocationTime_.location(),
                                                          buddyObs.obsLocationTime_.location());

  const float temporalResolutionSeconds =
    static_cast<float>(options.temporalResolution.value().toSeconds());

  // Estimate the speed and check if it is within the allowed range
  const float conservativeSpeedEstimate =
    (spatialDistance - options.spatialResolution) /
    (temporalDistance + temporalResolutionSeconds);

  results.speedCheckResult = TrackCheckUtils::CheckResult(conservativeSpeedEstimate <= maxSpeed);

  // Estimate the climb rate and check if it is within the allowed range
  if (options.maxClimbRate.value() != boost::none) {
    const float pressureDiff = std::abs(pressure_ - buddyObs.pressure_);
    const float conservativeClimbRateEstimate =
      pressureDiff / (temporalDistance + temporalResolutionSeconds);
    results.climbRateCheckResult =
        TrackCheckUtils::CheckResult
        (conservativeClimbRateEstimate <= *options.maxClimbRate.value());
  }

  const int resolutionMultiplier = options.distinctBuddyResolutionMultiplier;
  results.isBuddyDistinct =
      temporalDistance > resolutionMultiplier * temporalResolutionSeconds &&
      spatialDistance > resolutionMultiplier * options.spatialResolution;

  return;
}

void TrackCheck::TrackObservation::registerCheckResults(const CheckResults &results) {
  checkCounter_.registerCheckResult(results.speedCheckResult);
  checkCounter_.registerCheckResult(results.climbRateCheckResult);
}

void TrackCheck::TrackObservation::unregisterCheckResults(const CheckResults &results) {
  checkCounter_.unregisterCheckResult(results.speedCheckResult);
  checkCounter_.unregisterCheckResult(results.climbRateCheckResult);
}

void TrackCheck::TrackObservation::registerSweepOutcome(bool rejectedInSweep) {
  rejectedBeforePreviousSweep_ = rejected();
  rejectedInPreviousSweep_ = rejectedInSweep;
}

float TrackCheck::TrackObservation::getFailedChecksFraction() {
  return this->checkCounter_.failedChecksFraction();
}


TrackCheck::TrackCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), options_(parameters)
{
  oops::Log::debug() << "TrackCheck: config = " << options_ << std::endl;
}

// Required for the correct destruction of options_.
TrackCheck::~TrackCheck()
{}

void TrackCheck::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  // 3rd arg: recordsAreSingleObs = false for Track Check.
  ObsAccessor obsAccessor = TrackCheckUtils::createObsAccessor(options_.stationIdVariable,
                                                               obsdb_,
                                                               false);

  const std::vector<size_t> validObsIds =
    options_.ignoreExistingQCFlags ?
    obsAccessor.getValidObservationIds(apply) :
    obsAccessor.getValidObservationIds(apply, *flags_, filtervars);
  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  TrackCheckUtils::sortTracksChronologically(validObsIds, obsAccessor, splitter);

  PiecewiseLinearInterpolation maxSpeedByPressure = makeMaxSpeedByPressureInterpolation();
  ObsGroupPressureLocationTime obsPressureLocTime
    = collectObsPressuresLocationsTimes(obsAccessor, maxSpeedByPressure);

  std::vector<bool> isRejected(obsPressureLocTime.pressures.size(), false);
  for (auto track : splitter.multiElementGroups()) {
    identifyRejectedObservationsInTrack(track.begin(), track.end(), validObsIds,
                                        obsPressureLocTime, isRejected);
  }
  obsAccessor.flagRejectedObservations(isRejected, flagged);
}

TrackCheck::ObsGroupPressureLocationTime TrackCheck::collectObsPressuresLocationsTimes
(const ObsAccessor &obsAccessor,
 const PiecewiseLinearInterpolation &maxValidSpeedAtPressure) const {
  ObsGroupPressureLocationTime obsPressureLocTime;
  obsPressureLocTime.locationTimes = TrackCheckUtils::collectObservationsLocations(obsAccessor);
  obsPressureLocTime.pressures = obsAccessor.getFloatVariableFromObsSpace(options_.pressureGroup,
                                                                      options_.pressureCoord);
  // Compute speed from pressure and offset relative to epoch.
  const util::DateTime epoch(1970, 1, 1, 0, 0, 0);
  for (size_t jloc = 0; jloc < obsPressureLocTime.pressures.size(); ++jloc) {
    obsPressureLocTime.speeds.push_back(maxValidSpeedAtPressure
                                        (obsPressureLocTime.pressures[jloc]));
    obsPressureLocTime.timeOffsets.push_back((obsPressureLocTime.locationTimes.datetimes[jloc] -
                                          epoch).toSeconds());
  }

  return obsPressureLocTime;
}

PiecewiseLinearInterpolation TrackCheck::makeMaxSpeedByPressureInterpolation() const {
  const std::map<float, float> &maxSpeedInterpolationPoints =
      options_.maxSpeedInterpolationPoints.value();

  std::vector<double> pressures, maxSpeeds;
  pressures.reserve(maxSpeedInterpolationPoints.size());
  maxSpeeds.reserve(maxSpeedInterpolationPoints.size());

    // The interpolator needs to produce speeds in km/s rather than m/s because observation
    // locations are expressed in kilometers.
    const int metersPerKm = 1000;
  for (const auto &pressureAndMaxSpeed : maxSpeedInterpolationPoints) {
    pressures.push_back(pressureAndMaxSpeed.first);
    maxSpeeds.push_back(pressureAndMaxSpeed.second / metersPerKm);
  }

  return PiecewiseLinearInterpolation(std::move(pressures), std::move(maxSpeeds));
}

void TrackCheck::identifyRejectedObservationsInTrack(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const ObsGroupPressureLocationTime &obsPressureLocTime,
    std::vector<bool> &isRejected) const {

  std::vector<TrackObservation> trackObservations = collectTrackObservations(
        trackObsIndicesBegin, trackObsIndicesEnd, validObsIds, obsPressureLocTime);

  std::vector<float> workspace;
  while (sweepOverObservations(trackObservations, workspace) ==
         TrackCheckUtils::SweepResult::ANOTHER_SWEEP_REQUIRED) {
    // can't exit the loop yet
  }

  flagRejectedTrackObservations(trackObsIndicesBegin, trackObsIndicesEnd,
                                validObsIds, trackObservations, isRejected);
}

std::vector<TrackCheck::TrackObservation> TrackCheck::collectTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const ObsGroupPressureLocationTime &obsPressureLocTime) const {
  std::vector<TrackObservation> trackObservations;
  trackObservations.reserve(trackObsIndicesEnd - trackObsIndicesBegin);
  for (std::vector<size_t>::const_iterator it = trackObsIndicesBegin;
       it != trackObsIndicesEnd; ++it) {
    const size_t obsId = validObsIds[*it];
    trackObservations.emplace_back(obsPressureLocTime.locationTimes.latitudes[obsId],
                                   obsPressureLocTime.locationTimes.longitudes[obsId],
                                   obsPressureLocTime.locationTimes.datetimes[obsId],
                                   obsPressureLocTime.pressures[obsId],
                                   obsPressureLocTime.speeds[obsId],
                                   obsPressureLocTime.timeOffsets[obsId]);
  }
  return trackObservations;
}

TrackCheckUtils::SweepResult TrackCheck::sweepOverObservations(
    std::vector<TrackObservation> &trackObservations,
    std::vector<float> &workspace) const {

  TrackCheck::CheckResults results;

  std::vector<float> &failedChecksFraction = workspace;
  failedChecksFraction.assign(trackObservations.size(), 0.0f);

  for (int obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
    TrackObservation &obs = trackObservations[obsIdx];
    if (obs.rejected())
      continue;

    for (Direction dir : { FORWARD, BACKWARD}) {
      const bool firstSweep = obs.numNeighborsVisitedInPreviousSweep(dir) == NO_PREVIOUS_SWEEP;
      const int numNeighborsVisitedInPreviousSweep =
          firstSweep ? 0 : obs.numNeighborsVisitedInPreviousSweep(dir);
      int numNewDistinctBuddiesToVisit = firstSweep ? options_.numDistinctBuddiesPerDirection : 0;

      auto getNthNeighbor = [&trackObservations, obsIdx, dir](int n) -> const TrackObservation* {
        const int neighborObsIdx = obsIdx + (dir == FORWARD ? n : -n);
        if (neighborObsIdx < 0 || neighborObsIdx >= trackObservations.size())
          return nullptr;  // We've reached the end of the track
        else
          return &trackObservations[neighborObsIdx];
      };

      float minPressureBetween = obs.pressure();
      float maxSpeed = obs.speed();
      int neighborIdx = 1;
      const TrackObservation *neighborObs = getNthNeighbor(neighborIdx);
      for (; neighborIdx <= numNeighborsVisitedInPreviousSweep && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        // Strictly speaking, neighborObs->pressure should be disregarded if neighborObs has already
        // been rejected. However, that would force us to check each pair of observations anew
        // whenever an observation between them is rejected, whereas as things stand, we only
        // need to "undo" checks against rejected observations.

        if (minPressureBetween >= neighborObs->pressure())
          maxSpeed = neighborObs->speed();
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure());

        if (neighborObs->rejectedInPreviousSweep()) {
          obs.checkAgainstBuddy(*neighborObs, options_, maxSpeed, results);
          obs.unregisterCheckResults(results);
          if (results.isBuddyDistinct) {
            // The rejected distinct buddy needs to be replaced with another
            ++numNewDistinctBuddiesToVisit;
          }
        }
      }

      for (; numNewDistinctBuddiesToVisit > 0 && neighborObs != nullptr;
           neighborObs = getNthNeighbor(++neighborIdx)) {
        if (minPressureBetween >= neighborObs->pressure())
          maxSpeed = neighborObs->speed();
        minPressureBetween = std::min(minPressureBetween, neighborObs->pressure());

        if (!neighborObs->rejected()) {
          obs.checkAgainstBuddy(*neighborObs, options_, maxSpeed, results);
          obs.registerCheckResults(results);
          if (results.isBuddyDistinct)
            --numNewDistinctBuddiesToVisit;
        }
      }

      const int numNeighborsVisitedInThisSweep = neighborIdx - 1;
      assert(numNeighborsVisitedInThisSweep >= numNeighborsVisitedInPreviousSweep);
      obs.setNumNeighborsVisitedInPreviousSweep(dir, numNeighborsVisitedInThisSweep);
    }  // end of loop over directions

    failedChecksFraction[obsIdx] = obs.getFailedChecksFraction();
  }

  const float maxFailedChecksFraction = *std::max_element(failedChecksFraction.begin(),
                                                          failedChecksFraction.end());
  const float failedChecksThreshold = options_.rejectionThreshold * maxFailedChecksFraction;
  if (failedChecksThreshold <= 0)
    return TrackCheckUtils::SweepResult::NO_MORE_SWEEPS_REQUIRED;

  for (int obsIdx = 0; obsIdx < trackObservations.size(); ++obsIdx) {
    const bool rejected = failedChecksFraction[obsIdx] > failedChecksThreshold;
    trackObservations[obsIdx].registerSweepOutcome(rejected);
  }

  return TrackCheckUtils::SweepResult::ANOTHER_SWEEP_REQUIRED;
}

void TrackCheck::flagRejectedTrackObservations(
    std::vector<size_t>::const_iterator trackObsIndicesBegin,
    std::vector<size_t>::const_iterator trackObsIndicesEnd,
    const std::vector<size_t> &validObsIds,
    const std::vector<TrackObservation> &trackObservations,
    std::vector<bool> &isRejected) const {
  auto trackObsIndexIt = trackObsIndicesBegin;
  auto trackObsIt = trackObservations.begin();
  for (; trackObsIndexIt != trackObsIndicesEnd; ++trackObsIndexIt, ++trackObsIt)
    if (trackObsIt->rejected())
      isRejected[validObsIds[*trackObsIndexIt]] = true;
}

void TrackCheck::print(std::ostream & os) const {
  os << "TrackCheck: config = " << options_ << std::endl;
}

}  // namespace ufo
