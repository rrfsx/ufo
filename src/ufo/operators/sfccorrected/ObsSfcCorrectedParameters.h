/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
#define UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// enum type for surface correction type, and ParameterTraitsHelper for it
enum class SfcCorrectionType {
  UKMO, WRFDA, RRFS_GSL
};
struct SfcCorrectionTypeParameterTraitsHelper {
  typedef SfcCorrectionType EnumType;
  static constexpr char enumTypeName[] = "SfcCorrectionType";
  static constexpr util::NamedEnumerator<SfcCorrectionType> namedValues[] = {
    { SfcCorrectionType::UKMO, "UKMO" },
    { SfcCorrectionType::WRFDA, "WRFDA" },
    { SfcCorrectionType::RRFS_GSL, "RRFS_GSL" }
  };
};

}  // namespace ufo

namespace oops {

/// Extraction of SfcCorrectionType parameters from config
template <>
struct ParameterTraits<ufo::SfcCorrectionType> :
    public EnumParameterTraits<ufo::SfcCorrectionTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the SfcCorrected operator.
class ObsSfcCorrectedParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSfcCorrectedParameters, ObsOperatorParametersBase)

 public:
  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this operator is used as a
  /// component of the Composite operator. If `variables` is not set, the operator will simulate
  /// all ObsSpace variables. Please see the documentation of the Composite operator for further
  /// details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
     "variables",
     "List of variables to be simulated",
     this};

  oops::Parameter<SfcCorrectionType> correctionType{"da_sfc_scheme",
     "Scheme used for surface temperature correction (UKMO or WRFDA)",
     SfcCorrectionType::UKMO, this};

  /// Note: "height" default value has to be consistent with var_geomz defined
  /// in ufo_variables_mod.F90
  oops::Parameter<std::string> geovarGeomZ{"geovar_geomz",
     "Model variable for height of vertical levels",
     "height", this};

  /// Note: "surface_altitude" default value has to be consistent with var_sfc_geomz
  /// in ufo_variables_mod.F90
  oops::Parameter<std::string> geovarSfcGeomZ{"geovar_sfc_geomz",
     "Model variable for surface height",
     "surface_altitude", this};

  /// Note: "station_altitude" default value is "stationElevation"
  oops::Parameter<std::string> ObsHeightName{"station_altitude", "stationElevation", this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
