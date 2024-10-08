/*
 * (C) Copyright 2017-2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_
#define UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/rttovonedvarcheck/RTTOVOneDVarCheck.interface.h"
#include "ufo/filters/rttovonedvarcheck/RTTOVOneDVarCheckParameters.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {
//! \brief RTTOVOneDVarCheck
//!
//! \details RTTOVOneDVarCheck performs a 1D-Var minimization for satellite using the rttov
//! forward operator.  If a profile does not converge then all channels for this observation
//! are removed.  Some parameters (e.g. surface emissivity) are retrieved and saved in the
//! obsspace for use in 4D-Var.  The code is based on the Met Office 1D-Var scheme and thus
//! is predominently in Fortran.
//!
//! \author Met Office
//!
//! \date 09/06/2020
//!

class RTTOVOneDVarCheck : public FilterBase,
                          private util::ObjectCounter<RTTOVOneDVarCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef RTTOVOneDVarCheckParameters Parameters_;

  static const std::string classname() {return "ufo::RTTOVOneDVarCheck";}

  RTTOVOneDVarCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~RTTOVOneDVarCheck() override;

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::onedvar;}

  F90obfilter keyRTTOVOneDVarCheck_;
  std::vector<int> channels_;
  // TODO(AS): Danger: below should be ObsVariables but need to figure out GeoVals vs ObsDiagnostics
  oops::Variables retrieved_vars_;
  oops::Variables hoxdiags_retrieved_vars_;
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_
