!-------------------------------------------------------------------------------
! (C) British Crown Copyright 2022 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------

!> Fortran module to handle scatwind observations - Met Office neutral operator

module ufo_scatwind_neutralmetoffice_tlad_mod_c

  use fckit_configuration_module, only: fckit_configuration
  use ufo_scatwind_neutralmetoffice_tlad_mod
  use ufo_geovals_mod,    only: ufo_geovals
  use ufo_geovals_mod_c,  only: ufo_geovals_registry

  implicit none
  private
  
#define LISTED_TYPE ufo_scatwind_neutralmetoffice_tlad
  
  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
  
  !> Global registry
  type(registry_t) :: ufo_scatwind_neutralmetoffice_tlad_registry
  
  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  
! ------------------------------------------------------------------------------
  
subroutine ufo_scatwind_neutralmetoffice_tlad_setup_c(c_key_self,         &
                                                      surface_type_check, &
                                                      surface_type_sea,   &
                                                      c_obsvars,          &
                                                      c_geovars,          &
                                                      c_nchan,            &
                                                      c_channels,         &
                                                      c_conf) bind(c,name='ufo_scatwind_neutralmetoffice_tlad_setup_f90')
use oops_variables_mod
use obs_variables_mod
implicit none
integer(c_int), intent(inout)  :: c_key_self
logical(c_bool), intent(in)    :: surface_type_check !< Whether to check we are over sea before calculating H(x)
integer(c_int), intent(in)     :: surface_type_sea   !< Value for sea to be used in the surface type check
type(c_ptr), intent(in), value :: c_obsvars !< variables to be simulated
type(c_ptr), intent(in), value :: c_geovars !< variables requested from the model

type(c_ptr), value, intent(in) :: c_conf
integer(c_size_t), intent(in)  :: c_nchan
integer(c_int), intent(in)     :: c_channels(c_nchan)

type(fckit_configuration) :: f_conf

type(ufo_scatwind_NeutralMetOffice_tlad), pointer :: self

call ufo_scatwind_neutralmetoffice_tlad_registry%setup(c_key_self, self)
f_conf = fckit_configuration(c_conf)
self%obsvars = obs_variables(c_obsvars)
self%geovars = oops_variables(c_geovars)

call self%setup(c_channels, surface_type_check, surface_type_sea)

end subroutine ufo_scatwind_neutralmetoffice_tlad_setup_c
  
! ------------------------------------------------------------------------------
  
subroutine ufo_scatwind_neutralmetoffice_tlad_delete_c(c_key_self) bind(c,name='ufo_scatwind_neutralmetoffice_tlad_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
    
type(ufo_scatwind_NeutralMetOffice_tlad), pointer :: self

call ufo_scatwind_NeutralMetOffice_tlad_registry%delete(c_key_self,self)

end subroutine ufo_scatwind_neutralmetoffice_tlad_delete_c

! ------------------------------------------------------------------------------

subroutine ufo_scatwind_neutralmetoffice_tlad_settraj_c(c_key_self, c_key_geovals, c_obsspace) bind(c,name='ufo_scatwind_neutralmetoffice_tlad_settraj_f90')
implicit none
integer(c_int), intent(in)     :: c_key_self
integer(c_int), intent(in)     :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace

type(ufo_scatwind_NeutralMetOffice_tlad), pointer :: self
type(ufo_geovals),                        pointer :: geovals

character(len=*), parameter :: myname_="ufo_scatwind_neutralmetoffice_tlad_settraj_c"

call ufo_scatwind_neutralmetoffice_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%settraj(geovals, c_obsspace)

end subroutine ufo_scatwind_neutralmetoffice_tlad_settraj_c

! ------------------------------------------------------------------------------

subroutine ufo_scatwind_neutralmetoffice_simobs_tl_c(c_key_self, c_key_geovals, &
                                                     c_obsspace, c_nvars, c_nlocs, &
                                                     c_hofx) bind(c,name='ufo_scatwind_neutralmetoffice_simobs_tl_f90')

implicit none
integer(c_int), intent(in)     :: c_key_self
integer(c_int), intent(in)     :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_scatwind_NeutralMetOffice_tlad), pointer :: self
type(ufo_geovals),                        pointer :: geovals
character(len=*), parameter :: myname_="ufo_scatwind_neutralmetoffice_simobs_tl_c"

call ufo_scatwind_NeutralMetOffice_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_tl(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_scatwind_neutralmetoffice_simobs_tl_c

! ------------------------------------------------------------------------------

subroutine ufo_scatwind_neutralmetoffice_simobs_ad_c(c_key_self, c_key_geovals, &
                                                     c_obsspace, c_nvars, c_nlocs, &
                                                     c_hofx) bind(c,name='ufo_scatwind_neutralmetoffice_simobs_ad_f90')

implicit none
integer(c_int), intent(in)     :: c_key_self
integer(c_int), intent(in)     :: c_key_geovals
type(c_ptr), value, intent(in) :: c_obsspace
integer(c_int), intent(in)     :: c_nvars, c_nlocs
real(c_double), intent(inout)  :: c_hofx(c_nvars, c_nlocs)

type(ufo_scatwind_NeutralMetOffice_tlad), pointer :: self
type(ufo_geovals),                        pointer :: geovals
character(len=*), parameter :: myname_="ufo_scatwind_neutralmetoffice_simobs_ad_c"

call ufo_scatwind_NeutralMetOffice_tlad_registry%get(c_key_self, self)
call ufo_geovals_registry%get(c_key_geovals, geovals)

call self%simobs_ad(geovals, c_obsspace, c_nvars, c_nlocs, c_hofx)

end subroutine ufo_scatwind_neutralmetoffice_simobs_ad_c


! ------------------------------------------------------------------------------

end module ufo_scatwind_neutralmetoffice_tlad_mod_c
