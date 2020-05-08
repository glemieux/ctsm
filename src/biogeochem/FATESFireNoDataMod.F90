module FATESFireNoDataMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for FATES when not obtaining fire inputs from data
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod, only: errmsg => shr_log_errMsg
  use abortutils, only: endrun
  use clm_varctl, only: iulog
  use decompMod, only: bounds_type
  use CNFireMethodMod, only: cnfire_method_type
  use CNFireBaseMod, only: cnfire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_no_data_type
  !
  type, extends(cnfire_base_type) :: fates_fire_no_data_type

      ! !PRIVATE MEMBER DATA:

    contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: InitAccBuffer  ! Initialize accumulation processes
      procedure, public :: InitAccVars  ! Initialize accumulation variables
      procedure, public :: UpdateAccVars  ! Update/extract accumulations vars

  end type fates_fire_no_data_type

  character(len=*), parameter, private :: sourcefile = __FILE__
  !
  ! !PRIVATE MEMBER DATA:
  !-----------------------------------------------------------------------

  interface fates_fire_no_data_type
     ! initialize a new cnfire_base object
     module procedure constructor
  end interface fates_fire_no_data_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  type(fates_fire_no_data_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type fates_fire_no_data_type
    ! !ARGUMENTS:
    constructor%need_lightning_and_popdens = .false.
  end function constructor

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart
    ! file is read in and the accumulation buffer is obtained)
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! !DESCRIPTION:
    ! EMPTY subroutine for the no_data case.
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_no_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------


  end subroutine UpdateAccVars

end module FATESFireNoDataMod