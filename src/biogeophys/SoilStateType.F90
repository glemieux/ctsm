module SoilStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsoi, nlevgrnd, nlevlak, nlayer, nlevsno, nlevmaxurbgrnd
  use clm_varcon      , only : spval
  use clm_varctl      , only : use_hydrstress, use_cn, use_lch4, use_fates
  use clm_varctl      , only : iulog, hist_wrtch4diag
  use LandunitType    , only : lun                
  use ColumnType      , only : col                
  use PatchType       , only : patch                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilstate_type

     ! sand/ clay/ organic matter
     real(r8), pointer :: sandfrac_patch       (:)   ! patch sand fraction
     real(r8), pointer :: clayfrac_patch       (:)   ! patch clay fraction 
     real(r8), pointer :: mss_frc_cly_vld_col  (:)   ! col mass fraction clay limited to 0.20
     real(r8), pointer :: cellorg_col          (:,:) ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellsand_col         (:,:) ! sand value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellclay_col         (:,:) ! clay value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: bd_col               (:,:) ! col bulk density of dry soil material [kg/m^3] (CN)

     ! hydraulic properties
     real(r8), pointer :: hksat_col            (:,:) ! col hydraulic conductivity at saturation (mm H2O /s) 
     real(r8), pointer :: hksat_min_col        (:,:) ! col mineral hydraulic conductivity at saturation (hksat) (mm/s)
     real(r8), pointer :: hk_l_col             (:,:) ! col hydraulic conductivity (mm/s)
     real(r8), pointer :: smp_l_col            (:,:) ! col soil matric potential (mm)
     real(r8), pointer :: smpmin_col           (:)   ! col restriction for min of soil potential (mm) 
     real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)  
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity) 
     real(r8), pointer :: watdry_col           (:,:) ! col btran parameter for btran = 0
     real(r8), pointer :: watopt_col           (:,:) ! col btran parameter for btran = 1
     real(r8), pointer :: watfc_col            (:,:) ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd) 
     real(r8), pointer :: dsl_col              (:)   ! col dry surface layer thickness (mm)
     real(r8), pointer :: soilresis_col        (:)   ! col soil evaporative resistance S&L14 (s/m)
     real(r8), pointer :: soilbeta_col         (:)   ! col factor that reduces ground evaporation L&P1992(-)
     real(r8), pointer :: soilalpha_col        (:)   ! col factor that reduces ground saturated specific humidity (-)
     real(r8), pointer :: soilalpha_u_col      (:)   ! col urban factor that reduces ground saturated specific humidity (-) 
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: wtfact_col           (:)   ! col maximum saturated fraction for a gridcell
     real(r8), pointer :: porosity_col         (:,:) ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col     (:,:) ! col effective porosity = porosity - vol_ice (nlevgrnd) 
     real(r8), pointer :: gwc_thr_col          (:)   ! col threshold soil moisture based on clay content
!scs: vangenuchten
     real(r8), pointer :: msw_col              (:,:) ! col vanGenuchtenClapp "m"
     real(r8), pointer :: nsw_col              (:,:) ! col vanGenuchtenClapp "n"
     real(r8), pointer :: alphasw_col          (:,:) ! col vanGenuchtenClapp "nalpha"
     real(r8), pointer :: watres_col           (:,:) ! residual soil water content
     ! thermal conductivity / heat capacity
     real(r8), pointer :: thk_col              (:,:) ! col thermal conductivity of each layer [W/m-K] 
     real(r8), pointer :: tkmg_col             (:,:) ! col thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: tkdry_col            (:,:) ! col thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
     real(r8), pointer :: tksatu_col           (:,:) ! col thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: csol_col             (:,:) ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 

     ! roots
     real(r8), pointer :: rootr_patch          (:,:) ! patch effective fraction of roots in each soil layer (SMS method only) (nlevgrnd)
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (SMS method only) (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootfr_patch         (:,:) ! patch fraction of roots for water in each soil layer (nlevgrnd)
     real(r8), pointer :: crootfr_patch        (:,:) ! patch fraction of roots for carbon in each soil layer (nlevgrnd)
     real(r8), pointer :: root_depth_patch     (:)   ! root depth
     real(r8), pointer :: rootr_road_perv_col  (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: rootfr_road_perv_col (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: k_soil_root_patch    (:,:) ! patch soil-root interface conductance [mm/s]
     real(r8), pointer :: root_conductance_patch(:,:) ! patch root conductance [mm/s]
     real(r8), pointer :: soil_conductance_patch(:,:) ! patch soil conductance [mm/s]

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type soilstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%mss_frc_cly_vld_col  (begc:endc))                     ; this%mss_frc_cly_vld_col  (:)   = nan
    allocate(this%sandfrac_patch       (begp:endp))                     ; this%sandfrac_patch       (:)   = nan
    allocate(this%clayfrac_patch       (begp:endp))                     ; this%clayfrac_patch       (:)   = nan
    allocate(this%cellorg_col          (begc:endc,nlevsoi))             ; this%cellorg_col          (:,:) = nan 
    allocate(this%cellsand_col         (begc:endc,nlevsoi))             ; this%cellsand_col         (:,:) = nan 
    allocate(this%cellclay_col         (begc:endc,nlevsoi))             ; this%cellclay_col         (:,:) = nan 
    allocate(this%bd_col               (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan

    allocate(this%hksat_col            (begc:endc,nlevgrnd))            ; this%hksat_col            (:,:) = spval
    allocate(this%hksat_min_col        (begc:endc,nlevgrnd))            ; this%hksat_min_col        (:,:) = spval
    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan   
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan   
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%bsw_col              (begc:endc,nlevgrnd))            ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col           (begc:endc,nlevmaxurbgrnd))      ; this%watsat_col           (:,:) = nan
    allocate(this%watdry_col           (begc:endc,nlevgrnd))            ; this%watdry_col           (:,:) = spval
    allocate(this%watopt_col           (begc:endc,nlevgrnd))            ; this%watopt_col           (:,:) = spval
    allocate(this%watfc_col            (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%dsl_col              (begc:endc))                     ; this%dsl_col         (:)   = spval!nan   
    allocate(this%soilresis_col        (begc:endc))                     ; this%soilresis_col         (:)   = spval!nan   
    allocate(this%soilbeta_col         (begc:endc))                     ; this%soilbeta_col         (:)   = nan   
    allocate(this%soilalpha_col        (begc:endc))                     ; this%soilalpha_col        (:)   = nan
    allocate(this%soilalpha_u_col      (begc:endc))                     ; this%soilalpha_u_col      (:)   = nan
    allocate(this%soilpsi_col          (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%wtfact_col           (begc:endc))                     ; this%wtfact_col           (:)   = nan
    allocate(this%porosity_col         (begc:endc,nlayer))              ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col     (begc:endc,nlevgrnd))            ; this%eff_porosity_col     (:,:) = spval
    allocate(this%gwc_thr_col          (begc:endc))                     ; this%gwc_thr_col          (:)   = nan

    allocate(this%thk_col              (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%thk_col              (:,:) = nan
    allocate(this%tkmg_col             (begc:endc,nlevgrnd))            ; this%tkmg_col             (:,:) = nan
    allocate(this%tkdry_col            (begc:endc,nlevgrnd))            ; this%tkdry_col            (:,:) = nan
    allocate(this%tksatu_col           (begc:endc,nlevgrnd))            ; this%tksatu_col           (:,:) = nan
    allocate(this%csol_col             (begc:endc,nlevgrnd))            ; this%csol_col             (:,:) = nan

    allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootr_road_perv_col  (begc:endc,1:nlevgrnd))          ; this%rootr_road_perv_col  (:,:) = nan
    allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%crootfr_patch        (begp:endp,1:nlevgrnd))          ; this%crootfr_patch        (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan 
    allocate(this%rootfr_road_perv_col (begc:endc,1:nlevgrnd))          ; this%rootfr_road_perv_col (:,:) = nan
    allocate(this%k_soil_root_patch    (begp:endp,1:nlevsoi))           ; this%k_soil_root_patch (:,:) = nan
    allocate(this%root_conductance_patch(begp:endp,1:nlevsoi))          ; this%root_conductance_patch (:,:) = nan
    allocate(this%soil_conductance_patch(begp:endp,1:nlevsoi))          ; this%soil_conductance_patch (:,:) = nan
    allocate(this%msw_col              (begc:endc,1:nlevgrnd))          ; this%msw_col              (:,:) = nan
    allocate(this%nsw_col              (begc:endc,1:nlevgrnd))          ; this%nsw_col              (:,:) = nan
    allocate(this%alphasw_col          (begc:endc,1:nlevgrnd))          ; this%alphasw_col          (:,:) = nan
    allocate(this%watres_col           (begc:endc,1:nlevgrnd))          ; this%watres_col           (:,:) = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc


    if (use_lch4) then
       if (hist_wrtch4diag) then
          active = "active"
       else
          active = "inactive"
       end if
    else
       active = "inactive"
    end if

    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential (natural vegetated and crop landunits only)', &
         ptr_col=this%smp_l_col, set_spec=spval, l2g_scale_type='veg')

    this%root_conductance_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='KROOT', units='1/s', type2d='levsoi', &
         avgflag='A', long_name='root conductance each soil layer', &
         ptr_patch=this%root_conductance_patch, default='inactive')
    
    this%soil_conductance_patch(begp:endp,:) = spval
    call hist_addfld2d (fname='KSOIL', units='1/s', type2d='levsoi', &
         avgflag='A', long_name='soil conductance in each soil layer', &
         ptr_patch=this%soil_conductance_patch, default='inactive')

    if (use_cn) then
       this%bsw_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='bsw', units='unitless', type2d='levgrnd', &
            avgflag='A', long_name='clap and hornberger B', &
            ptr_col=this%bsw_col, default='inactive')
    end if

    if (use_cn) then
       this%rootr_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer (SMS method)', &
            ptr_patch=this%rootr_patch, l2g_scale_type='veg', default='inactive')
    end if

    if (use_cn .and. .not.(use_hydrstress)) then
       ! rootr_col isn't computed for use_hydrstress = .true. (In contrast, rootr_patch is
       ! still computed, albeit using the inconsistent Soil Moisture Stress (SMS) method.)
       ! (See also https://github.com/ESCOMP/CTSM/issues/812.)
       this%rootr_col(begc:endc,:) = spval
       call hist_addfld2d (fname='ROOTR_COLUMN', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer (SMS method)', &
            ptr_col=this%rootr_col, l2g_scale_type='veg', default='inactive')
       
    end if

    if (use_cn .or. use_fates) then
       this%soilpsi_col(begc:endc,:) = spval
       call hist_addfld2d (fname='SOILPSI', units='MPa', type2d='levgrnd', &
            avgflag='A', long_name='soil water potential in each soil layer', &
            ptr_col=this%soilpsi_col, default='inactive')
    end if

    this%thk_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%thk_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_TK', units='W/m-K', type2d='levsno', &
         avgflag='A', long_name='Thermal conductivity', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_TK_ICE', units='W/m-K', type2d='levsno', &
         avgflag='A', long_name='Thermal conductivity (ice landunits only)', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%hk_l_col(begc:endc,:) = spval
    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity (natural vegetated and crop landunits only)', &
         ptr_col=this%hk_l_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

    this%soilalpha_col(begc:endc) = spval
    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=this%soilalpha_col, set_urb=spval, default='inactive' )

    this%soilalpha_u_col(begc:endc) = spval
    call hist_addfld1d (fname='SoilAlpha_U',  units='unitless',  &
         avgflag='A', long_name='urban factor limiting ground evap', &
         ptr_col=this%soilalpha_u_col, set_nourb=spval, default='inactive')

    if (use_cn) then
       this%watsat_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='watsat', units='m^3/m^3', type2d='levgrnd', &
            avgflag='A', long_name='water saturated', &
            ptr_col=this%watsat_col, default='inactive')
    end if

    if (use_cn) then
       this%eff_porosity_col(begc:endc,:) = spval
       call hist_addfld2d (fname='EFF_POROSITY', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective porosity = porosity - vol_ice', &
            ptr_col=this%eff_porosity_col, default='inactive')
    end if

    if (use_cn) then
       this%watfc_col(begc:endc,:) = spval 
       call hist_addfld2d (fname='watfc', units='m^3/m^3', type2d='levgrnd', &
            avgflag='A', long_name='water field capacity', &
            ptr_col=this%watfc_col, default='inactive')
    end if

    this%soilresis_col(begc:endc) = spval
    call hist_addfld1d (fname='SOILRESIS',  units='s/m',  &
         avgflag='A', long_name='soil resistance to evaporation', &
         ptr_col=this%soilresis_col)

    this%dsl_col(begc:endc) = spval
    call hist_addfld1d (fname='DSL',  units='mm',  &
         avgflag='A', long_name='dry surface layer thickness', &
         ptr_col=this%dsl_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module soil state variables to reasonable values
    !
    ! !USES:
    use clm_varpar      , only : nlevgrnd
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    this%smp_l_col(bounds%begc:bounds%endc,1:nlevgrnd) = -1000._r8
    this%hk_l_col(bounds%begc:bounds%endc,1:nlevgrnd) = 0._r8

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    use spmdMod          , only : masterproc
    use RootBiophysMod   , only : init_vegrootfr
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c
    logical  :: readvar
    logical  :: readrootfr = .false.
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='DSL', xtype=ncd_double,  &
         dim1name='column', long_name='dsl thickness', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%dsl_col)
    
    call restartvar(ncid=ncid, flag=flag, varname='SOILRESIS', xtype=ncd_double,  &
         dim1name='column', long_name='soil resistance', units='s/m', &
         interpinic_flag='interp', readvar=readvar, data=this%soilresis_col)

    call restartvar(ncid=ncid, flag=flag, varname='SMP', xtype=ncd_double,  &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='soil matric potential', units='mm', &
         scale_by_thickness=.true., &
         interpinic_flag='interp', readvar=readvar, data=this%smp_l_col)

    call restartvar(ncid=ncid, flag=flag, varname='HK', xtype=ncd_double,  &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='hydraulic conductivity', units='mm/s', &
         scale_by_thickness=.true., &
         interpinic_flag='interp', readvar=readvar, data=this%hk_l_col)

     readrootfr = .false.
     if (flag=='read' .and. .not. readrootfr ) then
            if (masterproc) then
               write(iulog,*) "can't find rootfr in restart (or initial) file..."
               write(iulog,*) "Initialize rootfr to default"
            end if
            call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
            this%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd), 'water')
            call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
            this%crootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd), 'carbon')
         end if
    
  end subroutine Restart

end module SoilStateType
