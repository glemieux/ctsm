module SoilBiogeochemLittVertTranspMod

  !-----------------------------------------------------------------------
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog, use_c13, use_c14, spinup_state
  use clm_varcon                         , only : secspday
  use decompMod                          , only : bounds_type
  use abortutils                         , only : endrun
  use ActiveLayerMod                     , only : active_layer_type
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenFluxType     , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, use_soil_matrixcn
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  !
  implicit none
  private
  !
  public :: readParams
  public :: SoilBiogeochemLittVertTransp

  type, private :: params_type
     real(r8) :: som_diffus                 ! Soil organic matter diffusion
     real(r8) :: cryoturb_diffusion_k       ! The cryoturbation diffusive constant cryoturbation to the active layer thickness
     real(r8) :: max_altdepth_cryoturbation ! (m) maximum active layer thickness for cryoturbation to occur
  end type params_type

  type(params_type), private :: params_inst
  !
  real(r8), public :: som_adv_flux =  0._r8
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  real(r8) :: som_diffus                   ! [m^2/sec] = 1 cm^2 / yr
  real(r8) :: cryoturb_diffusion_k         ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
  real(r8) :: max_altdepth_cryoturbation   ! (m) maximum active layer thickness for cryoturbation to occur

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------  
  subroutine readParams ( ncid )
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'SoilBiogeochemLittVertTranspType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in parameters
    !

     tString='som_diffus'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     params_inst%som_diffus=tempr

     tString='cryoturb_diffusion_k'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     params_inst%cryoturb_diffusion_k=tempr

     tString='max_altdepth_cryoturbation'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     params_inst%max_altdepth_cryoturbation=tempr
    
   end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemLittVertTransp(bounds, num_bgc_soilc, filter_bgc_soilc,      &
       active_layer_inst, soilbiogeochem_state_inst,                     &
       soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst, &
       c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
       c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools 
    ! calculated in the CStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
    use clm_time_manager , only : get_step_size_real
    use clm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use clm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    use TridiagonalMod   , only : Tridiagonal
    use ColumnType       , only : col
    use clm_varctl       , only : use_bedrock

    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds 
    integer                                 , intent(in)    :: num_bgc_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)  ! filter for soil columns
    type(active_layer_type)                 , intent(in)    :: active_layer_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: diffus (bounds%begc:bounds%endc,1:nlevdecomp+1)                    ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! advective flux (m/s)  (includes spinup correction, if any)
    real(r8) :: aaa                                                                ! "A" function in Patankar
    real(r8) :: pe                                                                 ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1                                                         ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1                                                         ! Harmonic mean of diffusivity
    real(r8) :: a_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "r" vector for tridiagonal solution
    real(r8) :: d_p1_zp1(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: f_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)                       ! water flux for next j
    real(r8) :: f_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)                       ! water flux for previous j
    real(r8) :: pe_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)                      ! Peclet # for next j
    real(r8) :: pe_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)                      ! Peclet # for previous j
    real(r8) :: dz_node(1:nlevdecomp+1)                                            ! difference between nodes
    real(r8) :: epsilon_t (bounds%begc:bounds%endc,1:nlevdecomp+1,1:ndecomp_pools) !
    real(r8) :: conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1)                  !
    real(r8) :: a_p_0
    real(r8) :: deficit
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l,i                                                ! indices
    integer  :: jtop(bounds%begc:bounds%endc)                                      ! top level at each column
    real(r8) :: dtime                                                              ! land model time step (sec)
    integer  :: zerolev_diffus
    real(r8) :: spinup_term                                                        ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8) :: epsilon                                                            ! small number
    real(r8), pointer :: conc_ptr(:,:,:)                                           ! pointer, concentration state variable being transported
    real(r8), pointer :: source(:,:,:)                                             ! pointer, source term
    real(r8), pointer :: trcr_tendency_ptr(:,:,:)                                  ! poiner, store the vertical tendency (gain/loss due to vertical transport)
    real(r8), pointer :: matrix_input(:,:)                                  ! poiner, store the vertical tendency (gain/loss due to vertical transport)
    !-----------------------------------------------------------------------

    ! Set statement functions
    aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! A function from Patankar, Table 5.2, pg 95
  
    associate(                                                             &
         is_cwd           => decomp_cascade_con%is_cwd                  ,  & ! Input:  [logical (:)    ]  TRUE => pool is a cwd pool                                
         spinup_factor    => decomp_cascade_con%spinup_factor           ,  & ! Input:  [real(r8) (:)   ]  spinup accelerated decomposition factor, used to accelerate transport as well

         altmax           => active_layer_inst%altmax_col               ,  & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw                             
         altmax_lastyear  => active_layer_inst%altmax_lastyear_col      ,  & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw                  

         som_adv_coef     => soilbiogeochem_state_inst%som_adv_coef_col ,  & ! Output: [real(r8) (:,:) ]  SOM advective flux (m/s)                               
         som_diffus_coef  => soilbiogeochem_state_inst%som_diffus_coef_col,& ! Output: [real(r8) (:,:) ]  SOM diffusivity due to bio/cryo-turbation (m2/s)  
         tri_ma_vr        => soilbiogeochem_carbonflux_inst%tri_ma_vr &      ! Output: [real(r8) (:,:) ]  Vertical CN transfer rate in sparse matrix format (gC*m3)/(gC*m3*step))
         )

      !Set parameters of vertical mixing of SOM
      som_diffus                 = params_inst%som_diffus 
      cryoturb_diffusion_k       = params_inst%cryoturb_diffusion_k 
      max_altdepth_cryoturbation = params_inst%max_altdepth_cryoturbation 

      dtime = get_step_size_real()

      ntype = 2
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif
      spinup_term = 1._r8
      epsilon = 1.e-30

      !------ first get diffusivity / advection terms -------!
      ! use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
      do fc = 1, num_bgc_soilc
         c = filter_bgc_soilc (fc)
         if  (( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
              ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
            ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
            do j = 1,nlevdecomp+1
               if ( j <= col%nbedrock(c)+1 ) then
                  if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
                     som_diffus_coef(c,j) = cryoturb_diffusion_k 
                     som_adv_coef(c,j) = 0._r8
                  else
                     som_diffus_coef(c,j) = max(cryoturb_diffusion_k * & 
                       ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                       ( min(max_depth_cryoturb, zisoi(col%nbedrock(c)+1)) - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
                     som_adv_coef(c,j) = 0._r8
                  endif
               else
                  som_adv_coef(c,j) = 0._r8
                  som_diffus_coef(c,j) = 0._r8
               endif
            end do
         elseif (  max(altmax(c), altmax_lastyear(c)) > 0._r8 ) then
            ! constant advection, constant diffusion
            do j = 1,nlevdecomp+1
               if ( j <= col%nbedrock(c)+1 ) then
                  som_adv_coef(c,j) = som_adv_flux 
                  som_diffus_coef(c,j) = som_diffus
               else
                  som_adv_coef(c,j) = 0._r8
                  som_diffus_coef(c,j) = 0._r8
               endif
            end do
         else
            ! completely frozen soils--no mixing
            do j = 1,nlevdecomp+1
               som_adv_coef(c,j) = 0._r8
               som_diffus_coef(c,j) = 0._r8
            end do
         endif
      end do

      ! Set the distance between the node and the one ABOVE it   
      dz_node(1) = zsoi(1)
      do j = 2,nlevdecomp+1
         dz_node(j)= zsoi(j) - zsoi(j-1)
      enddo

      !------ loop over litter/som types
      do i_type = 1, ntype

         select case (i_type)
         case (1)  ! C
            conc_ptr          => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col
            source            => soilbiogeochem_carbonflux_inst%decomp_cpools_sourcesink_col
            trcr_tendency_ptr => soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col
            matrix_input      => soilbiogeochem_carbonflux_inst%matrix_Cinput%V
         case (2)  ! N
            conc_ptr          => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col
            source            => soilbiogeochem_nitrogenflux_inst%decomp_npools_sourcesink_col
            trcr_tendency_ptr => soilbiogeochem_nitrogenflux_inst%decomp_npools_transport_tendency_col
            matrix_input      => soilbiogeochem_nitrogenflux_inst%matrix_Ninput%V
         case (3)
            if ( use_c13 ) then
               ! C13
               conc_ptr          => c13_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col
               source            => c13_soilbiogeochem_carbonflux_inst%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c13_soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col
            else
               ! C14
               conc_ptr          => c14_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col
               source            => c14_soilbiogeochem_carbonflux_inst%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c14_soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col
            endif
         case (4)
            if ( use_c14 .and. use_c13 ) then
               ! C14
               conc_ptr          => c14_soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col
               source            => c14_soilbiogeochem_carbonflux_inst%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c14_soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col
            else
               write(iulog,*) 'error.  ncase = 4, but c13 and c14 not both enabled.'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            endif
         end select

         do s = 1, ndecomp_pools
            if ( .not. is_cwd(s) ) then
               if(.not. use_soil_matrixcn .or. s .eq. 1)then
                  do j = 1,nlevdecomp+1
                     do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)
                     !
                     if ( spinup_state >= 1 ) then
                        ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
                        spinup_term = spinup_factor(s)
                     else
                        spinup_term = 1._r8
                     endif

                     if (abs(spinup_term - 1._r8) > .000001_r8 ) then
                        spinup_term = spinup_term * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
                     endif

                     if ( abs(som_adv_coef(c,j)) * spinup_term < epsilon ) then
                        adv_flux(c,j) = epsilon
                     else
                        adv_flux(c,j) = som_adv_coef(c,j) * spinup_term
                     endif
                     !
                     if ( abs(som_diffus_coef(c,j)) * spinup_term < epsilon ) then
                        diffus(c,j) = epsilon
                     else
                        diffus(c,j) = som_diffus_coef(c,j) * spinup_term
                     endif
                     !
                     end do
                  end do

                  ! Set Pe (Peclet #) and D/dz throughout column

                  do fc = 1, num_bgc_soilc ! dummy terms here
                     c = filter_bgc_soilc (fc)
                     conc_trcr(c,0) = 0._r8
                     conc_trcr(c,col%nbedrock(c)+1:nlevdecomp+1) = 0._r8
                  end do


                  do j = 1,nlevdecomp+1
                     do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)

                     conc_trcr(c,j) = conc_ptr(c,j,s)
               
                     ! dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
                     ! dz_node_tracer is difference between cell centers 

                     ! Calculate the D and F terms in the Patankar algorithm
                     if (j == 1) then
                        d_m1_zm1(c,j) = 0._r8
                        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                        if ( diffus(c,j+1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                           d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(c,j) + w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                        else
                           d_p1 = 0._r8
                        endif
                        d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                        f_m1(c,j) = adv_flux(c,j)  ! Include infiltration here
                        f_p1(c,j) = adv_flux(c,j+1)
                        pe_m1(c,j) = 0._r8
                        pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                     elseif (j >= col%nbedrock(c)+1) then
                        ! At the bottom, assume no gradient in d_z (i.e., they're the same)
                        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                        if ( diffus(c,j) > 0._r8 .and. diffus(c,j-1) > 0._r8) then
                           d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(c,j) + w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                        else
                           d_m1 = 0._r8
                        endif
                        d_m1_zm1(c,j) = d_m1 / dz_node(j)
                        d_p1_zp1(c,j) = d_m1_zm1(c,j) ! Set to be the same
                        f_m1(c,j) = adv_flux(c,j)
                        !f_p1(c,j) = adv_flux(c,j+1)
                        f_p1(c,j) = 0._r8
                        pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                        pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                     else
                        ! Use distance from j-1 node to interface with j divided by distance between nodes
                        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                        if ( diffus(c,j-1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                           d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(c,j) + w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                        else
                           d_m1 = 0._r8
                        endif
                        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                        if ( diffus(c,j+1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                           d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(c,j) + w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                        else
                           d_p1 = (1._r8 - w_m1) * diffus(c,j) + w_p1 * diffus(c,j+1) ! Arithmetic mean of diffus
                        endif
                        d_m1_zm1(c,j) = d_m1 / dz_node(j)
                        d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                        f_m1(c,j) = adv_flux(c,j)
                        f_p1(c,j) = adv_flux(c,j+1)
                        pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                        pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                     end if
                     enddo ! fc
                  enddo ! j; nlevdecomp
               end if


               ! Calculate the tridiagonal coefficients
               do j = 0,nlevdecomp +1
                  do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)
                     ! g = cgridcell(c)

                     if (j > 0 .and. j < nlevdecomp+1) then
                        a_p_0 =  dzsoi_decomp(j) / dtime
                     endif

                     if (j == 0) then ! top layer (atmosphere)
                        a_tri(c,j) = 0._r8
                        b_tri(c,j) = 1._r8
                        c_tri(c,j) = -1._r8
                        r_tri(c,j) = 0._r8
                     elseif (j == 1) then
                        a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + max( f_m1(c,j), 0._r8)) ! Eqn 5.47 Patankar
                        c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + max(-f_p1(c,j), 0._r8))
                        b_tri(c,j) = - a_tri(c,j) - c_tri(c,j) + a_p_0
                        r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + (a_p_0 - adv_flux(c,j)) * conc_trcr(c,j)
                        if(s .eq. 1 .and. i_type .eq. 1 .and. use_soil_matrixcn )then !vertical matrix are the same for all pools
                           do i = 1,ndecomp_pools-1 !excluding cwd
                              tri_ma_vr(c,1+(i-1)*(nlevdecomp*3-2)) = (b_tri(c,j) - a_p_0) / dzsoi_decomp(j) * (-dtime)
                              tri_ma_vr(c,3+(i-1)*(nlevdecomp*3-2)) = c_tri(c,j) / dzsoi_decomp(j) * (-dtime)
                           end do
                        end if
                     elseif (j < nlevdecomp+1) then
                        a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + max( f_m1(c,j), 0._r8)) ! Eqn 5.47 Patankar
                        c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + max(-f_p1(c,j), 0._r8))
                        b_tri(c,j) = - a_tri(c,j) - c_tri(c,j) + a_p_0
                        r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + a_p_0 * conc_trcr(c,j)
                        if(s .eq. 1 .and. i_type .eq. 1 .and. use_soil_matrixcn )then                   
                           if(j .le. col%nbedrock(c))then
                              do i = 1,ndecomp_pools-1
                                 tri_ma_vr(c,j*3-4+(i-1)*(nlevdecomp*3-2)) = a_tri(c,j) / dzsoi_decomp(j) * (-dtime)
                                 if(j .ne. nlevdecomp)then
                                    tri_ma_vr(c,j*3  +(i-1)*(nlevdecomp*3-2)) = c_tri(c,j) / dzsoi_decomp(j) * (-dtime)
                                 end if
                                 tri_ma_vr(c,j*3-2+(i-1)*(nlevdecomp*3-2)) = (b_tri(c,j) - a_p_0) / dzsoi_decomp(j) * (-dtime)
                              end do
                           else
                              if(j .eq. col%nbedrock(c) + 1 .and. j .ne. nlevdecomp .and. j .gt. 1)then
                                 do i = 1,ndecomp_pools-1
                                    tri_ma_vr(c,(j-1)*3-2+(i-1)*(nlevdecomp*3-2)) = tri_ma_vr(c,(j-1)*3-2+(i-1)*(nlevdecomp*3-2)) &
                                                                                 + a_tri(c,j) / dzsoi_decomp(j-1)*(-dtime)
                                 end do
                              end if
                           end if
                        end if
                     else ! j==nlevdecomp+1; 0 concentration gradient at bottom
                        a_tri(c,j) = -1._r8
                        b_tri(c,j) = 1._r8
                        c_tri(c,j) = 0._r8 
                        r_tri(c,j) = 0._r8
                     endif
                  enddo ! fc; column
               enddo ! j; nlevdecomp

               do fc = 1, num_bgc_soilc
                  c = filter_bgc_soilc (fc)
                  jtop(c) = 0
               enddo

               ! subtract initial concentration and source terms for tendency calculation
               do fc = 1, num_bgc_soilc
                  c = filter_bgc_soilc (fc)
                  do j = 1, nlevdecomp
                     if (.not. use_soil_matrixcn) then
                        trcr_tendency_ptr(c,j,s) = 0.-(conc_trcr(c,j) + source(c,j,s))
                     else
                        trcr_tendency_ptr(c,j,s) = 0.0_r8
                    end if !soil_matrix 
                 end do
               end do

               if (.not. use_soil_matrixcn) then
                  ! Solve for the concentration profile for this time step
                  call Tridiagonal(bounds, 0, nlevdecomp+1, &
                    jtop(bounds%begc:bounds%endc), &
                    num_bgc_soilc, filter_bgc_soilc, &
                    a_tri(bounds%begc:bounds%endc, :), &
                    b_tri(bounds%begc:bounds%endc, :), &
                    c_tri(bounds%begc:bounds%endc, :), &
                    r_tri(bounds%begc:bounds%endc, :), &
                    conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1))
                  ! add post-transport concentration to calculate tendency term
                  do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)
                     do j = 1, nlevdecomp
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) + conc_trcr(c,j)
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) / dtime
                     end do
                  end do
               else  ! For matrix solution set the matrix input array
                  do j = 1,nlevdecomp
                     do fc =1,num_bgc_soilc
                        c = filter_bgc_soilc(fc)
                        matrix_input(c,j+(s-1)*nlevdecomp) = matrix_input(c,j+(s-1)*nlevdecomp) + source(c,j,s)
                     end do
                  end do
               end if  !soil_matrix
            else
               ! for CWD pools, just add
               do j = 1,nlevdecomp
                  do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)
                     if(.not. use_soil_matrixcn)then
                        conc_trcr(c,j) = conc_ptr(c,j,s) + source(c,j,s)
                     else
                        matrix_input(c,j+(s-1)*nlevdecomp) = matrix_input(c,j+(s-1)*nlevdecomp) + source(c,j,s)
                     end if
                     if (j > col%nbedrock(c) .and. source(c,j,s) > 0._r8) then 
                        write(iulog,*) 'source >0',c,j,s,source(c,j,s)
                     end if
                     if (j > col%nbedrock(c) .and. conc_ptr(c,j,s) > 0._r8) then
                        write(iulog,*) 'conc_ptr >0',c,j,s,conc_ptr(c,j,s)
                     end if
                  end do
               end do
            end if ! not CWD

            if (.not. use_soil_matrixcn) then
               do j = 1,nlevdecomp
                  do fc = 1, num_bgc_soilc
                     c = filter_bgc_soilc (fc)
                     conc_ptr(c,j,s) = conc_trcr(c,j) 
                     ! Correct for small amounts of carbon that leak into bedrock
                     if (j > col%nbedrock(c)) then 
                        conc_ptr(c,col%nbedrock(c),s) = conc_ptr(c,col%nbedrock(c),s) + &
                           conc_trcr(c,j) * (dzsoi_decomp(j) / dzsoi_decomp(col%nbedrock(c)))
                        conc_ptr(c,j,s) = 0._r8
                     end if
                  end do
               end do
            end if !not soil_matrix
         end do ! s (pool loop)

      end do  ! i_type
   
    end associate

  end subroutine SoilBiogeochemLittVertTransp
 
end module SoilBiogeochemLittVertTranspMod
