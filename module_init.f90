module soil_init

use machine, only : kind_phys

contains

  subroutine noahmp_soil_init(lsm_cold_start          , & ! in
                              im                      , & ! in
                              lsoil_lsm               , & ! in
                              lsoil_input             , & ! in
                              soil_depth_input        , & ! in
                              soil_depth_output       , & ! in
                              soil_moisture_input     , & ! in
                              soil_liquid_input       , & ! in
                              soil_temperature_input  , & ! in
                              soil_type               , & ! in
                              soil_moisture_output    , & ! out
                              soil_liquid_output      , & ! out
                              soil_temperature_output , & ! out
                              errmsg                  , & ! out
                              errflg                  )   ! out


      use noahmp_tables, only: smcref_table, smcdry_table, bexp_table, psisat_table, smcmax_table

      implicit none

      logical,                                          intent(in   ) :: lsm_cold_start
      integer,                                          intent(in   ) :: im
      integer,                                          intent(in   ) :: lsoil_lsm
      integer,                                          intent(in   ) :: lsoil_input
      real (kind=kind_phys), dimension(lsoil_input),    intent(in   ) :: soil_depth_input
      real (kind=kind_phys), dimension(lsoil_lsm),      intent(in   ) :: soil_depth_output
      real (kind=kind_phys), dimension(im,lsoil_input), intent(in   ) :: soil_moisture_input
      real (kind=kind_phys), dimension(im,lsoil_input), intent(in   ) :: soil_liquid_input
      real (kind=kind_phys), dimension(im,lsoil_input), intent(in   ) :: soil_temperature_input

      integer,               dimension(im),             intent(in   ) :: soil_type

      real (kind=kind_phys), dimension(im,lsoil_lsm),   intent(inout) :: soil_moisture_output
      real (kind=kind_phys), dimension(im,lsoil_lsm),   intent(inout) :: soil_liquid_output
      real (kind=kind_phys), dimension(im,lsoil_lsm),   intent(inout) :: soil_temperature_output

      character(len=*),                                 intent(out  ) :: errmsg
      integer,                                          intent(out  ) :: errflg

!> local

      real (kind=kind_phys), dimension(   0:lsoil_input+1)            :: interp_levels 
      real (kind=kind_phys), dimension(im,0:lsoil_input+1)            :: soil_moisture_interp
      real (kind=kind_phys), dimension(im,0:lsoil_input+1)            :: soil_liquid_interp
      real (kind=kind_phys), dimension(im,0:lsoil_input+1)            :: soil_temperature_interp
      real (kind=kind_phys), dimension(     lsoil_input)              :: level_bottom

      integer :: iloc, ilev, lhave, lwant
      real (kind=kind_phys) :: porosity, bexp, psisat, soil_matric_potential, supercool_water
      
      real (kind=kind_phys), parameter :: latent_heat_fusion = 0.3336e06
      real (kind=kind_phys), parameter :: temperature_freezing = 273.16
      real (kind=kind_phys), parameter :: gravity = 9.80616


! Initialize the CCPP error handling variables

      errmsg = ''
      errflg = 0

! interp_levels includes the top(0m) and bottom of the input soil column

      level_bottom(1) = 2.0 * soil_depth_input(1)
      do ilev = 2, lsoil_input
        level_bottom(ilev) = level_bottom(ilev-1) + 2.0 * (soil_depth_input(ilev) - level_bottom(ilev-1))
      end do

      interp_levels(0)             = 0.0
      interp_levels(1:lsoil_input) = soil_depth_input
      interp_levels(lsoil_input+1) = level_bottom(lsoil_input)
      
      soil_temperature_interp(:,1:lsoil_input) = soil_temperature_input
      soil_moisture_interp   (:,1:lsoil_input) = soil_moisture_input
      soil_liquid_interp     (:,1:lsoil_input) = soil_liquid_input

! Linear creation of top and bottom boundary temperature and moisture

      do iloc = 1 , im
        soil_temperature_interp(iloc,0) = soil_temperature_interp(iloc,1) + &
                 ( interp_levels(0) - interp_levels(1) ) * &
                 ( soil_temperature_interp(iloc,1) - soil_temperature_interp(iloc,2) ) /  &
                    ( interp_levels(1) - interp_levels(2) )
        soil_temperature_interp(iloc,lsoil_input+1) = soil_temperature_interp(iloc,lsoil_input) - &
                 ( interp_levels(lsoil_input) - interp_levels(lsoil_input+1) ) * &
                 ( soil_temperature_interp(iloc,lsoil_input-1) - soil_temperature_interp(iloc,lsoil_input) ) /  &
                    ( interp_levels(lsoil_input-1) - interp_levels(lsoil_input) )

        soil_moisture_interp(iloc,0) = soil_moisture_interp(iloc,1) + &
                 ( interp_levels(0) - interp_levels(1) ) * &
                 ( soil_moisture_interp(iloc,1) - soil_moisture_interp(iloc,2) ) /  &
                    ( interp_levels(1) - interp_levels(2) )
        soil_moisture_interp(iloc,lsoil_input+1) = soil_moisture_interp(iloc,lsoil_input) - &
                 ( interp_levels(lsoil_input) - interp_levels(lsoil_input+1) ) * &
                 ( soil_moisture_interp(iloc,lsoil_input-1) - soil_moisture_interp(iloc,lsoil_input) ) /  &
                    ( interp_levels(lsoil_input-1) - interp_levels(lsoil_input) )
      end do

! Interpolate temperature and moisture to wanted levels

      level_want1 : do lwant = 1, lsoil_lsm

        if( soil_depth_output(lwant) > interp_levels(lsoil_input+1) ) then   ! output_depth below input_depths
          do iloc = 1 , im
            soil_temperature_output(iloc,lwant) = soil_temperature_input(iloc,lsoil_input+1)
            soil_moisture_output   (iloc,lwant) = soil_moisture_input(iloc,lsoil_input+1)
          end do
          exit level_want1
        end if

        level_have1 : do lhave = 0 , lsoil_input
            if( ( soil_depth_output(lwant) >= interp_levels(lhave  ) ) .and. &
                ( soil_depth_output(lwant) <= interp_levels(lhave+1) ) ) then   ! output_depth between input_depths
               do iloc = 1 , im
                 soil_temperature_output(iloc,lwant) = &
                    ( soil_temperature_input(iloc,lhave  ) * ( interp_levels(lhave+1) - soil_depth_output(lwant) ) + &
                      soil_temperature_input(iloc,lhave+1) * ( soil_depth_output(lwant) - interp_levels(lhave)) )  / &
                                                          ( interp_levels(lhave+1) - interp_levels(lhave) )
                 soil_moisture_output(iloc,lwant) = &
                    ( soil_moisture_input(iloc,lhave  ) * ( interp_levels(lhave+1) - soil_depth_output(lwant) ) + &
                      soil_moisture_input(iloc,lhave+1) * ( soil_depth_output(lwant) - interp_levels(lhave)) )  / &
                                                          ( interp_levels(lhave+1) - interp_levels(lhave) )
               end do
               exit level_have1
            end if
        end do level_have1
      end do level_want1

! Initialize liquid soil moisture from total soil moisture and soil temperature

      do iloc = 1 , im
        porosity = smcmax_table(soil_type(iloc))
        bexp     = bexp_table(soil_type(iloc))
        psisat   = psisat_table(soil_type(iloc))
      do ilev = 1 , lsoil_lsm
        if(soil_temperature_output(iloc,ilev) >= temperature_freezing) then
          soil_liquid_output(iloc,ilev) = soil_moisture_output(iloc,ilev)
        else
          soil_matric_potential = latent_heat_fusion * (temperature_freezing - soil_temperature_output(iloc,ilev)) / &
                                   (gravity * soil_temperature_output(iloc,ilev))
                supercool_water = porosity*(soil_matric_potential/psisat)**(-1./bexp)
          soil_liquid_output(iloc,ilev) = supercool_water
        end if
      end do
      end do
      
  end subroutine noahmp_soil_init

end module soil_init
