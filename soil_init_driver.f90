program soil_init_driver

use soil_init
use noahmp_tables
use machine, only : kind_phys

implicit none

integer, parameter                              :: im = 1           ! number of spatial grids
integer, parameter                              :: lsoil_input = 3  ! input number of layers
integer, parameter                              :: lsoil_lsm   = 8  ! output number of layers
logical                                         :: lsm_cold_start = .true.

real (kind=kind_phys), dimension(lsoil_input)   :: soil_depth_input        ! input depths to interface [cm]
real (kind=kind_phys), dimension(lsoil_lsm)     :: soil_depth_output       ! output depths to interface [cm]
real (kind=kind_phys), dimension(lsoil_lsm)     :: soil_thickness_output   ! output thicknesses [m]
real (kind=kind_phys), dimension(im,lsoil_input):: soil_moisture_input     ! input soil moisture
real (kind=kind_phys), dimension(im,lsoil_input):: soil_temperature_input  ! input soil temperature
real (kind=kind_phys), dimension(im,lsoil_input):: soil_liquid_input       ! input soil liquid
integer,               dimension(im)            :: soil_type               ! soil type

real (kind=kind_phys), dimension(im,lsoil_lsm)  :: soil_moisture_output    ! output soil moisture
real (kind=kind_phys), dimension(im,lsoil_lsm)  :: soil_temperature_output ! output soil temperature
real (kind=kind_phys), dimension(im,lsoil_lsm)  :: soil_liquid_output      ! output soil liquid

character(len=128)                              :: errmsg
integer                                         :: errflg
integer                                         :: ilev

soil_type                    = 6
!soil_depth_input             = (/0.05, 0.25, 0.70, 1.50/)
!soil_depth_output            = (/0.05, 0.25, 0.70, 1.50/)
soil_depth_input             = (/0.5,1.5,2.5/)
soil_depth_output            = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/)
soil_moisture_input(1,:)     = (/0.5,0.5,0.5/)
soil_liquid_input(1,:)       = (/0.3,0.3,0.3/)
soil_temperature_input(1,:)  = (/100.0,200.0,300.0/)

call read_mp_table_parameters(errmsg, errflg)

call noahmp_soil_init (lsm_cold_start          , & ! in
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

print *, "input soil temperature:"
do ilev = 1, lsoil_input
 print *, soil_depth_input(ilev), soil_temperature_input(1,ilev)
end do

print *, "output soil temperature:"
do ilev = 1, lsoil_lsm
 print *, soil_depth_output(ilev), soil_temperature_output(1,ilev)
end do

print *
print *, "============================"
print *

print *, "input soil moisture:"
do ilev = 1, lsoil_input
 print *, soil_depth_input(ilev), soil_moisture_input(1,ilev)
end do

print *, "output soil moisture:"
do ilev = 1, lsoil_lsm
 print *, soil_depth_output(ilev), soil_moisture_output(1,ilev)
end do

print *
print *, "============================"
print *

print *, "input soil liquid:"
do ilev = 1, lsoil_input
 print *, soil_depth_input(ilev), soil_liquid_input(1,ilev)
end do

print *, "output soil liquid:"
do ilev = 1, lsoil_lsm
 print *, soil_depth_output(ilev), soil_liquid_output(1,ilev)
end do

end program
