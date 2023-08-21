program soil_init_driver

use soil_init
use noahmp_tables
use machine, only : kind_phys

implicit none

integer, parameter                             :: im = 1          ! number of spatial grids
integer, parameter                             :: lsoil = 4       ! input number of layers
integer, parameter                             :: lsoil_lsm = 4   ! output number of layers
logical                                        :: lsm_cold_start = .false.

real (kind=kind_phys), dimension(im)           :: tskin_lnd   ! surface temperature
real (kind=kind_phys), dimension(im)           :: tg3         ! deep soil temperature
real (kind=kind_phys), dimension(1:lsoil)      :: zsi         ! input depths to interface
real (kind=kind_phys), dimension(1:lsoil_lsm)  :: zso         ! output depths to interface
real (kind=kind_phys), dimension(1:lsoil_lsm)  :: dzso        ! output thicknesses
real (kind=kind_phys), dimension(im,lsoil)     :: smc         ! input soil moisture
real (kind=kind_phys), dimension(im,lsoil)     :: stc         ! input soil temperature
real (kind=kind_phys), dimension(im,lsoil)     :: slc         ! input soil liquid
integer,               dimension(im)           :: stype       ! input soil type
integer,               dimension(im)           :: vtype       ! input vegetation type

real (kind=kind_phys), dimension(im,lsoil_lsm) :: smois       ! output soil moisture
real (kind=kind_phys), dimension(im,lsoil_lsm) :: tslb        ! output soil temperature
real (kind=kind_phys), dimension(im,lsoil_lsm) :: sh2o        ! output soil liquid

character(len=128)                             :: errmsg
integer                                        :: errflg

integer, parameter :: me      = huge(1)
integer, parameter :: isot    = 1
integer, parameter :: ivegsrc = 1
integer, parameter :: nlunit  = huge(1)

stype     = 6
vtype     = -99999
tskin_lnd = 300.0
tg3       = 300.0
zsi       = (/ 5, 25, 70, 150/)
zso       = (/ 5, 25, 70, 150/)
dzso      = (/0.1, 0.3, 0.6, 1.0/)
smc(1,:)  = (/0.3,0.3,0.3,0.3/)
slc(1,:)  = (/0.3,0.3,0.3,0.3/)
stc(1,:)  = (/300.0,300.0,300.0,300.0/)

call read_mp_table_parameters(errmsg, errflg)

call noahmpsoilinit (lsm_cold_start, im, lsoil_lsm, lsoil,    & ! in
                     zsi,zso, dzso,tskin_lnd,tg3, smc, slc, stc,    & ! in
                     sh2o,tslb, smois, stype, vtype,          & ! in
                     errmsg, errflg)

print *, "input soil temperature:"
print *, stc
print *, "input surface temperature:"
print *, tskin_lnd
print *, "input deep soil temperature:"
print *, tg3

print *, "output soil temperature:"
print *, tslb

print *
print *, "============================"
print *

print *, "input soil moisture:"
print *, smc

print *, "output soil moisture:"
print *, smois

print *
print *, "============================"
print *

print *, "input soil liquid:"
print *, slc

print *, "output soil liquid:"
print *, sh2o

end program
