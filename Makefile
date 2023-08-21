# Makefile 
#

COMPILERF90 = gfortran -ffree-line-length-none -fdefault-real-8 

OBJS =  module_init.o \
        module_mp_soil_init.o \
        namelist_soilveg.o \
        set_soilveg.o \
        noahmp_tables.o


all:    soil_init_driver.exe

soil_init_driver.exe: $(OBJS)

	$(COMPILERF90) -o $(@) soil_init_driver.f90 $(OBJS)

module_init.o: module_init.f90

	$(COMPILERF90) -c $(*).f90 

module_mp_soil_init.o: module_mp_soil_init.f90

	$(COMPILERF90) -c $(*).f90 

noahmp_tables.o: noahmp_tables.f90

	$(COMPILERF90) -c $(*).f90 

namelist_soilveg.o: namelist_soilveg.f

	$(COMPILERF90) -c $(*).f

machine.o: machine.F

	$(COMPILERF90) -c $(*).F

set_soilveg.o: set_soilveg.f

	$(COMPILERF90) -c $(*).f

#
# This command cleans up object (etc) files:
#

clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#

module_init.o: module_mp_soil_init.o namelist_soilveg.o noahmp_tables.o machine.o
set_soilveg.o: namelist_soilveg.o
module_mp_soil_init.o: noahmp_tables.o
noahmp_tables.o: machine.o
