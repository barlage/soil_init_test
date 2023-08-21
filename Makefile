# Makefile 
#

COMPILERF90 = gfortran -ffree-line-length-none -fdefault-real-8 

OBJS =  module_init.o \
        noahmp_tables.o


all:    soil_init_driver.exe

soil_init_driver.exe: $(OBJS)

	$(COMPILERF90) -o $(@) soil_init_driver.f90 $(OBJS)

module_init.o: module_init.f90

	$(COMPILERF90) -c $(*).f90 

noahmp_tables.o: noahmp_tables.f90

	$(COMPILERF90) -c $(*).f90 

machine.o: machine.F

	$(COMPILERF90) -c $(*).F

#
# This command cleans up object (etc) files:
#

clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#

module_init.o: noahmp_tables.o machine.o
noahmp_tables.o: machine.o
