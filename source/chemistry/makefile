CHEMFILE=ion_elec_chem
FILE=fundconstants.f90 ${CHEMFILE}.f90
OBJFILE=${CHEMFILE}.o fundconstants.o

CHEMDRIVER: 
	gfortran -c ${FILE}
	gfortran ${OBJFILE} chem_driver.f90 -o chemdriver
clean:
	rm -rf *.o *.mod chemdriver
