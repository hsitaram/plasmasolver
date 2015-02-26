CHEMFILEPREFIX=ion_elec_chem
CHEMFILE    = ${CHEMFILEPREFIX}.f90
CHEMOBJFILE = ${CHEMFILEPREFIX}.o
SOURCEDIR = ../source
CHEMDIR  = ${SOURCEDIR}/chemistry

FORTFILES = ${SOURCEDIR}/fundconstants.f90 ${CHEMDIR}/${CHEMFILE} ${SOURCEDIR}/gmres.f90 ${SOURCEDIR}/convdiff.f90 ${SOURCEDIR}/solver.f90
OBJFILES = fundconstants.o gmres.o convdiff.o solver.o ${CHEMOBJFILE}
MAINFILE = ${SOURCEDIR}/main.f90

FORT = gfortran
#FORTFLAGS = -Wunused -fbounds-check
FORTFLAGS = -O3
#FORTFLAGS = 

fortranconvdiff:
	${FORT} -c ${FORTFLAGS} ${FORTFILES}
	${FORT} ${FORTFLAGS} ${OBJFILES} ${MAINFILE} -o plsolver_${CHEMFILEPREFIX}.exe

clean:
	rm -rf *.o *.mod *.exe *.dat
