PETSC_DIR=/opt/petsc-3.5.0/
#PETSC_ARCH=linux-gnu
CFLAGS = -I.
FFLAGS=
CPPFLAGS=-I.
FPPFLAGS=
LOCDIR=
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/test

topopt: main.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o chkopts
	rm -rf topopt
	-${CLINKER} -o topopt main.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o ${PETSC_SYS_LIB}
	${RM}  main.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o 
	rm -rf *.o

myclean:
	rm -rf topopt *.o output* binary* log* makevtu.pyc restore* 

