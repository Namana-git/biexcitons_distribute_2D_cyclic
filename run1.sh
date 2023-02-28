#!/bin/bash

# Netlib
# ------

#mpif90 -o pcad.o parallel_create_and_diagonalize.f -lscalapack -llapack -lblas

# Intel MKL
# ---------

LINK='  -L${MKLROOT}/lib/intel64 -L/home/namanav/hdf5/lib -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 '
COMP='-I${MKLROOT}/include -I/home/namanav/hdf5/include'

#mpiifort -o pcad1to1.o parallel_create_and_diagonalize1to1.f $LINK $COMP


mpiifort -c  hdf5_read_par.f90 $LINK $COMP
mpiifort  -c  mod_variables.f90 $LINK $COMP
mpiifort  -c  create_H1.f90  $LINK $COMP
mpiifort  -c  blacs_and_array.f90  $LINK $COMP
mpiifort  -c  main.f90 $LINK $COMP

mpiifort  main.o blacs_and_array.o create_H1.o hdf5_read_par.o mod_variables.o $LINK $COMP

#time mpirun -np 80 ./pcad.o

