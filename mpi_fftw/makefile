# MAKEFILE FOR TEST FFT_MPI

# DEFINE VARIABLES
# ---------------------------start-----
ob=fft_mpi_3d_r2c.f90
cc=mpif90
mpi_lc=-I/usr/local/include
mpi_lb_lc=-L/usr/local/include -lfftw3_mpi -lfftw3 -lm
run=mpirun -np 4 ./a.out
#----------------------------end-------


# MAKEFILE
# ---------------------------start----- 
ex:$(ob)
	$(cc) $(ob) $(mpi_lc) $(mpi_lb_lc)
	$(run)
#----------------------------end-------

# CLEANING
# ---------------------------start----- 
clean:
	rm a.out
	rm *.mod
	rm *.o
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
