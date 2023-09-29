# MAKEFILE FOR BURGERS

# DEFINE VARIABLES
# ---------------------------start-----
ob=sample2.f90
var=VAR_mod.f90
var_ob=VAR_mod.o
fft=FFT_mod.f90
fft_ob=FFT_mod.o
cc=mpif90
cc_lc= -I/usr/bin
mpi_lc=-L/usr/bin/mpif90
mpi_lb_lc=-L
lb_fftw=-L/usr/include -lfftw3 -lm
run=./ex
#----------------------------end-------


# MAKEFILE
# ---------------------------start----- 
ex:$(ob)
#	$(cc) $(cc_lc) -c $(fft) $(lb_fftw) 
#	$(cc) $(cc_lc) -c $(var) 
	$(cc) $(cc_lc) -c $(ob) 
	$(cc) $(cc_lc) $(ob) $(var_ob) $(fft_ob) $(lb_fftw) -o ex 
	$(run)
#----------------------------end-------

# CLEANING
# ---------------------------start----- 
clean:
	rm ex
	rm *.mod
	rm *.o
cl:
	rm *.mod
	rm *.o
#----------------------------end-------
