! PROGRAM THAT DOES FFT FOR A GIVEN 3D ARRAY FROM REAL TO COMPLEX
! ----------------------------------------------------------------------
module input_file
    implicit none
    integer(kind=4)::input_code
    double precision,parameter::two_pi=ATAN(1.0D0)*8.0D0
    double precision::data_temp,dx,junk
    contains
    subroutine input_function(i,j,k)
        integer(kind=4),intent(in)::i,j,k
        data_temp=2.0D0*(DSIN(dx*(7.0D0*(i-1)-18.0D0*(j-1)+4.0D0*(k-1)))+DCOS(dx*13.0D0*(j-1))-DSIN(dx*5.0D0*(k-1)))
    end
end module 

module output_file
    implicit none
    contains
    subroutine output_file_write(N,data_write)
        integer(kind=4),intent(in)::N
        double complex,dimension((N/2+1),N,N),intent(in)::data_write
        integer(kind=8)::i,j,k
        character(len=40)::file_dir,file_name
        file_dir='data/'
        file_name='fft_3d_r2c.dat'
        open(unit=1,file=trim(adjustl(file_dir))//trim(adjustl(file_name)))
        do i=1,(N/2+1)
        do j=1,N
        do k=1,N
             write(1,10,advance='no')i-1,j-1,k-1
             write(1,11,advance='yes')data_write(i,j,k)/(N**3.0D0)
        end do
        end do 
        end do
        close(1)
        10 format(3I8)
        11 format(2f24.12)
    end      
end module

program fft_3d_r2c

    ! HEADER FILES/MODULES INCLUSION
    ! ----------------------
    use,intrinsic::iso_c_binding ! Standard module which defines the equivalent of C types in fortran
    use input_file ! modules for input
    use output_file ! modules for output
    implicit none
    include 'fftw3.f03' ! Fortran interface files for all of the C routines for FFTW operation

    ! C VARIABLES DECLARATION
    ! ---------------------
    integer(C_INT),parameter::N=1024
    type(C_PTR)::plan,cdata_in,cdata_out ! all fftw plans are of this datatype in FORTRAN
    complex(C_DOUBLE_COMPLEX),pointer::data_out(:,:,:)
    real(C_DOUBLE),pointer::data_in(:,:,:)
    integer(C_INT)::i,j,k
    complex(C_DOUBLE_COMPLEX)::fout
    double precision::time_start,time_end
        
    ! ALLOCATE ARRAYS - DYNAMIC
    ! NOTE:- The array dimensions are in reverse order for FORTRAN
    ! ---------------------------------------
    call CPU_TIME(time_start)
    cdata_in=fftw_alloc_real(int(N*N*N,C_SIZE_T))
    call c_f_pointer(cdata_in,data_in,[N,N,N])
    cdata_out=fftw_alloc_complex(int((N/2+1)*N*N,C_SIZE_T))
    call c_f_pointer(cdata_out,data_out,[(N/2+1),N,N])
    
    ! PLAN FOR OUT-PLACE FORWARD DFT R2C
    ! -----------------------------------
    plan=fftw_plan_dft_r2c_3d(N,N,N,data_in,data_out,FFTW_ESTIMATE)

    ! INITIALIZE INPUT DATA
    ! ---------------------
    input_code=1
    if (input_code==1) then
        dx=two_pi/DFLOAT(N)
        do k=1,N
        do j=1,N
        do i=1,N
            call input_function(i,j,k)    
            data_in(i,j,k)=data_temp
        end do
        end do
        end do   
    end if
    ! EXECUTE DFT
    ! -----------
    call fftw_execute_dft_r2c(plan,data_in,data_out)
      
    ! WRITE OUTPUT TO FILE
    ! --------------------
   call CPU_TIME(time_end)
   print*,'FFT Execution time=',time_end-time_start
!    call output_file_write(N,data_out)
   
    ! DESTROY PLANS
    ! -------------
    call fftw_destroy_plan(plan)
    call fftw_free(cdata_in)
    call fftw_free(cdata_out)
 end

