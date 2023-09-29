!----------------------------------------------------------------------------------------
! PROGRAM TO COMPUTE FAST FOURIER TRANSFROM
! 3D ARRAY
! REAL TO COMPLEX (FORWARD)
! IN OF PLACE
! MPI BASED
! ---------------------------------------------------------------------------------------
! AUTHOR: SUGAN D MURUGAN
! DATE: MAY 2020
! ---------------------------------------------------------------------------------------
program fft_mpi_3d_r2c
    
    use,intrinsic::iso_c_binding
    implicit none
!   _____________________________________
!   HEADER FILE INCLUSION FOR MPI,FFT_MPI
    include 'mpif.h'
    include 'fftw3-mpi.f03'
! ---------------------------------------------------------------------------------------

!   _____________________
!   VARIABLES DECLARATION (FFT)
    integer(C_INTPTR_T),parameter::D=128
    type(C_PTR)::plan,cdata
    complex(C_DOUBLE_COMPLEX),pointer::data_fft(:,:,:)
    integer(C_INTPTR_T)::alloc_local,local_D_3,local_D_3_offset
    integer(C_INTPTR_T)::g_1,g_2,g_3
    integer::ierr,myid,nproc
    double precision::time_start,time_end
!   _____________________
!   VARIABLES DECLARATION (DATA RELATED)    
    double precision::delta
    integer(kind=8)::x_1,x_2,x_3,k_1,k_2,k_3
    character(len=40)::file_dir,file_name
    double precision,dimension(:,:,:),allocatable::vel_r
    double complex,dimension(:,:,:),allocatable::vel_k
! ---------------------------------------------------------------------------------------

!   _______________
!   INPUT FUNCTION
    allocate(vel_r(0:D-1,0:D-1,0:D-1),vel_k(0:D-1,0:D-1,0:D-1))
    delta=ATAN(1.0D0)*8.0D0/DBLE(D)
    do x_1=0,D-1
        do x_2=0,D-1
            do x_3=0,D-1
                vel_r(x_1,x_2,x_3)=2.0D0*(DSIN(delta*(7.0D0*x_1-18.0D0*x_2+4.0D0*x_3))+DCOS(delta*13.0D0*x_2)-DSIN(delta*5.0D0*x_3))
            end do
        end do
    end do
! ---------------------------------------------------------------------------------------

!   _____________________
!   INITIALIZING THE MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call FFTW_MPI_INIT()
    time_start=MPI_Wtime() ! Noting the start-time
! ---------------------------------------------------------------------------------------

!   ____________________________________
!   EACH PROCESSOR GETS LOCAL DATA ARRAY
    alloc_local=FFTW_MPI_LOCAL_SIZE_3D(D,D,D,MPI_COMM_WORLD,local_D_3,local_D_3_offset)
    cdata=FFTW_ALLOC_COMPLEX(alloc_local)
    call c_f_pointer(cdata,data_fft,[D,D,local_D_3])
!---------------------------------------------------------------------------------------

!   ____________________________________________
!   CREATE MPI PLAN FOR INPLACE FORWARD DFT
    plan=fftw_mpi_plan_dft_3d(D,D,D,data_fft,data_fft,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)
! ---------------------------------------------------------------------------------------

!   __________________________
!   INPUT ARRAY INITIALIZATION
    do g_3=1,local_D_3
        do g_2=1,D
            do g_1=1,D
                data_fft(g_1,g_2,g_3)=DCMPLX(vel_r((g_1-1),(g_2-1),local_D_3_offset+(g_3-1)),0.0D0)
            end do
        end do
    end do
! ---------------------------------------------------------------------------------------
      
!   ________________
!   EXECUTION OF FFT    
    call fftw_mpi_execute_dft(plan,data_fft,data_fft)
! ---------------------------------------------------------------------------------------
    
!   ________________
!   SAVE THE OUTPUT      
    do g_1=1,D
        do g_2=1,D
            do g_3=1,local_D_3
                vel_k(g_1-1,g_2-1,g_3-1+local_D_3_offset)=data_fft(g_1,g_2,g_3)/DBLE(D**3.0D0)
            end do
       end do
    end do
! ---------------------------------------------------------------------------------------
   
!   _____________________________
!   DEALLLOCATE AND DESTROY PLANS
    call FFTW_DESTROY_PLAN(plan)
    call FFTW_MPI_CLEANUP()
    call FFTW_FREE(cdata)
    time_end=MPI_Wtime()
    call MPI_FINALIZE(ierr)
    print*,'FFT_MPI Execution time=',time_end-time_start
! ---------------------------------------------------------------------------------------

!   _________________________
!   WRITE OUTPUT TO FILE
!    if (myid .EQ. 0) then
!        file_dir='data/'
!        file_name='fft_mpi_3d_r2c_output.dat'
!        call system('mkdir '//trim(adjustl(file_dir)))
!        open(unit=1,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name)))   
!        do k_1=0,D-1
!            do k_2=0,D-1
!                do k_3=0,D/2
!                    write(1,10,advance='no')k_1,k_2,k_3
!                    write(1,11,advance='yes')vel_k(k_1,k_2,k_3)
!                end do
!            end do
!        end do
!        close(1)
!        10 format(3I8)
!        11 format(2f24.12)
!    end if
    
end


    
