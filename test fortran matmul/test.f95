! cd "D:\_install\_Programming\parallel_computing\_matrix_product_tests\test fortran matmul"

! gfortran test.f95 -o test -O3 -ffast-math -lblas -llapack
! Time of operation (load A,B,C2) was    13.703 seconds
! Matrix product in fortran (matmul): C = A * B for the size (        2000 x        2000 )
! Time of operation was    1.02775 seconds
! Matrix product in fortran (dgemm): C = A * B for the size (        2000 x        2000 )
! Time of operation was  0.3133 seconds
! Norm c-cdgemm=    1.5695E-10 Time of operation was    1.6E-2  seconds
! Norm c-c2=    1.5695E-10 Time of operation was    3.1E-2  seconds


! gfortran  -g -fcheck=all -Wall -Wno-tabs test.f95 -o test -O3 -ffast-math -lblas -llapack

! gfortran  -g -fcheck=all -Wall -Wno-tabs test.f95 -o test -O3 -fexternal-blas "D:\_install\_Programming\parallel_computing\OpenBLAS\0.3.9 cyg\64\lib\cygopenblas_v0.3.9-gcc_9_3_0p-r0.3.9.a"
! Time of operation (load A,B,C2) was    13.655s
! Matrix product in fortran (matmul): C = A * B for the size (2000x2000)
! openblas_get_num_threads =            0
! openblas_get_num_threads =            0
! Time of operation was   0.2492seconds
! Matrix product in fortran (dgemm): C = A * B for the size (2000x2000 )
! Time of operation was   0.2479seconds
! Norm c-cdgemm=    0 Time of operation was    7E-03  seconds
! Norm c-c2=    0 Time of operation was    6E-03  seconds


! gfortran test.f95 -o test -O3 -ffast-math -fopenmp -fexternal-blas "D:\_install\_Programming\parallel_computing\OpenBLAS\0.3.9 openmpcyg\64\lib\cygopenblas_v0.3.9-gcc_9_3_0p-r0.3.9.a"

!openmp conflicts with openblas (disable it)

program product_test
!	Use omp_lib
	implicit none;

	INTEGER :: N,NA,NB,NC
!	INTEGER :: i,j,k
	INTEGER :: N_threads
	INTEGER :: ncout,N_count
	LOGICAL  :: isLoad
	INTEGER :: t_begin,t_end, clock_rate
!	REAL*8 :: time_end,time_begin
	REAL*8 :: t1
	REAL*8 :: tdgemm
	REAL(8),dimension (:,:), allocatable :: a,b
	REAL(8),dimension (:,:), allocatable :: c,c2
	REAL(8),dimension (:,:), allocatable :: cdgemm
!	REAL*8 :: t2
!	REAL*8,dimension (:,:), allocatable :: c3,c4
!	INTEGER,dimension (:), allocatable:: cs
	character(*), parameter :: pathstr="D:\_install\_Programming\parallel_computing\_matrix_product_tests\"

!	REAL :: x(5) = [1, 2, 3, 4, 5 ]
!	print *, 'norm2(x)', norm2(x)  ! = sqrt(55.) ~ 7.416

	isLoad = .TRUE.
!	isLoad = .FALSE.
	N = 2000
	N_count = 20

!	call CPU_TIME ( time_begin )
	call SYSTEM_CLOCK ( t_begin, clock_rate )
	IF (isLoad) THEN
		print *, 'Loading matrix A'
		CALL LoadMatrix(pathstr//"A.dat",N,NA,A)
		print *,'a[0,0]=',a(1,1)
		print *, 'Loading matrix B'
		CALL LoadMatrix(pathstr//"B.dat",N,NB,B)
		print *,'b[0,0]=',b(1,1)
		print *, 'Loading matrix C'
		CALL LoadMatrix(pathstr//"C.dat",N,NC,C2)
		IF ((NA/=NB) .OR. (NA/=NC)) THEN
			print *, 'Error in A,B and C matrices'
		ENDIF
	ELSE
		allocate ( a(N,N),b(N,N),c(N,N) ) 
	!	allocate ( c3(N,N) ) 
	!	allocate ( c4(N,N) ) 
		call RANDOM_NUMBER(a)
		call RANDOM_NUMBER(b)
	ENDIF
!	call CPU_TIME ( time_end )
!	print *, 'Time of operation was ', time_end - time_begin, ' seconds'
	call SYSTEM_CLOCK ( t_end, clock_rate )
	print *, 'Time of operation was ', (t_end - t_begin)/real(clock_rate), ' seconds'


!print *, 'num_procs=', omp_get_num_procs()
!print *, 'max_threads=', omp_get_max_threads()
!print *, 'num_threads=', omp_get_num_threads()

	print *
	print *, 'Matrix product in fortran (matmul): C = A * B for the size (',N,'x',N,')'


!	call openblas_get_num_threads(N_threads)
!	print *, 'openblas_get_num_threads = ', N_threads
!	CALL openblas_set_num_threads(1)
!	CALL openblas_set_num_threads(4)
!	call openblas_get_num_threads(N_threads)
!	print *, 'openblas_get_num_threads = ', N_threads

!	call CPU_TIME ( time_begin )
	call SYSTEM_CLOCK ( t_begin, clock_rate )
	do ncout=1,N_count
!!$OMP PARALLEL
!print *, 'num_threads=', omp_get_num_threads()
		c = MATMUL(a,b)
!!$OMP END PARALLEL
	end do
!	call CPU_TIME ( time_end )
!	t1 = (time_end - time_begin)/N_count
	call SYSTEM_CLOCK ( t_end, clock_rate )
	t1 = (t_end - t_begin)/real(clock_rate*N_count)
	print *, 'Time of operation was ', t1, ' seconds'

!subroutine dgemm 	( 	character  	TRANSA,
!		character  	TRANSB,
!		integer  	M,
!		integer  	N,
!		integer  	K,
!		double precision  	ALPHA,
!		double precision, dimension(lda,*)  	A,
!		integer  	LDA,
!		double precision, dimension(ldb,*)  	B,
!		integer  	LDB,
!		double precision  	BETA,
!		double precision, dimension(ldc,*)  	C,
!		integer  	LDC 
!	) 
!	C := alpha*op( A )*op( B ) + beta*C

	print *
	print *, 'Matrix product in fortran (dgemm): C = A * B for the size (',N,'x',N,')'
	allocate ( cdgemm(N,N) ) 

!	call CPU_TIME ( time_begin )
	call SYSTEM_CLOCK ( t_begin, clock_rate )
	do ncout=1,N_count
!!$OMP PARALLEL
		call dgemm('N','N',N,N,N,1.d0,A,N,B,N,0.d0,Cdgemm,N)
!!$OMP END PARALLEL
	end do
!	call CPU_TIME ( time_end )
!	tdgemm = (time_end - time_begin)/N_count
	call SYSTEM_CLOCK ( t_end, clock_rate )
	tdgemm = (t_end - t_begin)/real(clock_rate*N_count)
	print *, 'Time of operation was ', tdgemm, ' seconds'

	print *
!	call CPU_TIME ( time_begin )
	call SYSTEM_CLOCK ( t_begin, clock_rate )
	print *, 'Norm c-cdgemm= ', NORM2(c-Cdgemm)
!	call CPU_TIME ( time_end )
!	print *, 'Time of operation was ', time_end - time_begin, ' seconds'
	call SYSTEM_CLOCK ( t_end, clock_rate )
	print *, 'Time of operation was ', (t_end - t_begin)/real(clock_rate), ' seconds'



!print *, 'SIZEOF(REAL(10))',SIZEOF(REAL(10))
!print *, 'SIZEOF(a(1,1))',SIZEOF(a(1,1))
!print *, 'SIZEOF(a)',SIZEOF(a)

!	allocate ( c2(N,N) ) 
!	print *, 'Matrix Multiplication: C = A * B for size (',N,',',N,')'
!	call CPU_TIME ( time_begin )
!	do ncout=1,N_count
!		do i=1,N
!			do j=1,N
!!				c2(i,j) = sum(a(i,:)*b(:,j))
!				c2(i,j) = 0
!				do k=1,N
!					c2(i,j) = c2(i,j) + a(i,k)*b(k,j)
!				end do
!			end do
!		end do
!	end do
!	call CPU_TIME ( time_end )
!	t2 = (time_end - time_begin)/N_count
!	print *, 'Time of operation was ', t2, ' seconds'
!	print *, '(ratio = ', t1/t2,')'

!
	IF (isLoad) THEN
		print *
!		call CPU_TIME ( time_begin )
		call SYSTEM_CLOCK ( t_begin, clock_rate )
		print *, 'Norm c-c2= ', NORM2(c-c2)
!		call CPU_TIME ( time_end )
!		print *, 'Time of operation was ', time_end - time_begin, ' seconds'
		call SYSTEM_CLOCK ( t_end, clock_rate )
		print *, 'Time of operation was ', (t_end - t_begin)/real(clock_rate), ' seconds'
	ENDIF




!! elemet-by-element
!	print *
!
!	print *, 'Matrix Multiplication: C3 = A * B elemet-by-element'
!	call CPU_TIME ( time_begin )
!	do i=1,N
!		do j=1,N
!			c3(i,j) = a(i,j)*b(i,j)
!		end do
!	end do
!	call CPU_TIME ( time_end )
!	print *, 'Time of operation was ', time_end - time_begin, ' seconds'
!
!	print *, 'Matrix Multiplication: C4 = A * B elemet-by-element'
!	call CPU_TIME ( time_begin )
!
!	c4 = a*b
!
!	call CPU_TIME ( time_end )
!	print *, 'Time of operation was ', time_end - time_begin, ' seconds'
!
!	print *, 'Norm c4-c3= ', NORM2(c4-c3)


contains
SUBROUTINE LoadMatrix(fname,n_rows,n_cols,A)
	INTEGER :: n_rows,n_cols,IO_Status, Allocate_Status
	INTEGER :: nLast
	CHARACTER (LEN=*) :: fname
	REAL(8),dimension (:,:), allocatable :: A
!	CHARACTER(1) :: dummyrow
	CHARACTER(100000) :: dummyrow
!	REAL(8) :: dummycol

	open(10, file=fname)

	! -- Count the rows
	n_rows = 0
	n_cols = 0
	DO
		READ( 10, '(a)', IOSTAT = IO_Status ) dummyrow
		IF ( IO_Status < 0 ) THEN
			EXIT
		ELSE IF ( IO_Status > 0 ) THEN
!			...Process error....
		END IF

		if (n_rows==0) then
			nLast = LNBLNK(dummyrow)
			n = nLast
			DO WHILE (n>2)
				nLast = LNBLNK(dummyrow(1:n))
				if (nLast<n) then
					n_cols = n_cols + 1
					n = nLast
				else
					n = n-1
				end if
			ENDDO
		endif

		n_rows = n_rows + 1
	ENDDO
	print *, 'n_rows=',n_rows
	print *, 'n_cols=',n_cols
	REWIND(10)



!	allocate ( A(n_rows,n_rows) , STAT = Allocate_Status )
	allocate ( A(n_rows,n_cols) , STAT = Allocate_Status )
	IF ( Allocate_Status /= 0 ) THEN
!	...Process error...
	END IF

	! read in values
	read(10,*) A
	A = transpose(A)
!	read(10,*) ((A(i,j), j=1,n_cols), i=1,n_rows)

	close(10)


	RETURN
END
end

