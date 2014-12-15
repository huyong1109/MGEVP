! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
module domain 
  implicit none
  
  integer, parameter :: nlev = 2
  integer, parameter :: mlev = 2
  integer, parameter :: nn = 4
  integer, parameter :: mm = 4
  integer, parameter :: choice = 2

  integer, parameter :: ngrid = nn**nlev+1  ! grid num on finest level 
  integer, parameter :: mgrid = mm**mlev+1  ! grid num on finest level 
  integer, parameter :: tngrid = nn*(nn**nlev -1)/3 +nlev ! total grids all level 
  integer, parameter :: tmgrid = mm*(mm**mlev -1)/3 +mlev ! total grids all level 
  real*8,dimension(tngrid,tmgrid) :: ax,bx,cc,ay,by
  real*8,dimension(tngrid -ngrid +2,tmgrid-mgrid+2,mm-1,mm-1) :: rinv
  real*8,parameter :: rtol = 1.0e-5
  

  private
  save

  public :: grid_num


contains 
  
subroutine grid_num(lev,nm,lgrid,tgrid,lowgrid)
  implicit none
  integer,intent(in)  :: lev
  integer,intent(inout) :: lgrid,tgrid,lowgrid ! lev grid 
                                             ! total grids from lev 1
    
  lgrid = nm**lev + 1
  tgrid  = nm*(nm**lev -1)/3 +lev
  lowgrid = nm**(lev-1) + 1

  end subroutine 
   
end module 
program main 
  implicit none

  ! local variable
  real*8 :: s,ntmp,mtmp
  integer :: i,j,ii,jj,kn,km,kn1,km1,hn,hm,n,tn,ln,m,tm,lm,iter
  integer :: ngrid2,mgrid2

  print *, '2D Multigrid EVP solver'





  subroutine inverse(a,c,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed 
  ! during the calculation
  !===========================================================
  implicit none 
  integer,intent(in) ::  n
  real*8,intent(inout) :: a(n,n)
  real*8,intent(inout) :: c(n,n)
  real*8 :: L(n,n), U(n,n), b(n), d(n), x(n)
  real*8 :: coeff
  integer::  i, j, k

  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 aloows such operations on matrices
  L=0.0
  U=0.0
  b=0.0

  ! step 1: forward elimination
  do k=1, n-1
    do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
        a(i,j) = a(i,j)-coeff*a(k,j)
      end do
    end do
  end do

  ! Step 2: prepare L and U matrices 
  ! L matrix is a matrix of the
  ! elimination coefficient
  ! + the diagonal elements are 1.0
  do i=1,n
    L(i,i) = 1.0
  end do
  ! U matrix is the upper triangular
  ! part of A
  do j=1,n
    do i=1,j
      U(i,j) = a(i,j)
    end do
  end do

  ! Step 3: compute
  ! columns of the inverse
  ! matrix C
  do k=1,n
    b(k)=1.0
    d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
      d(i)=b(i)
      do j=1,i-1
        d(i)=d(i)-L(i,j)*d(j)
      end do
    end do
    ! Step 3b:Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
    end do
    ! Step  3c: fill the solutions  x(n) into column  k  of  C
    do  i=1,n
      c(i,k) = x(i)
    end do

    b(k)=0.0
  end do 
  end subroutine





