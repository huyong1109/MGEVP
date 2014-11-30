! 1D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real,parameter :: rtol = 0.0001
  integer, parameter :: nlev = 2
  integer, parameter :: ngrid = (2*4**nlev+7)/3  ! grid num on finest level 
  integer, parameter :: tgrid = (8*(4**nlev -1) +21*nlev)/9 ! total grids all level 
  !integer, parameter :: tingrid = (8*(4**(nlev-1) -1) &     ! total inner grid 
  !    +21*(nlev-1))/9  +1-2*(nlev- 1)
  real,dimension(tgrid) :: r, f
  real,dimension(tgrid) :: ax,cx,bx
  real,dimension(tgrid -ngrid ) :: rinv

  ! function 
  !integer :: grid_num

  ! local variable
  real :: s,tmp
  integer :: i,j,k,n,tn,ln


  print *, '1D Multigrid EVP solver'

  ! set problem 
  j = 1
  k =1
  do i = nlev,1,-1
    call grid_num(i,n,tn,ln)
    tmp = 1.0/real(4**(2*(nlev-i)))
    ax(j+1:j+n-2) = tmp
    bx(j+1:j+n-2) = tmp
    cx(j+1:j+n-2) = -2.0*tmp
    if ( i .lt. nlev ) then 
    ax(j+1) = 2.0*tmp
    cx(j+1) = -3.0*tmp
    bx(j+n-2) = 2.0*tmp
    cx(j+n-2) = -3.0*tmp
    end if 

    call lev_rep(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),rinv(k+1:k+ln-2),n,ln)
    j = j+n 
    k = k +ln
    print *, 'lev : ',i, ' ax_end : ', j, ' rinv_end :', k
  end do 

  ! input right side
  do i = 1, ngrid
    r(i) = 1.0 
  end do 
  s = (1.0/(n-1))**2
  do i = 1, ngrid 
    r(i) = s*r(i)
  end do 

  
  ! initalization 
  do i = 1, tgrid
    f(i) = 0.0
  end do 


end program

subroutine calc_rhs(ax,cx,bx,x,f,r,n)
  implicit none
  integer,intent(in) :: n
  real,dimension(n),intent(in) :: ax,cx,bx
  real,dimension(n),intent(in) :: x,f
  real,dimension(n),intent(inout) :: r

  !local 
  integer :: i
  do i = 2, n-1
    r(i) = f(i) -ax(i)*x(i-1)-cx(i)*x(i)-bx(i)*x(i+1)
  end do
end subroutine 

subroutine lev_rep(ax,cx,bx,rinv,n,ln)
  implicit none
  integer,intent(in) :: n,ln! grid number of current and lower (coarser) levels
  real,dimension(n),intent(in) :: ax,cx,bx
  real,dimension(ln),intent(in) :: rinv

  !local 
  integer :: i,j
  j = 0
  do i = 1, n, 4
    j = j+1
    call rep(ax(i:i+4), cx(i:i+4),bx(i:i+4),rinv(j))
  end do
end subroutine 

! solve each block in one level
subroutine lev_evp(ax,cx,bx,x,r,n)
  implicit none
  integer,intent(in) :: n
  real,dimension(n),intent(in) :: ax,cx,bx
  real,dimension(n),intent(in) :: r
  real,dimension(n),intent(inout) :: x

  !local 
  integer :: i
end subroutine 

! 5 points REP
subroutine rep(ax,cx,bx,rinv,r)
  implicit none
  integer,parameter :: n = 5
  real,dimension(n),intent(in) :: ax,cx,bx,r

  real, intent(inout) :: rinv

  !local 
  integer :: i
  real,dimension(n) :: y

  y(:) = 0.0
  y(2) = 1.0
  do i = 2, n-1
    y(i+1) = (r(i) -ax(i)*y(i-1)-cx(i)*y(i))/bx(i)
  end do

  rinv = -1.0/y(5)
end subroutine 
! 5 points EVP 
subroutine evp(ax,cx,bx,rinv,x,r)
  implicit none
  integer,parameter :: n = 5
  real,dimension(n),intent(in) :: ax,cx,bx
  real,dimension(n),intent(in) :: r
  real,dimension(n),intent(inout) :: x

  real, intent(in) :: rinv

  !local 
  integer :: i
  real,dimension(n) :: y

  y(:) = x(:) 
  do i = 2, n-1
    y(i+1) = (r(i) -ax(i)*y(i-1)-cx(i)*y(i))/bx(i)
  end do

  x(2) = x(2) + (y(n)-x(n))/rinv
  do i = 2, n-1
    x(i+1) = (r(i) -ax(i)*x(i-1)-cx(i)*x(i))/bx(i)
  end do
end subroutine 

subroutine grid_num(lev, lgrid,tgrid,lowgrid)
  implicit none
  integer,intent(in)  :: lev
  integer,intent(inout) :: lgrid,tgrid,lowgrid ! lev grid 
                                             ! total grids from lev 1
                                             !to level lev
    
  lgrid = 2*4**lev + 7
  tgrid  = (8*(4**lev -1) +21*lev)/9
  if (lev .gt. 1) then 
    lowgrid = (2*4**(lev-1) + 7)/3
  else if (lev .eq. 1 ) then 
    lowgrid = 1
  else 
    print *, '!!!Error: LEV should be positive integer '
    stop
  end if 

  if (mod(lgrid, 3) .ne. 0) then 
    print *, 'Error in grid number',lgrid
  end if 
  lgrid = lgrid/3

  end subroutine 





