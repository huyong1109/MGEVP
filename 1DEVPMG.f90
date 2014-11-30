! 1D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real,parameter :: rtol = 0.0001
  integer, parameter :: nlev = 2
  integer, parameter :: ngrid = (2*4**nlev+7)/3
  integer, parameter :: tgrid = (8*(4**nlev -1) +21*nlev)/9
  real,dimension(tgrid) :: r, f
  real,dimension(tgrid) :: ax,cx,bx

  ! function 
  !integer :: grid_num

  ! local variable
  real :: s
  integer :: i,j,k,n,tn



  ! set problem 
  !do i = 1 ,nlev
  !  call grid_num(i,n,tn)
  !  print *, n,tn
  !end do 

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

subroutine grid_num(lev, lgrid,tgrid)
  implicit none
  integer,intent(in)  :: lev
  integer,intent(out) :: lgrid,tgrid
    
  lgrid = 2*4**lev + 7
  tgrid  = (8*(4**lev -1) +21*lev)/9
  if (mod(lgrid, 3) .ne. 0) then 
    print *, 'Error in grid number',lgrid
  end if 
  lgrid = lgrid/3

  end subroutine 





