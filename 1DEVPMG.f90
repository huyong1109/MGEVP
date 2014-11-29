! 1D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real,dimension(518) :: r, f
  real :: rtol = 0.0001
  integer, parameter :: nlev = 6
  integer :: i,j,k


  ! set problem 
  do i = 1 ,nlev
    n = grid_num(k)
    print *, n
  end do 

  ! input right side
  do i = 1, n
    r(i) = 1.0 
  end do 
  s = (1.0/(n-1))**2
  do i = 1, n 
    r(i) = s*r(i)
  end do 

  
  ! initalization 
  do i = 1, n
    f(i) = 0.0
  end do 


end program



function grid_num(lev)
  integer,intent(in)  :: lev
  integer,intent(out) :: grid_num
    
  grid_num = 8*4**((lev)-1) + 7
  if (mod(grid_num, 3) .ne. 0) then 
    print *, 'Error in grid number'
  end if 
  grid_num = grid_num/3
  return 
end function 





