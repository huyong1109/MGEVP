program main 
implicit none

  integer :: i,j
  real :: a

  !do i = 1, 10 
  !  j = j+i
  !end do 

  !print *, 'Without increment : i = ', i
  !do i = 1, 10 ,1
  !  j = j+i
  !end do 

  !print *, 'With an increment : i = ', i

  do i = 1 ,10
    a = (8*4**(real(i)-1) + 7)/3
    print *, a
  end do 

end 
