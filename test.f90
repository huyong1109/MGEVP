program main 
implicit none

  integer :: i,j
  real :: a(10)

  !do i = 1, 10 
  !  j = j+i
  !end do 

  !print *, 'Without increment : i = ', i
  !do i = 1, 10 ,1
  !  j = j+i
  !end do 

  !print *, 'With an increment : i = ', i

  !do i = 1 ,10
  !  a = (8*4**(real(i)-1) + 7)/3
  !  print *, a
  !end do 

  do i = 1 ,10
    a(i)  = i
  end do 
  write(*,*)'before', a
  call testprint(a(3:4), 2)
  write(*,*) 'after', a
end 
subroutine testprint(a, n)
implicit none
integer,intent(in) :: n 
real,intent(inout) :: a(n)

!local 
integer :: i
do i = 1, n   
  print*, a(i) 
  a(i) = 111
end do 

end subroutine 
