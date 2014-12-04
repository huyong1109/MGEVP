program main 
implicit none

  integer :: i,j
  real :: a(5,5)

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
  do j = 1 ,10
    a(i,j)  = 1
  end do 
  end do 
  write(*,*)'before', a
  call testprint(a(2:4,2:4), 3)
  write(*,*) 'after', a
end 
subroutine testprint(a, n)
implicit none
integer,intent(in) :: n 
real,intent(inout) :: a(n,n)

!local 
integer :: i,j
do i = 1, n   
do j = 1, n   
  print*, a(i,j) 
  a(i,j) = 2
end do 
end do 

end subroutine 
