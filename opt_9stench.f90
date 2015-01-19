program main
  
  integer, parameter :: n =10, m=10
  real*8,dimension(n,m) :: cc,ns,ew,ne,ine
  real,dimension(n,m) :: f


  integer :: i,j,k
  integer, parameter :: iter = 10000
  real*8,dimension(n,m) :: y
  real*8,dimension(n) ::ry
  real*8 :: start, stopt


  ! initalize 
  cc(:,:) = 10.0
  ns(:,:) = -0.5
  ew(:,:) = 0.5
  ne(:,:) = -2.5
  ine(:,:) = 0.1
  f(:,:) = 99.99
  call cpu_time(start)
  do k = 1, iter
    do j = 2, m-1
      do i = 2, n-1
        y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
        -  ns(i,j)     * y(i,j+1)    &
        -  ns(i,j-1)   * y(i,j-1)    &
        -  ew(i,j)     * y(i+1,j)    &
        -  ew(i-1,j)   * y(i-1,j)    &
        -  ne(i,j-1)   * y(i+1,j-1)  &
        -  ne(i-1,j)   * y(i-1,j+1)  & 
        -  ne(i-1,j-1) * y(i-1,j-1) ) *ine(i,j) 
      end do
    end do
  end do
  call cpu_time(stopt)
  write(*,*) 'Orignl time : ', stopt -start
  write(*,*) y(:,m)

  call cpu_time(start)
  do k = 1, iter
    do j = 2, m-1
      ry(1:2) = y(1:2,j+1)
      ry(3:n) = f(2:n-1,j)                       &
      -  ne(2:n-1,j-1)   * y(3:n,j-1)    &
      -  ne(1:n-2,j-1)   * y(1:n-2,j-1)  &
      -  ew(2:n-1,j)     * y(3:n,j)      &
      -  ew(1:n-2,j)     * y(1:n-2,j)  
      do i = 2, n-1
        ry(i+1) = (ry(i+1)              &
        -  cc(i,j)     * y(i,j)    &
        -  ns(i,j-1)   * y(i,j-1)    &
        -  ns(i,j)  * ry(i)       & 
        -  ne(i-1,j)* ry(i-1))*ine(i,j) 
      end do
      y(:,j+1) = ry(:)
    end do
  end do
        
  call cpu_time(stopt)
  write(*,*) 'Optimt time : ', stopt -start

  write(*,*) y(:,m)

  call cpu_time(start)
  do k = 1, iter
    do j = 2, m-1
      ry(1:2) = y(1:2,j+1)
      ry(3:n) = f(2:n-1,j)                       &
      -  ne(2:n-1,j-1)   * y(3:n,j-1)    &
      -  ne(1:n-2,j-1)   * y(1:n-2,j-1)  &
      -  ew(2:n-1,j)     * y(3:n,j)      &
      -  ew(1:n-2,j)     * y(1:n-2,j)    &
      -  cc(2:n-1,j)     * y(2:n-1,j)    &
      -  ns(2:n-1,j-1)   * y(2:n-1,j-1)    
      do i = 2, n-1
        ry(i+1) = (ry(i+1)              &
        -  ns(i,j)  * ry(i)       & 
        -  ne(i-1,j)* ry(i-1))*ine(i,j) 
      end do
      y(:,j+1) = ry(:)
    end do
  end do
        
  call cpu_time(stopt)
  write(*,*) 'Optimt time : ', stopt -start
  write(*,*) y(:,m)
end program
