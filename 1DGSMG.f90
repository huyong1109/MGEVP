! 1D multigrid example from MUDPACK
program main 
  real,dimension(518) :: r, f
  real :: d0,d1,s,t,u
  integer :: i,it,j,k,l,ll,m,n,nl

  ! set problem 
  k = 6 
  n = 2**k
  it = 4
  t = 0.0001
  u = 0.7 

  ! input right side
  do i = 1, n
    r(i) = 1.0 
  end do 

  ! initialize 
  s = (1.0/n)**2
  do i = 1, n 
    r(i) = s*r(i)
  end do 

  m = n 
  l =1 
  ll = n 
  nl = n+n+k-2
  do i = 1, nl
    f(i) = 0.0
  end do 

  d1 = 0.0 
  j = 0

  ! Gauss-Seidel Iteration 
  do 
  d0 = d1
  d1 = 0.0
  j = j + 1 

  do i = 1, ll-1
    s = 0.5*(f(i-1)+f(i+1)+r(i))
    d1 = d1 +abs(s-f(i))
    f(i) = s
  end do 








end
