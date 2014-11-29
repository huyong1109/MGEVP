! 1D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real,dimension(518) :: r, f
  real :: d0,d1,s,t,u
  integer :: i,it,j,k,l,ll,m,n,nl,lev, iter,nc


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
  s = (1.0/n)**2
  do i = 1, n 
    r(i) = s*r(i)
  end do 

  m = n+1  ! total point
  l =1  ! first point of temporary  mesh ( left boundary) 
  ll = n  ! end point of temporary mesh ( the point next to right bound)
  nl = n+n+k-2 ! total points of all level mesh
  
  ! initalization 
  do i = 1, nl
    f(i) = 0.0
  end do 

  nc = 0

  do  while(nc .le. 10 )! for multigrid V cycle
    nc = nc +1
    iter = 0
    do while(iter .le. 10) ! GS and then coarse mesh
      iter = iter +1
      ! Gauss-Seidel Iteration 
      j  = 0
      d0 = 1.0e-12 ! very small inital value. 
      d1 = 1.0

      ! first iterate 4 times, then keep iterating until 
      ! convergence rate is slower than u(=0.7)
      do while (((j .lt. it ) .or. ((d1/d0) .lt. u)) .and. (d1 .ge. t))
        d0 = d1
        d1 = 0.0
        j = j+1
        do i = l+1, ll
          s = 0.5*(f(i-1)+f(i+1)+r(i))
          d1 = d1 +abs(s-f(i))
          f(i) = s
        end do 
        write(*,*) 'DIFF: ', d1 
      end do  ! end GS iteration 
      if (d1 .lt. t) exit


      ! Coasrse mesh 
      if (n .ne. 2) then  
        i = ll +2
        do j = l+2, ll,2
          i = i+1
          l = l+2
          f(i) = 0.0
          r(i) = 4.0*(r(j)+f(j-1)-2*f(j)+f(j+1))
        end do 
        l = l+2
        i = i +1
        n = n/2
        write(*,*) 'NUMBER OF MESH INTERVALS : ',n
        ll = ll + n+ 1
        l =l+1
      end if 

    end do  ! end coarsing 

    ! if back to the finer mesh, then exit
    if(n .eq. m-1 ) exit

    ! finer mesh 
    i = l-3 
    do j = ll, l+1,-1
      f(i) = f(i) +f(j)
      f(i+1) = f(i+1) +0.5*(f(j)+f(j+1))
      i = i-2
    end do
    !write(*,*) 'j,l :',j,l
    n = n+n
    write(*,*) 'NUMBER OF MESH INTERVALS : ',n
    ll = l-2 
    l = i
    write(*,*) 'i,j : ',i,j
    f(i+1) = f(i+1) +0.5*(f(j)+f(j+1))

  end do


  do i = 1,m
    j = i-1
    s = j/float(n)
    write(*,*) j,s,f(i)
  end do 

end program
