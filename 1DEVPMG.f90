! 1D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*8,parameter :: rtol = 1.0e-9
  integer*8, parameter :: nlev = 3
  integer*8, parameter :: choice = 2
  integer*8, parameter :: ngrid = 4**nlev+1  ! grid num on finest level 
  integer*8, parameter :: tgrid = 4*(4**nlev -1)/3 +nlev ! total grids all level 
  real*8  :: tol 
  real*8,dimension(tgrid) :: r, f,x,xt
  real*8,dimension(tgrid) :: ax,cx,bx
  real*8,dimension(tgrid -ngrid +2) :: rinv


  ! local variable
  real*8 :: s,tmp
  integer*8 :: i,j,k,n,tn,ln
  integer*8 :: ngrid2

  tol = rtol/dble(ngrid)
  print *, '1D Multigrid EVP solver'

  ! set problem 
  j = 1
  k =1
  do i = nlev,1,-1
    call grid_num(i,n,tn,ln)
    ngrid2 = 2*(i)
    ngrid2 = 4**ngrid2
    tmp = dble(ngrid2)
    ax(j+1:j+n-2) = tmp
    bx(j+1:j+n-2) = tmp
    cx(j+1:j+n-2) = -2.0*tmp

    write(*,'(A,I4.1,a,i11.1,a,i12.1,a,i10.1,a,i10.1)') 'lev',i, '  n',n,'  tn',tn,'  axs', j, '  rinvs', k
    call lev_pre(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),rinv(k:k+ln-2),n,ln)

    !write(*,*) ax(j:j+n-1)
    !write(*,*) cx(j:j+n-1)
    !write(*,*) bx(j:j+n-1)
    !write(*,*) rinv(k:k+ln-1)
    j = j+n 
    k = k +ln
  end do 

  call analytic_init(f,ngrid,choice)

  
  ! initalization 
  do i = 1, tgrid
    r(i) = 0.0
    x(i) = 0.0
    xt(i) = 0.0
  end do 

  ! traditional iteration 
  call grid_num(nlev,n,tn,ln)
  call MG(ax,cx,bx,rinv,x,f,nlev,n,tn,ln,tol)

  !write(*,'(4f11.4)') x(:) 
  call analytic_check(x(1:n),ngrid,choice)



end program

subroutine analytic_init(f,n,choice)
  integer*8,intent(in):: n
  integer,intent(in):: choice
  real*8, dimension(n), intent(inout) :: f
  
  ! local 
  real*8 :: s,h,h2
  real*8,parameter :: pi = 4.*atan(1.) 

  ! input right side
  f(:) = 0.0
  h = 1.0/dble(n-1)
  h2= h**2

  select case(choice)
  case (1)
    do i = 2, n-1
      f(i) = -1.0
    end do 
  case (2)
    do i = 2, n-1
      s = (i-1)*h
      f(i) = 2.0-6.0*s
    end do 
  case (3)
  kn = 1
  ln =1
    do i = 2, n-1
      s = (i-1)*h
      f(i) = -4.0*pi**2*sin(2.0*pi*s)
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select

  !write(*,*) 'init f(:)'
  !write(*,*) f(:)
  
end subroutine
subroutine analytic_check(x,n,choice)
  integer*8,intent(in):: n
  integer,intent(in):: choice
  real*8, dimension(n), intent(in) :: x
  
  ! local 
  real*8 :: h,s,e,maxe
  real*8,parameter :: pi = 4.*atan(1.) 

  
  maxe = 0.0
  h = 1.0/dble(n-1)
  select case(choice)
  case (1)
  !=========== for u'' = -1=========
    h = 1.0/dble(n-1)
    do i = 1, n
      s = (i-1)*h
      e = x(i) -0.5*s*(1.0-s)
      maxe= max(abs(e),maxe) 
      write(*,*) i, s, x(i), e
    end do 

  case(2)
  !=========== for u'' = -1=========
    h = 1.0/dble(n-1)
    do i = 1, n
      s = (i-1)*h
      e = x(i) -s*s*(1.0-s)
      maxe= max(abs(e),maxe) 
      write(*,'(I6.2,3f11.4)') i, s, x(i), e
    end do 

  case(3)
  !=========== for u'' = -1=========
    do i = 1, n
      s = (i-1)*h
      e = x(i) -sin(2.0*pi*s)
      maxe= max(abs(e),maxe) 
      write(*,'(I6.2,3f11.4)') i, s, x(i), e
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select
      write(*,*) 'MAX ERROR(relative to analytic): ',maxe

end subroutine


recursive subroutine MG(ax,cx,bx,rinv,x,r,lev,n,tn,ln,tol)
  ! traditional iteration 
  integer*8, intent(in) :: lev
  integer*8, intent(in) :: n,tn,ln
  !real*8, intent(inout) :: rr
  real*8,dimension(tn),intent(in) :: ax,cx,bx
  real*8,dimension(tn -n +2),intent(in) :: rinv
  real*8,dimension(tn),intent(inout) :: r,x
  real*8,intent(in) :: tol

  integer*8 :: i,j,k,j1,k1,j1n,k1n,n1,tn1,ln1
  integer :: iter
  real*8 :: rr
  real*8,dimension(tn) :: f,xt
  
  j=1
  k =1

  xt(:) = 0.
  write(*,*) 'r1'
  write(*,'(4f11.2)') r(j:j+n-1) 
  call lev_evp(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),rinv(k:k+ln-2),xt(j:j+n-1),r(j:j+n-1),n,ln)
  f(j:j+n-1) = r(j:j+n-1)
  x(j:j+n-1) = x(j:j+n-1)+xt(j:j+n-1)
  call calc_rhs(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),xt(j:j+n-1),f(j:j+n-1),r(j:j+n-1),n,rr)
  write(*,*) 'r2'
  write(*,'(4f11.2)') r(j:j+n-1) 
 
   write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV',lev, 'rr= ', rr, 'tol= ',tol
  !write(*,*) 'r'
  !write(*,'(17f7.4)') r(j:j+n-1) 
  !write(*,*) 'x'
  !write(*,'(17f7.4)') x(j:j+n-1) 
   
  if (lev .ne. 1) then
    call grid_num(lev-1,n1,tn1,ln1)
    j1 = j +n
    j1n = j1 +tn1
    k1 = k +ln
    k1n = k1 +tn1-n1+2
    iter  = 0
    !do while ((rr .gt. tol) .and. (iter < 5) )
      !iter = iter +1
      ! coarsing 
  !write(*,'(17f7.2)') r(j:j+n-1) 
      call coarse_rhs(r(j:j+n+ln-1), n,ln)
  !write(*,'(5f7.2)') r(j+n:j+n+ln-1) 
      xt(:) = 0.
      call MG(ax(j1:j1n),cx(j1:j1n),bx(j1:j1n),rinv(k1:k1n),xt(j1:j1n),r(j1:j1n),lev-1,n1,tn1,ln1,tol)
      ! relaxing
      call finer_x(xt(1:n+ln),n,ln) 
      call lev_evp(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),rinv(k:k+ln-2),xt(j:j+n-1),r(j:j+n-1),n,ln)
      f(j:j+n-1) = r(j:j+n-1)
      x(j:j+n-1) = x(j:j+n-1)+xt(j:j+n-1)
      call calc_rhs(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),xt(j:j+n-1),f(j:j+n-1),r(j:j+n-1),n,rr)
      write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV',lev, 'rr= ', rr, 'tol= ',tol
      !write(*,*) 'r'
      !write(*,'(17f7.4)') r(j:j+n-1) 
      !write(*,*) 'x'
      !write(*,'(17f7.4)') x(j:j+n-1) 

    !end do 
  end if 
end subroutine

subroutine coarse_rhs(r,n,ln)
  implicit none
  integer*8,intent(in) :: n,ln
  real*8,dimension(n+ln),intent(inout) :: r

  !local 
  integer*8 :: i,j
  j = 2
  r(n+1:n+ln) = 0.0
  do i = 5, n-1, 4
    r(n+j) = 0.25*r(i)
    j = j+1
  end do
end subroutine 

subroutine finer_x(x,n,ln)
  implicit none
  integer*8,intent(in) :: n,ln
  real*8,dimension(n+ln),intent(inout) :: x

  !local 
  integer*8 :: i,j
  j = 2
  do i = 5, n-1, 4
    x(i) = x(i) +x(n+j)
    j = j+1
  end do
end subroutine 

subroutine calc_rhs(ax,cx,bx,x,f,r,n,rr)
  implicit none
  integer*8,intent(in) :: n
  real*8, intent(inout) :: rr
  real*8,dimension(n),intent(in) :: ax,cx,bx
  real*8,dimension(n),intent(in) :: x,f
  real*8,dimension(n),intent(inout) :: r

  !local 
  integer*8 :: i
  rr = 0.0
  do i = 2, n-1
    r(i) = f(i) -ax(i)*x(i-1)-cx(i)*x(i)-bx(i)*x(i+1)
    rr = rr + abs(r(i)**2)
  end do
  rr = sqrt(rr/dble(n))
end subroutine 

subroutine lev_pre(ax,cx,bx,rinv,n,ln)
  implicit none
  integer*8,intent(in) :: n,ln! grid number of current and lower (coarser) levels
  real*8,dimension(n),intent(in) :: ax,cx,bx
  real*8,dimension(ln-1),intent(inout) :: rinv

  !local 
  integer*8 :: i,j
  j = 0
  do i = 1, n-1, 4
    j = j+1
    call pre(ax(i:i+4), cx(i:i+4),bx(i:i+4),rinv(j))
  end do

end subroutine 

! solve each block in one level
subroutine lev_evp(ax,cx,bx,rinv,x,r,n,ln)
  implicit none
  integer*8,intent(in) :: n,ln
  real*8,dimension(n),intent(in) :: ax,cx,bx
  real*8,dimension(n),intent(in) :: r
  real*8,dimension(ln-1),intent(in) :: rinv
  real*8,dimension(n),intent(inout) :: x

  !local 
  integer*8 :: i,j
  j = 0
  do i = 1, n-1, 4
    j = j+1
   call evp(ax(i:i+4), cx(i:i+4),bx(i:i+4),rinv(j),x(i:i+4),r(i:i+4))

  end do
end subroutine 

! 5 points EVP 
subroutine evp(ax,cx,bx,rinv,x,r)
  implicit none
  integer,parameter :: n = 5
  real*8,dimension(n),intent(in) :: ax,cx,bx
  real*8,dimension(n),intent(in) :: r
  real*8,dimension(n),intent(inout) :: x

  real*8, intent(in) :: rinv

  !local 
  integer :: i
  real*8,dimension(n) :: y

  y(:) = x(:) 
  do i = 2, n-1
    y(i+1) = (r(i) -ax(i)*y(i-1)-cx(i)*y(i))/bx(i)
  end do

  x(2) = x(2) + (y(n)-x(n))*rinv
  do i = 2, n-2
    x(i+1) = (r(i) -ax(i)*x(i-1)-cx(i)*x(i))/bx(i)
  end do
end subroutine 

! 5 points PRE
subroutine pre(ax,cx,bx,rinv)
  implicit none
  integer,parameter :: n = 5
  real*8,dimension(n),intent(in) :: ax,cx,bx

  real*8, intent(inout) :: rinv

  !local 
  integer :: i
  real*8,dimension(n) :: y

  y(:) = 0.0
  y(2) = 1.0
  do i = 2, n-1
    y(i+1) = ( -ax(i)*y(i-1)-cx(i)*y(i))/bx(i)
  end do

  rinv = -1.0/y(5)
end subroutine 

subroutine grid_num(lev, lgrid,tgrid,lowgrid)
  implicit none
  integer*8,intent(in)  :: lev
  integer*8,intent(inout) :: lgrid,tgrid,lowgrid ! lev grid 
                                             ! total grids from lev 1
                                             !to level lev
  lgrid = 4**lev + 1
  tgrid  = 4*(4**lev -1)/3 +lev
  lowgrid = 4**(lev-1) + 1

end subroutine 
