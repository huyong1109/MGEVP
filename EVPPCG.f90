! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*8,parameter :: rtol = 1.0e-8
  integer, parameter :: nb = 16 ! x direct
  integer, parameter :: mb = 16 ! y direct
  integer, parameter :: nn = 4 ! block size 
  integer, parameter :: choice = 2
  integer, parameter :: solver = 1
  integer, parameter :: n = nb*(nn-1)+1
  integer, parameter :: m = mb*(nn-1)+1
  !integer, parameter :: tingrid = (8*(4**(nlev-1) -1) &     ! total inner grid 
  !    +21*(nlev-1))/9  +1-2*(nlev- 1)
  real*8  :: tol ,rr,mtmp,ntmp
  real*8,dimension(n,m) :: r,f,u,ut
  real*8,dimension(n,m) :: ax,bx,cc,ay,by
  real*8,dimension(nb,mb,nn-2,nn-2) :: rinv

  ! function 
  !integer :: grid_num

  ! local variable
  integer :: ii,jj,iter

  tol = rtol
  print *, '2D Multigrid EVP solver'

  ! set problem 
   
    ntmp = dble(n -1)**2
    mtmp = dble(m -1)**2
    do ii = 1, n
      do jj = 1,m
        ax(ii,jj) = ntmp
        bx(ii,jj) = ntmp
        ay(ii,jj) = mtmp
        by(ii,jj) = mtmp
        cc(ii,jj) = -ax(ii,jj)-bx(ii,jj)-ay(ii,jj)-by(ii,jj)
      end do 
    end do 

  call lev_pre(ax,bx,cc,ay,by,rinv(:,:,:,:),n,nb,m,mb,nn)

  call analytic_init(f,n,m,choice)


  ! initalization 
  do ii = 1, n
    do jj = 1, m
      r(ii,jj) = 0.0
      u(ii,jj) = 0.0
      ut(ii,jj) = 0.0
    end do 
  end do 

  ! traditional iteration 
  write(*,*) 'f'
  write(*,'(17f7.3)') f(:,:)
  rr = 1.0
  iter = 0
  if (solver == 1) then 
   ! do while ((rr  > tol) .and. (iter <8))
   !   iter =iter +1
   !   call lev_evp(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn)
   !   call calc_rhs(ax(1:n,1:m),bx(1:n,1:m),cc(1:n,1:m),ay(1:n,1:m),by(1:n,1:m),u,f,r,n,m,rr)
   !   write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'CG1 iter',iter, 'rr= ', rr, 'tol= ',tol
   !   write(*,*) 'r'
   !   write(*,'(17f7.3)') r(:,:)

   !   !call pcg(ax,bx,cc,ay,by,u,f,n,m)
   !   call evppcg(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn,tol)
   !   call calc_rhs(ax(1:n,1:m),bx(1:n,1:m),cc(1:n,1:m),ay(1:n,1:m),by(1:n,1:m),u,f,r,n,m,rr)
   !   write(*,*) 'r'
   !   write(*,'(17f7.3)') r(:,:)

   !   write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'CG2 iter',iter, 'rr= ', rr, 'tol= ',tol
   ! end do 
    call evppcg(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn,tol)
    write(*,*) 'evppcg end'
  else if (solver == 2 ) then

    call pcg(ax,bx,cc,ay,by,u,f,n,m,tol)
    write(*,*) 'pcg end'
  else if (solver == 3 ) then
  
    call diagpcg(ax,bx,cc,ay,by,u,f,n,m,tol)
    write(*,*) 'diagpcg end'
  end if 
  !call evppcg(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn)
  call analytic_check(u,n,m,choice)



  end program

  subroutine analytic_init(f,n,m,choice)
  implicit none
  integer,intent(in):: n,m
  integer,intent(in):: choice
  real*8, dimension(n,m), intent(inout) :: f

  ! local 
  integer :: i,j
  real*8 :: x,y,hm,hn
  real*8,parameter :: pi = 4.*atan(1.) 

  ! input right side
  f(:,:) = 0.0

  select case(choice)
  case (1)
  hn = 0.01/dble(n-1)
  hm = 1.0/dble(m-1)
    do i = 2, n-1
        x = dble(i-1)*hn
      do j = 2, m-1
        y = dble(j-1)*hm
        f(i,j) = -10000*(x*(0.01-x)+y*(1.0-y))
      end do 
    end do 
  case (2)
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)
    do i = 2, n-1
        x = dble(i-1)*hn
      do j = 2, m-1
        y = dble(j-1)*hm
        f(i,j) = -(x*(1.0-x)+y*(1.0-y))
      end do 
    end do 
  case (3)
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)

    do i = 2, n-1
        x = dble(i-1)*hn
      do j = 2, m-1
        y = dble(j-1)*hm
        f(i,j) = -8.0*pi**2*sin(2*pi*x)*sin(2*pi*y)
      end do 
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select


  end subroutine
  subroutine analytic_check(u,n,m,choice)
  implicit none 
  integer,intent(in):: n,m
  integer,intent(in):: choice
  real*8, dimension(n,m), intent(in) :: u

  ! local 
  integer :: i,j
  real*8 :: hm,hn,x,y,tu,e,maxe
  real*8,parameter :: pi = 4.*atan(1.) 


  maxe = 0.0
  select case(choice)
  case (1)
    !=========== for u'' = -1=========
  hn = 0.01/dble(n-1)
  hm = 1.0/dble(m-1)
    do i = 1, n
        x = dble(i-1)*hn
      do j = 1, m
        y = dble(j-1)*hm
        tu = 10000*0.5*x*(0.01-x)*y*(1.0-y)
        e = u(i,j) -tu
        maxe= max(abs(e),maxe) 
        write(*,'(I4.1,I4.1,E13.5,E13.5,E13.5)') i, j, u(i,j), tu, e
      end do 
    end do 

  case (2)
    !=========== for u'' = -1=========
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)
    do i = 1, n
        x = dble(i-1)*hn
      do j = 1, m
        y = dble(j-1)*hm
        tu = 0.5*x*(1.0-x)*y*(1.0-y)
        e = u(i,j) -tu
        maxe= max(abs(e),maxe) 
        write(*,'(I4.1,I4.1,E13.5,E13.5,E13.5)') i, j, u(i,j), tu, e
      end do 
    end do 

  case (3)
    !=========== for u'' = -1=========
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)
    do i = 1, n
        x = dble(i-1)*hn
      do j = 1, m
        y = dble(j-1)*hm
        tu = sin(2.0*pi*x)*sin(2.0*pi*y)
        e = u(i,j) -tu
        maxe= max(abs(e),maxe) 
        write(*,'(I4.1,I4.1,E13.5,E13.5,E13.5)') i, j, u(i,j), tu, e
      end do 
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select
  write(*,*) 'MAX ERROR(relative to analytic): ',maxe

  end subroutine



  subroutine calc_rhs(ax,bx,cc,ay,by,u,f,r,n,m,rr)
  implicit none
  integer,intent(in) :: n,m
  real*8, intent(inout) :: rr
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: u,f
  real*8,dimension(n,m),intent(inout) :: r
  real*8,dimension(n,m) :: rx,ry

  !local 
  integer :: i,j
  rr = 0.0
  rx(:,:) = 0.
  ry(:,:) = 0.
  do i = 2, n-1
    do j = 2, m-1
      rx(i,j)=ax(i,j)*u(i-1,j)+bx(i+1,j)*u(i+1,j)+0.5*cc(i,j)*u(i,j) 
      ry(i,j)=0.5*cc(i,j)*u(i,j)+ay(i,j)*u(i,j-1)+by(i,j)*u(i,j+1)
      r(i,j) = f(i,j) -rx(i,j)-ry(i,j)
      rr = rr + abs(r(i,j)**2)
    end do
  end do
  !rr = sqrt(rr/dble(n*m))
  rr = rr/dble(n*m)
  end subroutine 

! 5 points PRE
subroutine pre(ax,bx,cc,ay,by,rinv,nn,mm)
  implicit none
  integer :: nn,mm
  real*8,dimension(nn,mm),intent(in) :: ax,bx,cc,ay,by

  real*8,dimension(nn-2,nn-2),intent(inout) :: rinv

  !local 
  integer :: i,j,k,ii,info
  real*8,dimension(nn,mm) :: y
  real*8,dimension(nn-2,nn-2) :: work,ipiv
  real*8,dimension(nn-2,nn-2) :: rin

  y(:,:) = 0.0
  if (ax(2,2) .le. ay(2,2)) then
    do ii = 2, nn-1
      y(ii,2) = 1.0
      do j = 2, mm-1
        do i = 2, nn-1
          y(i,j+1) = ( -ax(i,j)*y(i-1,j)-bx(i,j)*y(i+1,j)-cc(i,j)*y(i,j)-ay(i,j)*y(i,j-1))/by(i,j)
        end do
      end do
      y(ii,2) = 0.0
      rinv(ii-1,:) = -y(2:nn-1,mm)
    end do 
  write(*,*) 'rin'
  write(*,'(3f11.5)') rinv(:,:)
  rin(:,:) = rinv(:,:)
  work(:,:) = rinv(:,:)
  call inverse(rin,rinv,nn-2)

  else 
    do ii = 2, mm-1
      y(2,ii) = 1.0
      do i = 2, nn-1
        do j = 2, mm-1
          y(i+1,j) = ( -ax(i,j)*y(i-1,j)-cc(i,j)*y(i,j)-ay(i,j)*y(i,j-1)-by(i,j)*y(i,j+1))/bx(i,j)
        end do
      end do
      y(2,ii) = 0.0
      rinv(ii-1,:) = -y(mm,2:mm-1)
    end do 
  write(*,*) 'rin'
  write(*,'(3f11.5)') rinv(:,:)
  rin(:,:) = rinv(:,:)
  rin(:,:) = rinv(:,:)
  work(:,:) = rinv(:,:)
  call inverse(rin,rinv,mm-2)
  endif 
  
  !! check pre rinv 
  write(*,*) 'rinv'
  write(*,'(3f11.5)') rinv(:,:)
  rin(:,:) = 0.0
  do j = 1,nn-2
    do i = 1,nn-2
      do k = 1,nn-2
        rin(i,j) = rin(i,j) + rinv(i,k)*work(k,j)
      end do 
      if (i == j ) then 
        if (abs(rin(i,j) -1.0) > 1.0e-10 ) then 
          write(*,*) 'fail in pre'
        endif 
      else
        if (abs(rin(i,j) -0.0) > 1.0e-10 ) then 
          write(*,*) 'fail in pre'
        endif 
      endif 
    end do 
  end do 
  !write(*,*) 'work'
  !write(*,*) rin(:,:)

end subroutine 

subroutine lev_pre(ax,bx,cc,ay,by,rinv,n,nb,m,mb,nn)
  implicit none
  integer,intent(in) :: n,nb,m,mb
  integer:: nn 
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(nb,mb,nn-2,nn-2),intent(inout) :: rinv

  !local 
  integer :: i,j,ie,je,is,js

  !write(*,*) 'in lev_evp_y'
  do i = 1, nb
    is = (i-1)*(nn-1)+1
    ie = i*(nn-1)+1
  do j = 1, mb
    js = (j-1)*(nn-1)+1
    je = j*(nn-1)+1
    write(*,*) 'j=',js,je,'i=',is,ie
   call pre(ax(is:ie,js:je), bx(is:ie,js:je),cc(is:ie,js:je),ay(is:ie,js:je),by(is:ie,js:je),rinv(i,j,:,:),nn,nn)

  end do
  end do

end subroutine 


! solve each block in one level
subroutine lev_evp(ax,bx,cc,ay,by,rinv,u,r,n,nb,m,mb,nn)
  implicit none
  integer,intent(in) :: n,nb,m,mb,nn
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(nb,mb,nn-2,nn-2),intent(in) :: rinv
  real*8,dimension(n,m),intent(inout) :: u

  !local 
  integer :: i,j,ie,je,is,js

  !write(*,*) 'in lev_evp_y'
  do i = 1, nb
    is = (i-1)*(nn-1)+1
    ie = i*(nn-1)+1
  do j = 1, mb
    js = (j-1)*(nn-1)+1
    je = j*(nn-1)+1
    !write(*,*) 'j=',j,'i=',i,'i+4',i+4
   call evp(ax(is:ie,js:je), bx(is:ie,js:je),cc(is:ie,js:je),ay(is:ie,js:je),by(is:ie,js:je),rinv(i,j,:,:),u(is:ie,js:je),r(is:ie,js:je),nn,nn)

  end do
  end do

end subroutine 

! 5 points EVP 
subroutine evp(ax,bx,cc,ay,by,rinv,u,r,n,m)
  implicit none
  integer:: n,m
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(n,m),intent(inout) :: u

  real*8,dimension(n-2,n-2),intent(in) :: rinv

  !local 
  integer :: i,j
  real*8,dimension(n,m) :: y,ry
  real*8 :: rr

  write(*,*) 'before evp'
  write(*,'(5f11.5)') u
  write(*,*) 'r'
  write(*,'(5f11.5)') r
  write(*,*) 'ax by cc ay by'
  write(*,'(5f11.5)') ax(2,2),bx(2,2),cc(2,2),ay(2,2),ay(2,2)
  y(:,:) = u(:,:) 

  if (ax(2,2) .le. ay(2,2)) then
  do j = 2, m-1
  do i = 2, n-1
    y(i,j+1) = (r(i,j) -ax(i,j)*y(i-1,j)-bx(i,j)*y(i+1,j)-cc(i,j)*y(i,j)-ay(i,j)*y(i,j-1))/by(i,j)
  end do
  end do

  do i = 2, n-1
  do j = 2, m-1
  u(i,2) = u(i,2) +(y(j,m)-u(j,m))*rinv(i-1,j-1)
  end do
  end do
  do j = 2, m-2
  do i = 2, n-1
    u(i,j+1) = (r(i,j) -ax(i,j)*u(i-1,j)-bx(i,j)*u(i+1,j)-cc(i,j)*u(i,j)-ay(i,j)*u(i,j-1))/by(i,j)
  end do
  end do

  else 

  do i = 2, n-1
  do j = 2, m-1
    y(i+1,j) = (r(i,j) -ax(i,j)*y(i-1,j)-by(i,j)*y(i,j+1)-cc(i,j)*y(i,j)-ay(i,j)*y(i,j-1))/bx(i,j)
  end do
  end do

  do j = 2, m-1
  do i = 2, n-1
  u(2,j) = u(2,j) +(y(n,i)-u(n,i))*rinv(i-1,j-1)
  end do
  end do
  do i = 2, n-2
  do j = 2, m-1
    u(i+1,j) = (r(i,j) -ax(i,j)*u(i-1,j)-by(i,j)*u(i,j+1)-cc(i,j)*u(i,j)-ay(i,j)*u(i,j-1))/bx(i,j)
  end do
  end do
end if


  write(*,*) 'after evp'
  write(*,'(5f11.4)') u

  call calc_rhs(ax,bx,cc,ay,by,u,r,y,n,m,rr)
  write(*,*) 'rr in evp  :', rr

end subroutine 
subroutine cgstep(ax,bx,cc,ay,by,u,r,n,m)

implicit none 
integer :: n,m
real*8,dimension(n,m) :: ax,bx,cc,ay,by,u,r
! local 

real*8,dimension(n,m) :: p,ap,s
real*8 :: mu,nu,rho

p(:,:) = r(:,:)
call vectornorm(r,r,mu,n,m)
call matrixmultiple(ax,bx,cc,ay,by,p,ap,n,m)
call vectornorm(p,ap,rho,n,m)
write(*,*) 'mu,rho,mu/rho',mu,rho,mu/rho
u(:,:) = u(:,:) + mu/rho*p(:,:)
end

subroutine pcg(ax,bx,cc,ay,by,u,f,n,m,tol)

implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by,f
real*8,dimension(n,m),intent(inout) :: u
real*8,intent(in) :: tol
! local 
integer:: iter
real*8,dimension(n,m) :: p,ap,s,au,r
real*8 :: mu,nu,rho,alpha

call matrixmultiple(ax,bx,cc,ay,by,u,au,n,m)
r(:,:) = f(:,:) - au(:,:)
p(:,:) = r(:,:) 
call vectornorm(r,r,mu,n,m)

iter = 0
do while ((mu > tol) .and. (iter < 80))
  call matrixmultiple(ax,bx,cc,ay,by,p,ap,n,m)
  write(*,*) ap
  call vectornorm(p,ap,rho,n,m)
  alpha = mu/rho
  write(*,*) 'mu,rho,alpha',mu,rho,alpha
  u(:,:) = u(:,:) + 1.*alpha*p(:,:)
  r(:,:) = r(:,:) - alpha*ap(:,:)
  call  vectornorm(r,r,nu,n,m)
  p(:,:) = r(:,:)+nu/mu*p(:,:)
  mu = nu
  iter = iter +1
  write(*,*) 'CG iter', iter,mu
end do 
end
subroutine diagpcg(ax,bx,cc,ay,by,u,f,n,m,tol)

implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by,f
real*8,dimension(n,m),intent(inout) :: u
real*8,intent(in) :: tol
! local 
integer:: iter
real*8,dimension(n,m) :: p,ap,s,au,r,z
real*8 :: mu,nu,rho,alpha,rr

call matrixmultiple(ax,bx,cc,ay,by,u,au,n,m)
r(:,:) = f(:,:) - au(:,:)
z(:,:) = r(:,:)/cc(:,:) 
p(:,:) = z(:,:) 
call vectornorm(r,z,mu,n,m)
iter = 0
rr = 1.0
do while ((rr > tol) .and. (iter < 80))
  call matrixmultiple(ax,bx,cc,ay,by,p,ap,n,m)
  write(*,*) ap
  call vectornorm(p,ap,rho,n,m)
  alpha = mu/rho
  write(*,*) 'mu,rho,alpha',mu,rho,alpha
  u(:,:) = u(:,:) + alpha*p(:,:)
  r(:,:) = r(:,:) - alpha*ap(:,:)
  z(:,:) = r(:,:)/cc(:,:)
  call  vectornorm(r,r,rr,n,m)
  !call lev_evp(ax,bx,cc,ay,by,rinv,z,r,n,nb,m,mb,nn)
  call  vectornorm(z,r,nu,n,m)
  call  vectornorm(r,r,rr,n,m)
  p(:,:) = z(:,:)+nu/mu*p(:,:)
  iter = iter +1
  mu = nu
  write(*,*) 'CG ',iter,mu,rr,tol
end do 
end
subroutine evppcg(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn,tol)
  implicit none 
  integer :: n,m,nb,mb,nn
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by,f
  real*8,dimension(nb,mb,nn-2,nn-2),intent(in) :: rinv
  real*8,dimension(n,m),intent(inout) :: u
  real*8,intent(in) :: tol
  ! local 
  integer:: iter
  real*8,dimension(n,m) :: p,ap,s,au,r,z
  real*8 :: mu,nu,rho,alpha,rr

  call lev_evp(ax,bx,cc,ay,by,rinv,u,f,n,nb,m,mb,nn)
  call matrixmultiple(ax,bx,cc,ay,by,u,au,n,m)
  r(:,:) = f(:,:) - au(:,:)
  z(:,:) = r(:,:) 
  call lev_evp(ax,bx,cc,ay,by,rinv,z,r,n,nb,m,mb,nn)
  p(:,:) = z(:,:) 
  call vectornorm(r,z,mu,n,m)
  iter = 0
  rr = 1.0
  do while ((rr > tol) .and. (iter < 80))
    call matrixmultiple(ax,bx,cc,ay,by,p,ap,n,m)
    write(*,*) ap
    call vectornorm(p,ap,rho,n,m)
    alpha = mu/rho
    write(*,*) 'mu,rho,alpha',mu,rho,alpha
    u(:,:) = u(:,:) + alpha*p(:,:)
    r(:,:) = r(:,:) - alpha*ap(:,:)
    z(:,:) = r(:,:)
    call  vectornorm(r,r,rr,n,m)
    call lev_evp(ax,bx,cc,ay,by,rinv,z,r,n,nb,m,mb,nn)
    call  vectornorm(z,r,nu,n,m)
    call  vectornorm(r,r,rr,n,m)
    p(:,:) = z(:,:)+nu/mu*p(:,:)
    iter = iter +1
    mu = nu
    write(*,*) 'CG ',iter,mu,rr,tol
  end do 
end

subroutine vectornorm(u,v,rr,n,m)
implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: u,v
real*8, intent(inout) :: rr

! local 
integer :: i,j
rr = 0.
do i = 2, n-1 
  do j = 2, m-1
    rr = rr + u(i,j)*v(i,j) 
   end do 
   end do 
   end subroutine 
subroutine matrixmultiple(ax,bx,cc,ay,by,u,s,n,m)
implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by,u
real*8,dimension(n,m),intent(inout) :: s

! local 
integer :: i,j
s(:,:) = 0.
do i = 2, n-1 
  do j = 2, m-1
     s(i,j) =  ax(i,j)*u(i-1,j)+by(i,j)*u(i,j+1)+cc(i,j)*u(i,j)+ay(i,j)*u(i,j-1)+bx(i,j)*u(i+1,j)
   end do 
   end do 
   end subroutine 


  subroutine inverse(a,c,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed 
  ! during the calculation
  !===========================================================
  implicit none 
  integer,intent(in) ::  n
  real*8,intent(inout) :: a(n,n)
  real*8,intent(inout) :: c(n,n)
  real*8 :: L(n,n), U(n,n), b(n), d(n), x(n)
  real*8 :: coeff
  integer::  i, j, k

  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 aloows such operations on matrices
  L=0.0
  U=0.0
  b=0.0

  ! step 1: forward elimination
  do k=1, n-1
    do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
        a(i,j) = a(i,j)-coeff*a(k,j)
      end do
    end do
  end do

  ! Step 2: prepare L and U matrices 
  ! L matrix is a matrix of the
  ! elimination coefficient
  ! + the diagonal elements are 1.0
  do i=1,n
    L(i,i) = 1.0
  end do
  ! U matrix is the upper triangular
  ! part of A
  do j=1,n
    do i=1,j
      U(i,j) = a(i,j)
    end do
  end do

  ! Step 3: compute
  ! columns of the inverse
  ! matrix C
  do k=1,n
    b(k)=1.0
    d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
      d(i)=b(i)
      do j=1,i-1
        d(i)=d(i)-L(i,j)*d(j)
      end do
    end do
    ! Step 3b:Solve Ux=d using the back substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
    end do
    ! Step  3c: fill the solutions  x(n) into column  k  of  C
    do  i=1,n
      c(i,k) = x(i)
    end do

    b(k)=0.0
  end do 
  end subroutine





