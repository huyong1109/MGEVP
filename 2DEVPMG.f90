! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*8,parameter :: rtol = 1.0e-8
  integer*8, parameter :: nlev = 2
  integer*8, parameter :: mlev = 2
  integer*8, parameter :: choice = 1
  integer*8, parameter :: ngrid = 4**nlev+1  ! grid num on finest level 
  integer*8, parameter :: mgrid = 4**mlev+1  ! grid num on finest level 
  integer*8, parameter :: tngrid = 4*(4**nlev -1)/3 +nlev ! total grids all level 
  integer*8, parameter :: tmgrid = 4*(4**mlev -1)/3 +mlev ! total grids all level 
  !integer, parameter :: tingrid = (8*(4**(nlev-1) -1) &     ! total inner grid 
  !    +21*(nlev-1))/9  +1-2*(nlev- 1)
  real*8  :: tol 
  real*8,dimension(tngrid,tmgrid) :: r, f,u,ut
  real*8,dimension(tngrid,tmgrid) :: ax,bx,cc,ay,by
  real*8,dimension(tngrid -ngrid +2,tmgrid-mgrid+2,3,3) :: rinv

  ! function 
  !integer :: grid_num

  ! local variable
  real*8 :: s,ntmp,mtmp
  integer*8 :: i,j,ii,jj,kn,km,kn1,km1,hn,hm,n,tn,ln,m,tm,lm
  integer*8 :: ngrid2,mgrid2

  tol = rtol/dble(ngrid)**2
  print *, '2D Multigrid EVP solver'

  ! set problem 
  kn = 1
  hn =1

  do i = nlev,1,-1

    call grid_num(i,n,tn,ln)
    ngrid2 = 2*(nlev-i)
    ngrid2 = 4**ngrid2
    mtmp = 1.0/dble(ngrid2)
    write(*,'(A,I4.1,a,i11.1,a,i12.1,a,i10.1,a,i10.1)') 'nlev',i, '  n',n,'  tn',tn,'  axs', kn, '  rinvs', hn

    km = 1
    hm =1
    do j = mlev,1,-1
      call grid_num(j,m,tm,lm)
      mgrid2 = 2*(mlev-j)
      mgrid2 = 4**mgrid2
      ntmp = 1.0/dble(mgrid2)
      !write(*,*) 'tmp', tmp
      kn1 = kn+n-1
      km1 = km+m-1
      do ii = kn, kn1
        do jj = km, km1
          ax(ii,jj) = mtmp
          bx(ii,jj) = mtmp
          ay(ii,jj) = ntmp
          by(ii,jj) = ntmp
          cc(ii,jj) = -ax(ii,jj)-bx(ii,jj)-ay(ii,jj)-by(ii,jj)
        end do 
      end do 

      write(*,'(A,I4.1,a,i11.1,a,i12.1,a,i10.1,a,i10.1)') 'mlev',j, '  m',m,'  tm',tm,'  axs', km, '  rinvs', hm
      write(*,*) kn, kn1, km,km1,hn,ln,hm,lm
      call lev_pre(ax(kn:kn1,km:km1),bx(kn:kn1,km:km1),cc(kn:kn1,km:km1),ay(kn:kn1,km:km1),by(kn:kn1,km:km1),rinv(hn:hn+ln-2,hm:hm+lm-2,:,:),n,ln,m,lm)

      !write(*,*) ax(j:j+n-1)
      !write(*,*) cx(j:j+n-1)
      !write(*,*) bx(j:j+n-1)
      !write(*,*) rinv(k:k+ln-1)
      km = km+m 
      hm = hm +lm
    end do 
    kn = kn+n 
    hn = hn +ln
  end do 

  call analytic_init(f,ngrid,mgrid,choice)


  ! initalization 
  do i = 1, tmgrid
    do j = 1, tngrid
      r(i,j) = 0.0
      u(i,j) = 0.0
      ut(i,j) = 0.0
    end do 
  end do 

  ! traditional iteration 
  call grid_num(nlev,n,tn,ln)
  call grid_num(mlev,m,tm,lm)
  call MG_x(ax,bx,cc,ay,by,rinv,u,f,nlev,n,tn,ln,mlev,m,tm,tol)
  call analytic_check(u(1:n,1:m),ngrid,mgrid,choice)



  end program

  subroutine analytic_init(f,n,m,choice)
  implicit none
  integer*8,intent(in):: n,m
  integer,intent(in):: choice
  real*8, dimension(n,m), intent(inout) :: f

  ! local 
  integer*8 :: i,j
  real*8 :: x,y,hm,hn,hm2,hn2
  real*8,parameter :: pi = 4.*atan(1.) 

  ! input right side
  f(:,:) = 0.0
  hm = 1.0/dble(n-1)
  hn = 1.0/dble(m-1)
  hm2= hm**2
  hn2= hn**2

  select case(choice)
  case (1)
    do i = 2, n-1
      do j = 2, m-1
        x = (i-1)*hm
        y = (j-1)*hn
        f(i,j) = -(x*(1.0-x)+y*(1.0-y))
      end do 
    end do 
  case (2)

    do i = 2, n-1
      do j = 2, m-1
        x = (i-1)*hm
        y = (j-1)*hn
        f(i,j) = -8.0*pi**2*sin(2*pi*x)*sin(2*pi*y)
      end do 
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select

  !write(*,*) 'init f(:)'
  !write(*,*) f(:)

  end subroutine
  subroutine analytic_check(u,n,m,choice)
  implicit none 
  integer*8,intent(in):: n,m
  integer,intent(in):: choice
  real*8, dimension(n,m), intent(in) :: u

  ! local 
  integer*8 :: i,j
  real*8 :: hm,hn,x,y,e,maxe
  real*8,parameter :: pi = 4.*atan(1.) 


  maxe = 0.0
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)
  select case(choice)
  case (1)
    !=========== for u'' = -1=========
    do i = 1, n
      do j = 1, m
        x = (i-1)*hn
        y = (j-1)*hm
        e = u(i,j) -0.5*x*(1.0-x)*y*(1.0-y)
        maxe= max(abs(e),maxe) 
        !write(*,*) i, s, x(i), e
      end do 
    end do 

  case (2)
    !=========== for u'' = -1=========
    do i = 1, n
      do j = 1, m
        x = (i-1)*hn
        y = (j-1)*hm
        e = u(i,j) -sin(2.0*pi*x)*sin(2.0*pi*y)
        maxe= max(abs(e),maxe) 
        !write(*,*) i, s, x(i), e
      end do 
    end do 
  case default
    write(*,*) 'Unsupported initial case!'
  end select
  write(*,*) 'MAX ERROR(relative to analytic): ',maxe

  end subroutine


  recursive subroutine MG_x(ax,bx,cc,ay,by,rinv,u,r,nlev,n,tn,ln,mlev,m,tm,lm,tol)
  implicit none
  ! traditional iteration 
  integer*8, intent(in) :: nlev,mlev
  integer*8, intent(in) :: n,tn,ln,m,tm,lm
  !real*8, intent(inout) :: rr
  integer*8, parameter :: nn = 3
  real*8,dimension(tn,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(tn-n+2,tm-m+2,nn,nn),intent(in) :: rinv
  real*8,dimension(tn,tm),intent(inout) :: r,u
  real*8,intent(in) :: tol

  integer*8 :: i,j,k,j1,k1,j1n,k1n,n1,tn1,ln1
  integer :: iter
  real*8 :: rr
  real*8,dimension(tn,tm) :: f,ut

  j=1
  k =1
  !f(j:j+n-1) = r(j:j+n-1)
  !call calc_rhs(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),x(j:j+n-1),f(j:j+n-1),r(j:j+n-1),n,rr)

  ut(:,:) = 0.
  call lev_evp_x(ax(j:j+n-1,:),bx(j:j+n-1,:),cc(j:j+n-1,:),ay(j:j+n-1,:),by(j:j+n-1,:),rinv(k:k+ln-2,:,:,:),ut(j:j+n-1,:),r(j:j+n-1,:),n,ln,mlev,m,tm,lm,tol)
  f(j:j+n-1,:) = r(j:j+n-1,:)
  u(j:j+n-1,:) = u(j:j+n-1,:)+ut(j:j+n-1,:)
  call calc_rhs(ax(j:j+n-1,:),bx(j:j+n-1,:),cc(j:j+n-1,:),ay(j:j+n-1,:),by(j:j+n-1,:),ut(j:j+n-1,:),f(j:j+n-1,:),r(j:j+n-1,:),n,m,rr)

  write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_x',nlev, 'rr= ', rr, 'tol= ',tol
  !write(*,*) 'r'
  !write(*,'(17f7.4)') r(j:j+n-1) 
  !write(*,*) 'x'
  !write(*,'(17f7.4)') x(j:j+n-1) 

  if (nlev .ne. 1) then
    call grid_num(nlev-1,n1,tn1,ln1)
    j1 = j +n
    j1n = j1 +tn1
    k1 = k +ln
    k1n = k1 +tn1-n1+2
    iter  = 0
    do while ((rr .gt. tol) .and. (iter < 2) )
      !do while (rr .gt. tol)
      iter = iter +1
      ! coarsing 
      call coarse_rhs_x(r(j:j+n+ln-1,:), n,ln,m)
      ut(:,:) = 0.
      write(*,*) 'MGX :',nlev
      call MG_x(ax(j1:j1n,:),bx(j1:j1n,:),cc(j1:j1n,:),ay(j1:j1n,:),by(j1:j1n,:),rinv(k1:k1n,:,:,:),ut(j1:j1n,:),r(j1:j1n,:),nlev-1,n1,tn1,ln1,mlev,m,tm,lm,tol)
      ! relaxing
      call finer_x(ut(1:n+ln,:),n,ln,m) 
      call lev_evp_x(ax(j:j+n-1,:),bx(j:j+n-1,:),cc(j:j+n-1,:),ay(j:j+n-1,:),by(j:j+n-1,:),rinv(k:k+ln-2,:,:,:),ut(j:j+n-1,:),r(j:j+n-1,:),n,ln,mlev,m,tm,lm,tol)
      f(j:j+n-1,:) = r(j:j+n-1,:)
      u(j:j+n-1,:) = u(j:j+n-1,:)+ut(j:j+n-1,:)
      call calc_rhs(ax(j:j+n-1,:),bx(j:j+n-1,:),cc(j:j+n-1,:),ay(j:j+n-1,:),by(j:j+n-1,:),ut(j:j+n-1,:),f(j:j+n-1,:),r(j:j+n-1,:),n,m,rr)
      write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_x',nlev, 'rr= ', rr, 'tol= ',tol
      !write(*,*) 'r'
      !write(*,'(17f7.4)') r(j:j+n-1) 
      !write(*,*) 'x'
      !write(*,'(17f7.4)') x(j:j+n-1) 

    end do 
  end if 



  end subroutine


  recursive subroutine MG_y(ax,bx,cc,ay,by,rinv,u,r,mlev,m,tm,lm,tol)
  implicit none
  ! traditional iteration 
  integer*8, intent(in) :: mlev
  integer*8, intent(in) :: m,tm,lm
  !real*8, intent(inout) :: rr
  integer*8,parameter :: n = 5
  real*8,dimension(n,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(tm-m+2,n-2,n-2),intent(in) :: rinv
  real*8,dimension(n,tm),intent(inout) :: r,u
  real*8,intent(in) :: tol

  integer*8 :: i,j,k,j1,k1,j1n,k1n,m1,tm1,lm1
  integer :: iter
  real*8 :: rr
  real*8,dimension(n,tm) :: f,ut

  j=1
  k =1
  !f(j:j+n-1) = r(j:j+n-1)
  !call calc_rhs(ax(j:j+n-1),cx(j:j+n-1),bx(j:j+n-1),x(j:j+n-1),f(j:j+n-1),r(j:j+n-1),n,rr)

  ut(:,:) = 0.
  call lev_evp_y(ax(:,j:j+m-1),bx(:,j:j+m-1),cc(:,j:j+m-1),ay(:,j:j+m-1),by(:,j:j+m-1),rinv(k:k+lm-2,:,:),ut(:,j:j+m-1),r(:,j:j+m-1),m,lm)
  f(:,j:j+m-1) = r(:,j:j+m-1)
  u(:,j:j+m-1) = u(:,j:j+m-1)+ut(:,j:j+m-1)
  call calc_rhs(ax(:,j:j+m-1),bx(:,j:j+m-1),cc(:,j:j+m-1),ay(:,j:j+m-1),by(:,j:j+m-1),rinv(k:k+lm-2,:,:),ut(:,j:j+m-1),f(:,j:j+m-1),r(:,j:j+m-1),m,rr)

  write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_y',mlev, 'rr= ', rr, 'tol= ',tol
  !write(*,*) 'r'
  !write(*,'(17f7.4)') r(j:j+n-1) 
  !write(*,*) 'x'
  !write(*,'(17f7.4)') x(j:j+n-1) 

  if (mlev .ne. 1) then
    call grid_num(mlev-1,m1,tm1,lm1)
    j1 = j +m
    j1n = j1 +tm1
    k1 = k +lm
    k1n = k1 +tm1-m1+2
    iter  = 0
    do while ((rr .gt. tol) .and. (iter < 2) )
      !do while (rr .gt. tol)
      iter = iter +1
      ! coarsing 
      call coarse_rhs_y(r(:,j:j+m+lm-1), m,lm)
      ut(:,:) = 0.
      write(*,*) 'MGY ',mlev
      call MG_y(ax(:,j1:j1n),bx(:,j1:j1n),cc(:,j1:j1n),ay(:,j1:j1n),by(:,j1:j1n),rinv(k1:k1n,:,:),ut(:,j1:j1n),r(:,j1:j1n),mlev-1,m1,tm1,lm1,tol)
      ! relaxing
      call finer_y(ut(:,1:m+lm),m,lm) 
      call lev_evp_y(ax(:,j:j+m-1),bx(:,j:j+m-1),cc(:,j:j+m-1),ay(:,j:j+m-1),by(:,j:j+m-1),rinv(k:k+lm-2,:,:),ut(:,j:j+m-1),r(:,j:j+m-1),m,lm)
      f(:,j:j+m-1) = r(:,j:j+m-1)
      u(:,j:j+m-1) = u(:,j:j+m-1)+ut(:,j:j+m-1)
      call calc_rhs(ax(:,j:j+m-1),bx(:,j:j+m-1),cc(:,j:j+m-1),ay(:,j:j+m-1),by(:,j:j+m-1),rinv(k:k+lm-2,:,:),ut(:,j:j+m-1),f(:,j:j+m-1),r(:,j:j+m-1),m,rr)
      write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_y',mlev, 'rr= ', rr, 'tol= ',tol
      !write(*,*) 'r'
      !write(*,'(17f7.4)') r(j:j+n-1) 
      !write(*,*) 'x'
      !write(*,'(17f7.4)') x(j:j+n-1) 

    end do 
  end if 



  end subroutine


  subroutine coarse_rhs_x(r,n,ln,m)
  implicit none
  integer*8,intent(in) :: n,ln,m
  real*8,dimension(n+ln,m),intent(inout) :: r

  !local 
  integer*8 :: i,j
  j = 2
  r(n+1:n+ln,:) = 0.0
  do i = 5, n-1, 4
    r(n+j,:) = 0.25*r(i,:)
    j = j+1
  end do
  end subroutine 

  subroutine coarse_rhs_y(r,n,m,lm)
  implicit none
  integer*8,intent(in) :: n,lm,m
  real*8,dimension(n,m+lm),intent(inout) :: r

  !local 
  integer*8 :: i,j
  j = 2
  r(:,m+1:m+lm) = 0.0
  do i = 5, n-1, 4
    r(:,m+j) = 0.25*r(:,i)
    j = j+1
  end do
  end subroutine 

  subroutine finer_x(u,n,ln,m)
  implicit none
  integer*8,intent(in) :: n,ln,m
  real*8,dimension(n+ln,m),intent(inout) :: u

  !local 
  integer*8 :: i,j
  j = 2
  do i = 5, n-1, 4
    u(i,:) = u(i,:) +u(n+j,:)
    j = j+1
  end do
  end subroutine 
  subroutine finer_y(u,n,m,lm)
  implicit none
  integer*8,intent(in) :: n,lm,m
  real*8,dimension(n,m+lm),intent(inout) :: u

  !local 
  integer*8 :: i,j
  j = 2
  do i = 5, m-1, 4
    u(:,i) = u(:,i) +u(:,m+j)
    j = j+1
  end do
  end subroutine 

  subroutine calc_rhs(ax,bx,cc,ay,by,u,f,r,n,m,rr)
  implicit none
  integer*8,intent(in) :: n,m
  real*8, intent(inout) :: rr
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: u,f
  real*8,dimension(n,m),intent(inout) :: r

  !local 
  integer*8 :: i,j
  rr = 0.0
  do i = 2, n-1
    do j = 2, m-1
      r(i,j) = f(i,j) -ax(i,j)*u(i-1,j)-bx(i+1,j)*u(i+1,j)-cc(i,j)*u(i,j)-ay(i,j)*u(i,j-1)-by(i,j)*u(i,j+1)
      rr = rr + abs(r(i,j)**2)
    end do
  end do
  rr = sqrt(rr/dble(n*m))
  end subroutine 

! 5 points PRE
subroutine pre(ax,bx,cc,ay,by,rinv)
  implicit none
  integer,parameter :: nn = 5
  integer,parameter :: mm = 5
  real*8,dimension(nn,mm),intent(in) :: ax,bx,cc,ay,by

  real*8,dimension(nn-2,nn-2),intent(inout) :: rinv

  !local 
  integer :: i,j,ii,info
  real*8,dimension(nn,mm) :: y
  real*8,dimension(nn-2) :: work,ipiv
  real*8,dimension(nn-2,nn-2) :: rin

  y(:,:) = 0.0
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

  write(*,*) 'rinv',rinv(:,:)
  rin(:,:) = rinv(:,:)
  call inverse(rin,rinv,nn-2)

  !call ZGETRF(nn-2,nn-2,rinv,nn-2,IPIV,info)
  !if(info .eq. 0) then
  !  write(*,*)"succeded"
  !else
  !  write(*,*)"failed"
  !end if

  !call ZGETRI(nn-2,rinv,nn-2,IPIV,WORK,nn-2,info)
  !if(info .eq. 0) then
  !  write(*,*)"succeded"
  !else
  !  write(*,*)"failed"
  !end if
  
  !! check pre rinv 
  write(*,*) 'rinv',rinv(:,:)
  work(:) = 0.0
  do j = 1,nn-2
    do i = 1,nn-2
      work(j) = work(j) + rinv(i,j)*y(i+1,mm)
    end do 
  end do 
  write(*,*) 'work',work(:)

end subroutine 
subroutine lev_pre(ax,bx,cc,ay,by,rinv,n,ln,m,lm)
  implicit none
  integer*8,intent(in) :: n,ln,m,lm! grid number of current and lower (coarser) levels
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(ln-1,lm-1,3,3),intent(inout) :: rinv

  !local 
  integer*8 :: i,j,i4,j4,k,l
  k = 0
  do i = 1, n-1, 4
    k = k+1
    l = 0
    do j = 1, m-1, 4
      l = l+1
      write(*,*) i,j,k,l
      i4 =i+4 
      j4 =j+4
      call pre(ax(i:i4,j:j4), bx(i:i4,j:j4),cc(i:i4,j:j4),ay(i:i4,j:j4),by(i:i4,j:j4),rinv(k,l,:,:))
  end do
  end do

end subroutine 

! solve each block in one level
subroutine lev_evp_x(ax,bx,cc,ay,by,rinv,u,r,n,ln,mlev,m,tm,lm,tol)
  implicit none
  integer*8,intent(in) :: n,ln,mlev,m,tm,lm
  real*8,dimension(n,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,tm),intent(in) :: r
  integer*8, parameter :: nn = 3
  real*8,dimension(n,tm-m+2,nn,nn),intent(in) :: rinv
  real*8,dimension(n,tm),intent(inout) :: u
  real*8,intent(in) :: tol

  !local 
  integer*8 :: i,j
  j = 0
  do i = 1, n-1, 4
    j = j+1
    !write(*,*) 'j=',j,'i=',i,'i+4',i+4
   call MG_y(ax(i:i+4,:), bx(i:i+4,:),cc(i:i+4,:),ay(i:i+4,:),by(i:i+4,:),rinv(j,:,:,:),u(i:i+4,:),r(i:i+4,:),mlev,m,tm,lm,tol)

  end do
  !write(*,'(17f7.4)') x(1:n-1) 
  !do i = 5, n-1, 4
  ! !x(i) = 0.25*(2.0*x(i) + x(i-1)+x(i+1))
  ! x(i) = 0.5*( x(i-1)+x(i+1))
  !end do
  !write(*,*) 'average'
  !write(*,'(17f7.4)') x(1:n-1) 
end subroutine 

! solve each block in one level
subroutine lev_evp_y(ax,bx,cc,ay,by,rinv,x,r,m,lm)
  implicit none
  integer*8,intent(in) :: m,lm
  integer, parameter :: n =5
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(lm-1,n-2,n-2),intent(in) :: rinv
  real*8,dimension(n,m),intent(inout) :: x

  !local 
  integer*8 :: i,j
  j = 0
  do i = 1, n-1, 4
    j = j+1
    !write(*,*) 'j=',j,'i=',i,'i+4',i+4
   call evp(ax(:,i:i+4), bx(:,i:i+4),cc(:,i:i+4),ay(:,i:i+4),by(:,i:i+4),rinv(j,:,:),x(:,i:i+4),r(:,i:i+4))

  end do
  !write(*,'(17f7.4)') x(1:n-1) 
  !do i = 5, n-1, 4
  ! !x(i) = 0.25*(2.0*x(i) + x(i-1)+x(i+1))
  ! x(i) = 0.5*( x(i-1)+x(i+1))
  !end do
  !write(*,*) 'average'
  !write(*,'(17f7.4)') x(1:n-1) 
end subroutine 

! 5 points EVP 
subroutine evp(ax,bx,cc,ay,by,rinv,x,r)
  implicit none
  integer,parameter :: n = 5
  integer,parameter :: m = 5
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(n,m),intent(inout) :: x

  real*8,dimension(n-2,n-2),intent(in) :: rinv

  !local 
  integer :: i,j
  real*8,dimension(n,m) :: y

  y(:,:) = x(:,:) 
  do j = 2, m-1
  do i = 2, n-1
    y(i,j+1) = (r(i,j) -ax(i,j)*y(i-1,j)-bx(i,j)*y(i+1,j)-cc(i,j)*y(i,j)-ay(i,j)*y(i,j-1))/by(i,j)
  end do
  end do

  do i = 1, n-2
  do j = 1, n-2
  x(i,2) = x(i,2) + (y(i,m)-x(i,m))*rinv(j,i)
  end do
  end do
  do j = 2, m-1
  do i = 2, n-1
    x(i,j+1) = (r(i,j) -ax(i,j)*x(i-1,j)-bx(i,j)*x(i+1,j)-cc(i,j)*x(i,j)-ay(i,j)*x(i,j-1))/by(i,j)
  end do
  end do

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
  integer i, j, k

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
    ! Step 3a: Solve
    ! Ld=b using the
    ! forward
    ! substitution
    do i=2,n
      d(i)=b(i)
      do j=1,i-1
        d(i)=d(i)-L(i,j)*d(j)
      end do
    end do
    ! Step
    ! 3b:
    ! Solve
    ! Ux=d
    ! using
    ! the
    ! back
    ! substitution
    x(n)=d(n)/U(n,n)
    do i = n-1,1,-1
      x(i) = d(i)
      do j=n,i+1,-1
        x(i)=x(i)-U(i,j)*x(j)
      end do
      x(i) = x(i)/u(i,i)
    end do
    ! Step
    ! 3c:
    ! fill
    ! the
    ! solutions
    ! x(n)
    ! into
    ! column
    ! k
    ! of
    ! C
    do  i=1,n
      c(i,k) = x(i)
    end do

    b(k)=0.0
  end do 
  end subroutine





