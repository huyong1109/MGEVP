! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*8,parameter :: rtol = 1.0e-8
  integer, parameter :: nlev = 2
  integer, parameter :: mlev = 2
  integer, parameter :: choice = 1
  integer, parameter :: ngrid = 4**nlev+1  ! grid num on finest level 
  integer, parameter :: mgrid = 4**mlev+1  ! grid num on finest level 
  integer, parameter :: tngrid = 4*(4**nlev -1)/3 +nlev ! total grids all level 
  integer, parameter :: tmgrid = 4*(4**mlev -1)/3 +mlev ! total grids all level 
  !integer, parameter :: tingrid = (8*(4**(nlev-1) -1) &     ! total inner grid 
  !    +21*(nlev-1))/9  +1-2*(nlev- 1)
  real*8  :: tol 
  real*8,dimension(ngrid,mgrid) :: r,f,u,ut
  real*8,dimension(tngrid,tmgrid) :: ax,bx,cc,ay,by
  real*8,dimension(tngrid -ngrid +2,tmgrid-mgrid+2,3,3) :: rinv

  ! function 
  !integer :: grid_num

  ! local variable
  real*8 :: s,ntmp,mtmp
  integer :: i,j,ii,jj,kn,km,kn1,km1,hn,hm,n,tn,ln,m,tm,lm
  integer :: ngrid2,mgrid2

  tol = rtol/dble(ngrid)**2
  print *, '2D Multigrid EVP solver'

  ! set problem 
  kn = 1
  hn =1

  do i = nlev,1,-1

    call grid_num(i,n,tn,ln)
    ngrid2 = 2*i
    ngrid2 = 4**ngrid2
    ntmp = dble(ngrid2)
    write(*,'(A,I4.1,a,i11.1,a,i12.1,a,i10.1,a,i10.1)') 'nlev',i, '  n',n,'  tn',tn,'  axs', kn, '  rinvs', hn

    km = 1
    hm =1
    do j = mlev,1,-1
      call grid_num(j,m,tm,lm)
      mgrid2 = 2*j
      mgrid2 = 4**mgrid2
      mtmp = dble(mgrid2)
     ! write(*,*) 'tmp', ntmp,mtmp
     ! write(*,*) 'grid', ngrid2, mgrid2
      kn1 = kn+n-1
      km1 = km+m-1
      do ii = kn, kn1
        do jj = km, km1
          ax(ii,jj) = ntmp
          bx(ii,jj) = ntmp
          ay(ii,jj) = mtmp
          by(ii,jj) = mtmp
          cc(ii,jj) = -ax(ii,jj)-bx(ii,jj)-ay(ii,jj)-by(ii,jj)
        end do 
      end do 

      write(*,'(A,I4.1,a,i11.1,a,i12.1,a,i10.1,a,i10.1)') 'mlev',j, '  m',m,'  tm',tm,'  axs', km, '  rinvs', hm
      !write(*,*) kn, kn1, km,km1,hn,ln,hm,lm
      call lev_pre(ax(kn:kn1,km:km1),bx(kn:kn1,km:km1),cc(kn:kn1,km:km1),ay(kn:kn1,km:km1),by(kn:kn1,km:km1),rinv(hn:hn+ln-2,hm:hm+lm-2,:,:),n,ln,m,lm)
     ! write(*,*) 'ax,ay'
     ! write(*,*) ax
     ! write(*,*) bx
     ! write(*,*) cc
     ! write(*,*) ay
     ! write(*,*) by
     ! write(*,*) rinv
      km = km+m 
      hm = hm +lm
    end do 
    kn = kn+n 
    hn = hn +ln
  end do 

  call analytic_init(f,ngrid,mgrid,choice)


  ! initalization 
  do i = 1, ngrid
    do j = 1, mgrid
      r(i,j) = 0.0
      u(i,j) = 0.0
      ut(i,j) = 0.0
    end do 
  end do 

  ! traditional iteration 
  call grid_num(nlev,n,tn,ln)
  call grid_num(mlev,m,tm,lm)
  !write(*,*) 'u'
  !write(*,*) u
  write(*,*) 'f'
  write(*,'(5f7.3)') f(:,:)
  !write(*,*) f(:,:)
  call MG_x(ax,bx,cc,ay,by,rinv,u,f,nlev,n,tn,ln,mlev,m,tm,lm,tol)
  call analytic_check(u(1:n,1:m),ngrid,mgrid,choice)



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
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)

  select case(choice)
  case (1)
    do i = 2, n-1
        x = dble(i-1)*hn
      do j = 2, m-1
        y = dble(j-1)*hm
        f(i,j) = -(x*(1.0-x)+y*(1.0-y))
      end do 
    end do 
  case (2)

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
  hn = 1.0/dble(n-1)
  hm = 1.0/dble(m-1)
  select case(choice)
  case (1)
    !=========== for u'' = -1=========
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

  case (2)
    !=========== for u'' = -1=========
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


  recursive subroutine MG_x(ax,bx,cc,ay,by,rinv,u,r,nlev,n,tn,ln,mlev,m,tm,lm,tol)
  implicit none
  integer, intent(in) :: nlev,mlev
  integer, intent(in) :: n,tn,ln,m,tm,lm
  integer, parameter :: nn = 3
  real*8,dimension(tn,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(tn-n+2,tm-m+2,nn,nn),intent(in) :: rinv
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(n,m),intent(inout) :: u
  real*8,intent(in) :: tol
  
  !local 
  integer :: i,j,k,j1,jn1,k1,j1n,k1n,n1,tn1,ln1
  integer :: iter
  real*8 :: rr
  real*8,dimension(n,m) :: f,ut,rt
  real*8,dimension(ln,m) :: lr,lut

  j=1
  k =1

  ut(:,:) = u(:,:)
  jn1=j+n-1

  call lev_evp_x(ax(j:jn1,:),bx(j:jn1,:),cc(j:jn1,:),ay(j:jn1,:),by(j:jn1,:),rinv(k:k+ln-1,:,:,:),ut,r,n,ln,mlev,m,tm,lm,tol)

  f(:,:) = r(:,:)
  u(:,:) = ut(:,:)
  
  call calc_rhs(ax(j:jn1,1:m),bx(j:jn1,1:m),cc(j:jn1,1:m),ay(j:jn1,1:m),by(j:jn1,1:m),ut,f,rt,n,m,rr)

  write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_x_<',nlev, 'rr= ', rr, 'tol= ',tol

  !write(*,*) 'rt'
  !write(*,'(17f7.3)') rt

  if (nlev .ne. 1) then
    call grid_num(nlev-1,n1,tn1,ln1)
    j1 = j +n
    j1n = j1 +tn1-1
    k1 = k +ln
    k1n = k1 +tn1-n1+1
    iter  = 0
    !do while ((rr .gt. tol) .and. (iter < 2) )
      !do while (rr .gt. tol)
      iter = iter +1
      ! coarsing 
      call coarse_rhs_x(f,lr,n,ln,m)
     ! write(*,*) 'lr'
     ! write(*,'(5f8.3)') lr
      lut(:,:) = 0.
     write(*,*) 'MGX :',nlev
      call MG_x(ax(j1:j1n,:),bx(j1:j1n,:),cc(j1:j1n,:),ay(j1:j1n,:),by(j1:j1n,:),rinv(k1:k1n,:,:,:),lut,lr,nlev-1,n1,tn1,ln1,mlev,m,tm,lm,tol)
      ! relaxing
      call finer_x(ut,lut,n,ln,m) 
      call lev_evp_x(ax(j:jn1,:),bx(j:jn1,:),cc(j:jn1,:),ay(j:jn1,:),by(j:jn1,:),rinv(k:k+ln-1,:,:,:),ut,f,n,ln,mlev,m,tm,lm,tol)
      u(:,:) = u(:,:)+ut(:,:)
      f(:,:) = rt(:,:)
      call calc_rhs(ax(j:jn1,1:m),bx(j:jn1,1:m),cc(j:jn1,1:m),ay(j:jn1,1:m),by(j:jn1,1:m),ut,f,rt,n,m,rr)
      write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_x_>',nlev, 'rr= ', rr, 'tol= ',tol
      !write(*,*) 'r'
      !write(*,'(17f7.4)') r(j:j+n-1) 
      !write(*,*) 'x'
      !write(*,'(17f7.4)') x(j:j+n-1) 

    !end do 
  end if 

  end subroutine


  recursive subroutine MG_y(ax,bx,cc,ay,by,rinv,u,r,mlev,m,tm,lm,tol)
  implicit none
  ! traditional iteration 
  integer, intent(in) :: mlev
  integer, intent(in) :: m,tm,lm
  !real*8, intent(inout) :: rr
  integer,parameter :: nn = 5
  real*8,dimension(nn,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(tm-m+2,nn-2,nn-2),intent(in) :: rinv
  real*8,dimension(nn,m),intent(in) :: r
  real*8,dimension(nn,m),intent(inout) :: u
  real*8,intent(in) :: tol

  integer :: i,j,k,j1,k1,jm1,j1n,k1n,m1,tm1,lm1
  integer :: iter
  real*8 :: rr
  real*8,dimension(nn,m) :: f,ut,rt
  !real*8,dimension(n,lm) :: ry
  real*8,dimension(nn,lm) :: lr,lut

  !write(*,*) 'in MG_y'
  j=1
  k =1

  ut(:,:) = u(:,:)
  jm1 = j+m-1

  write(*,*) 'MGy f1 before lev_evp'
  write(*,'(5f11.5)') r

  call lev_evp_y(ax(:,j:jm1),bx(:,j:jm1),cc(:,j:jm1),ay(:,j:jm1),by(:,j:jm1),rinv(k:k+lm-1,:,:),ut,r,m,lm)
  u(:,:) = ut(:,:)
  f(:,:) = r(:,:)
  call calc_rhs(ax(:,j:jm1),bx(:,j:jm1),cc(:,j:jm1),ay(:,j:jm1),by(:,j:jm1),ut,f,rt,nn,m,rr)
  

  !write(*,*) 'MGy f1 after lev_evp'
  !write(*,'(5f11.5)') rt
  !write(*,*) 'MGy ut1'
  !write(*,'(5f11.5)') ut
  write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_y_<',mlev, 'rr= ', rr, 'tol= ',tol

  if (mlev .ne. 1) then

    call grid_num(mlev-1,m1,tm1,lm1)
    j1 = j +m
    j1n = j1 +tm1-1
    k1 = k +lm
    k1n = k1 +tm1-m1+1
    iter  = 0
    !do while ((rr .gt. tol) .and. (iter < 2) )
      !do while (rr .gt. tol)
      iter = iter +1
      ! coarsing 
      call coarse_rhs_y(f,lr,nn,m,lm)
      lut(:,:) = 0.
      !write(*,*) 'MGY ',mlev
      !write(*,*) 'lr'
      !write(*,'(5f11.4)') lr

      call MG_y(ax(:,j1:j1n),bx(:,j1:j1n),cc(:,j1:j1n),ay(:,j1:j1n),by(:,j1:j1n),rinv(k1:k1n,:,:),lut,lr,mlev-1,m1,tm1,lm1,tol)
      ! relaxing
      write(*,*) 'lut'
      write(*,'(5f7.4)') lut
      call finer_y(ut,lut,nn,m,lm) 
      write(*,*) 'finer ut'
      write(*,'(5f7.4)') ut

      call lev_evp_y(ax(:,j:jm1),bx(:,j:jm1),cc(:,j:jm1),ay(:,j:jm1),by(:,j:jm1),rinv(k:k+lm-1,:,:),ut,f,m,lm)

      u(:,:) = u(:,:)+ut(:,:)
      write(*,*) 'u2'
      write(*,'(5f7.4)') u
      f(:,:) = rt(:,:)
      call calc_rhs(ax(:,j:jm1),bx(:,j:jm1),cc(:,j:jm1),ay(:,j:jm1),by(:,j:jm1),ut,f,rt,nn,m,rr)

      write(*,'(A,I4.1,2X,A,e12.5,2X,A,e12.5)')  'LEV_y_>',mlev, 'rr= ', rr, 'tol= ',tol
      !write(*,*) 'r'
      !write(*,'(17f7.4)') r(j:j+n-1) 
      !write(*,*) 'x'
      !write(*,'(17f7.4)') x(j:j+n-1) 

   ! end do 
  end if 



  end subroutine


  subroutine coarse_rhs_x(r,lr,n,ln,m)
  implicit none
  integer,intent(in) :: n,ln,m
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(ln,m),intent(inout) :: lr

  !local 
  integer :: i,j
  j = 2
  lr(:,:) = 0.0
  do i = 5, n-1, 4
    lr(j,:) =r(i,:)
    j = j+1
  end do
  end subroutine 

  subroutine coarse_rhs_y(r,lr,n,m,lm)
  implicit none
  integer,intent(in) :: n,lm,m
  real*8,dimension(n,m),intent(in) :: r
  real*8,dimension(n,lm),intent(inout) :: lr

  !local 
  integer :: i,j
  j = 2
  lr(:,:) = 0.0
  do i = 5, m-1, 4
    lr(:,j) = r(:,i) !-ry(:,i) +0.25*(ry(:,i-1)+ry(:,i-1)+ry(:,i-2)+ry(:,i-3))
    j = j+1
  end do
  end subroutine 

  subroutine finer_x(u,lu,n,ln,m)
  implicit none
  integer,intent(in) :: n,ln,m
  real*8,dimension(ln,m),intent(in) :: lu
  real*8,dimension(n,m),intent(inout) :: u

  !local 
  integer :: i,j
  j = 2
  do i = 5, n-1, 4
    u(i,:) = u(i,:) +lu(j,:)
    j = j+1
  end do
  end subroutine 

  subroutine finer_y(u,lu,n,m,lm)
  implicit none
  integer,intent(in) :: n,lm,m
  real*8,dimension(n,lm),intent(in) :: lu
  real*8,dimension(n,m),intent(inout) :: u

  !local 
  integer :: i,j,k
  j = 2
  do i = 5, m-1, 4
    do k = 2,n-1
      u(k,i) = u(k,i) +lu(k,j)
    end do
    j = j+1
  end do
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

subroutine lev_pre(ax,bx,cc,ay,by,rinv,n,ln,m,lm)
  implicit none
  integer,intent(in) :: n,ln,m,lm! grid number of current and lower (coarser) levels
  integer,parameter :: nn = 5
  real*8,dimension(n,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(ln-1,lm-1,3,3),intent(inout) :: rinv

  !local 
  integer :: i,j,i4,j4,k,l,nn1
  k = 0
  nn1 = nn-1
  do i = 1, n-1, nn1
    k = k+1
    l = 0
    do j = 1, m-1, nn1
      l = l+1
      write(*,*) i,j,k,l
      i4 =i+nn1
      j4 =j+nn1
      call pre(ax(i:i4,j:j4),bx(i:i4,j:j4),cc(i:i4,j:j4),ay(i:i4,j:j4),by(i:i4,j:j4),rinv(k,l,:,:))
  end do
  end do

end subroutine 

! solve each block in one level
subroutine lev_evp_x(ax,bx,cc,ay,by,rinv,u,r,n,ln,mlev,m,tm,lm,tol)
  implicit none
  integer,intent(in) :: n,ln,mlev,m,tm,lm
  real*8,dimension(n,tm),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(n,m),intent(in) :: r
  integer, parameter :: nn = 3
  real*8,dimension(ln,tm-m+2,nn,nn),intent(in) :: rinv
  real*8,dimension(n,m),intent(inout) :: u
  real*8,intent(in) :: tol

  !local 
  integer :: i,j,i4

  j = 0
  
  do i = 1, n-1, nn+1
    j = j+1
    i4 = i+4
    !write(*,*) 'j=',j,'i=',i,'i+4',i+4
   call MG_y(ax(i:i4,:), bx(i:i4,:),cc(i:i4,:),ay(i:i4,:),by(i:i4,:),rinv(j,:,:,:),u(i:i4,:),r(i:i4,:),mlev,m,tm,lm,tol)
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
subroutine lev_evp_y(ax,bx,cc,ay,by,rinv,u,r,m,lm)
  implicit none
  integer,intent(in) :: m,lm
  integer, parameter :: nn =5
  real*8,dimension(nn,m),intent(in) :: ax,bx,cc,ay,by
  real*8,dimension(nn,m),intent(in) :: r
  real*8,dimension(lm,nn-2,nn-2),intent(in) :: rinv
  real*8,dimension(nn,m),intent(inout) :: u

  !local 
  integer :: i,j

  !write(*,*) 'in lev_evp_y'
  j = 0
  do i = 1, m-1, 4
    j = j+1
    !write(*,*) 'j=',j,'i=',i,'i+4',i+4
   call evp(ax(:,i:i+4), bx(:,i:i+4),cc(:,i:i+4),ay(:,i:i+4),by(:,i:i+4),rinv(j,:,:),u(:,i:i+4),r(:,i:i+4))

  end do

end subroutine 

! 5 points EVP 
subroutine evp(ax,bx,cc,ay,by,rinv,u,r)
  implicit none
  integer,parameter :: n = 5
  integer,parameter :: m = 5
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


subroutine grid_num(lev, lgrid,tgrid,lowgrid)
  implicit none
  integer,intent(in)  :: lev
  integer,intent(inout) :: lgrid,tgrid,lowgrid ! lev grid 
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





