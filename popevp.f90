! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*16,parameter :: rtol = 1.0e-9
  integer, parameter :: n = 44
  integer, parameter :: m = 36
  integer, parameter :: nn = 5
  integer, parameter :: mm = 5
  integer, parameter :: miter = 100
  integer, parameter :: choice = 2
  integer, parameter :: solver = 1
  character*2 :: grd  = '00'
  real*16  :: tol,rr,mtmp,ntmp
  real*16,dimension(n,m) :: f,u0,u,r
  real*16,dimension(n,m) :: cc,ns,ew,ne
  real*8,dimension(n,m) :: fcc,fns,few,fne,ff,fu0
  real*16,dimension(nn,mm) :: tcc,tns,tew,tne,tf,tu0,tu
  real*16,dimension(n-2,n-2) :: rinv

  ! function 
  !integer :: grid_num

  ! local variable
  integer :: i,j,iter,cnt

  print *, '2D EVP preconditioned PCG solver'

 !  set problem 
 ! do i = 0,95
  do i = 9,9
    write(GRD,'(I2.2)') i
    write(*,*) 'block ',i
    open(888,file='./popdata/ocnproc_'//GRD,form='unformatted')
    rewind(888)
    read(888) fcc
    cc = fcc
    read(888) fns
    ns = fns
    read(888) few
    ew = few
    read(888) fne
    ne = fne
    read(888) ff
    f = ff
    close(888)
    open(999,file='./popdata/ocnproc_x_'//GRD,form='unformatted')
    read(999) fu0
    u0 = fu0
    close(999)
    !write(*,*) 'cc'
    !write(*,'(4f18.5)') cc
    !call geoplot(cc,n,m,cnt)
    !write(*,*) 'ns'
    !write(*,'(4f18.5)') ns
    !call geoplot(ns,n,m,cnt)
    !write(*,*) 'ew'
    !write(*,'(4f18.5)') ew
    !call geoplot(ew,n,m,cnt)
    !write(*,*) 'ne'
    !write(*,'(4f18.5)') ne
    !call geoplot(ne,n,m,cnt)
    !write(*,*) 'f'
    !write(*,'(4f18.5)') f
    !call geoplot(f,n,m,cnt)
    !write(*,*) 'x'
    !write(*,'(4f18.5)') u0
    !call geoplot(u0,n,m,cnt)
  end do  
  !========================= pop test case ====== 
  !!u0(10,:) = 1.0 
  !!u0(9+nn,:) = 1.0 
  !!u0(:,10) = 1.0 
  !!u0(:,9+mm) = 1.0 
  !!u0(:,:) = -u0(:,:)/100000.
  !call matrixmultiple(cc,ns,ew,ne,u0,r,n,m)
  !do i = 1, nn
  !  do j = 1, mm
  !    tcc(i,j) = cc(0+i,j+0)
  !    tns(i,j) = ns(0+i,j+0)
  !    tew(i,j) = ew(0+i,j+0)
  !    tne(i,j) = ne(0+i,j+0)
  !    tu0(i,j) = u0(0+i,j+0)
  !    tf (i,j) = r (0+i,j+0)
  !  end do 
  !end do 
  !!tcc(:,:) = cc(1:nn,1:mm)
  !!tns(:,:) = ns(1:nn,1:mm)
  !!tew(:,:) = ew(1:nn,1:mm)
  !!tne(:,:) = ne(1:nn,1:mm)
  !!tu(:,:)  = u0(1:nn,1:mm)
  !!tf(:,:)  = r(1:nn,1:mm)
  !tu(:,:) = tu0
  !tu(2:nn-1,2:mm-1) = 0.
  !tf(1,:) = 0.
  !tf(nn,:) = 0.
  !tf(:,1) = 0.
  !tf(:,mm) = 0.
  !call calc_rhs(tcc,tns,tew,tne,tu0,tf,f,nn,mm,rr)
  !write(*,*) 'rr =', rr
  !call calc_rhs(tcc,tns,tew,tne,tu,tf,f,nn,mm,rr)
  !write(*,*) 'rr =', rr
  !=========================pop test case ====== 
  
  !=========================ideal test case ====== 
  tcc(:,:) = 1
  tns(:,:) = -0.05
  tew(:,:) = -0.05
  tne(:,:) = -0.25
  !call analytic_init(tu0,nn,mm,2)
  tu0 = 1.0
  tf(:,:) = 0.
  call matrixmultiple(tcc,tns,tew,tne,tu0,tf,nn,mm)
  write(*,*) 'tf'
  write(*,'(5f18.5)') tf
  tu(:,:) = tu0(1:nn,1:mm)
  tu(2:nn-1,2:mm-1) = 0.
  write(*,*) 'inital tu'
  write(*,'(5e18.5)') tu(1:nn,1:mm)
  !=========================ideal test case ====== 

  !call testninevp(tcc,tns,tew,tne,tu,tf,nn,mm)
  call testevp(tcc,tns,tew,tne,tu,tf,nn,mm)

  
  write(*,*) 'tcc'
  write(*,'(5f15.5)') tcc
  write(*,*) 'tns'
  write(*,'(5f15.5)') tns
  write(*,*) 'tew'
  write(*,'(5f15.5)') tew
  write(*,*) 'tne'
  write(*,'(5f15.5)') tne
  write(*,*) 'tu'
  write(*,'(5f15.5)') tu
  write(*,*) 'tf'
  write(*,'(5f15.5)') tf
  tol = rtol
  !call diagpcg(tcc,tns,tew,tne,tu,tf,nn,mm,tol,miter)
  !call pcg(tcc,tns,tew,tne,tu,tf,nn,mm,tol,miter)

  write(*,*) 'tu0'
  write(*,'(5e18.5)') tu0(1:nn,1:mm)
  write(*,*) 'tu'
  write(*,'(5e18.5)') tu(1:nn,1:mm)
  write(*,*) 'tu -tu0'
  write(*,'(5e18.5)') tu(1:nn,1:mm) -tu0(1:nn,1:mm)
  write(*,*) 'max(tu -tu0)'
  write(*,'(e18.5)') maxval(abs(tu(1:nn,1:mm) -tu0(1:nn,1:mm)))

  tu(:,:) = tu0
  tu(2:nn-1,2:mm-1) = 0.
  call testevp_diag(tcc,tns,tew,tne,tu,tf,nn,mm)
  write(*,*) 'tu0'
  write(*,'(5e18.5)') tu0(1:nn,1:mm)
  write(*,*) 'tu'
  write(*,'(5e18.5)') tu(1:nn,1:mm)
  write(*,*) 'tu -tu0'
  write(*,'(5e18.5)') tu(1:nn,1:mm) -tu0(1:nn,1:mm)
  write(*,*) 'max(tu -tu0)'
  write(*,'(e18.5)') maxval(abs(tu(1:nn,1:mm) -tu0(1:nn,1:mm)))


  !u(:,:) = 0.
  !call vectornorm(r,r,tol,n,m)
  !tol = tol*rtol
  !write(*,*) 'tol' ,tol
   
  !u(:,:) = u0(:,:)
  !u(2:n-1,2:m-1) = 0.
  !call pcg(cc,ns,ew,ne,u,r,n,m,tol,miter)
  !write(*,*) 'u'
  !write(*,'(4f18.5)') u



  end program

  subroutine testevp(cc,ns,ew,ne,u,r,n,m)
  integer,intent(in) :: n,m

  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,r
  real*16,dimension(n,m),intent(inout) :: u

  real*16,dimension(n+m-5,n+m-5) :: rinv
  write(*,*) 'r'
  write(*,'(5f18.5)') r
  write(*,*) 'pre'
  call exppre(cc,ns,ew,ne,rinv,n,m)
  write(*,*) 'evp'
  call expevp(cc,ns,ew,ne,rinv,u,r,n,m)


  end subroutine

  subroutine testevp_diag(cc,ns,ew,ne,u,r,n,m)
  integer,intent(in) :: n,m

  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,r
  real*16,dimension(n,m),intent(inout) :: u

  real*16,dimension(n,m) :: dcc,dns,dew,dne,dse,tr
  real*16,dimension(n+m-5,n+m-5) :: rinv

  ! test diagprocondition
  dcc(:,:) = cc(:,:)
  dns(:,:) = ns(:,:)
  dew(:,:) = ew(:,:)
  dne(:,:) = ne(:,:)
  dse(:,:) = 0.0
  call diagprecond(dcc,dns,dew,dne,dse,n,m)

  write(*,*) 'dcc'
  write(*,'(5f15.5)') dcc
  write(*,*) 'dns'
  write(*,'(5f15.5)') dns
  write(*,*) 'dew'
  write(*,'(5f15.5)') dew
  write(*,*) 'dne'
  write(*,'(5f15.5)') dne
  write(*,*) 'dse'
  write(*,'(5f15.5)') dse

  write(*,*) 'diag preconditioning'
  call exppre_diag(dcc,dns,dew,dne,dse,rinv,n,m)
  where(cc(:,:) /= 0.0 ) 
    u(:,:) = u(:,:)/sqrt(-cc(:,:))
    tr(:,:) = r(:,:)/sqrt(-cc(:,:))
  end where
  call expevp_diag(dcc,dns,dew,dne,dse,rinv,u,tr,n,m)
  where(cc(:,:) /= 0.0 ) 
    u(:,:) = u(:,:)*sqrt(-cc(:,:))
  end where
  

  end subroutine

  subroutine analytic_init(f,n,m,choice)
  implicit none
  integer,intent(in):: n,m
  integer,intent(in):: choice
  real*16, dimension(n,m), intent(inout) :: f

  ! local 
  integer :: i,j
  real*16 :: x,y,hm,hn
  real*16,parameter :: pi = 4.*atan(1.) 

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
  real*16, dimension(n,m), intent(in) :: u

  ! local 
  integer :: i,j
  real*16 :: hm,hn,x,y,tu,e,maxe
  real*16,parameter :: pi = 4.*atan(1.) 


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


  subroutine calc_rhs_diag(cc,ns,ew,ne,se,u,f,r,n,m,rr)
  implicit none
  integer,intent(in) :: n,m
  real*16, intent(inout) :: rr
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,se
  real*16,dimension(n,m),intent(in) :: u,f
  real*16,dimension(n,m),intent(inout) :: r

  !local 
  integer :: i,j
  rr = 0.0
  r(:,:) = 0.0
  call matrixmultiple_diag(cc,ns,ew,se,ne,u,r,n,m)
  write(*,*) r
  r(:,:) = f(:,:) -r(:,:)
  do i = 2, n-1
    do j = 2, m-1
     rr = rr + r(i,j)*r(i,j)
    end do
  end do
  end subroutine 

  subroutine calc_rhs(cc,ns,ew,ne,u,f,r,n,m,rr)
  implicit none
  integer,intent(in) :: n,m
  real*16, intent(inout) :: rr
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*16,dimension(n,m),intent(in) :: u,f
  real*16,dimension(n,m),intent(inout) :: r

  !local 
  integer :: i,j
  rr = 0.0
  r(:,:) = 0.0
  call matrixmultiple(cc,ns,ew,ne,u,r,n,m)
  write(*,*) r
  r(:,:) = f(:,:) -r(:,:)
  do i = 2, n-1
    do j = 2, m-1
     rr = rr + r(i,j)*r(i,j)
    end do
  end do
  end subroutine 

! 9 points PRE
subroutine pre(cc,ns,ew,ne,rinv,nn,mm)
  implicit none
  integer :: nn,mm
  real*16,dimension(nn,mm),intent(in) :: cc,ns,ew,ne

  real*16,dimension(nn-2,nn-2),intent(inout) :: rinv

  !local 
  integer :: i,j,k,ii,info
  real*16,dimension(nn,mm) :: y
  real*16,dimension(nn-2,nn-2) :: work,ipiv
  real*16,dimension(nn-2,nn-2) :: rin
  real*16,dimension(nn-2) :: r

  y(:,:) = 0.0
  do ii = 2, nn-1
    y(ii,2) = 1.0
    do j = 2, mm-2
      do i = 2, nn-1

        r(i-1) = -  cc(i,j)     * y(i,j )     & 
                 -  ns(i,j-1)   * y(i,j-1)    &
                 -  ew(i,j)     * y(i+1,j)    &
                 -  ew(i-1,j)   * y(i-1,j)    &
                 -  ne(i,j-1)   * y(i+1,j-1)  &
                 -  ne(i-1,j-1) * y(i-1,j-1) 
      end do
      call tridiagsolver(ns(:,j),ne(:,j),y(:,j+1),r(:),nn)
    end do
    ! get F error
    j = mm-1
    do i = 2, nn-1
      r(i-1) = -  cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  ne(i,j)     * y(i+1,j+1)  &
               -  ne(i,j-1)   * y(i+1,j-1)  &
               -  ne(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    rinv(ii-1,:) = -r(:)
    y(ii,2) = 0.0
  end do 

  write(*,*) 'rin'
  write(*,'(3e18.5)') rinv(:,:)
  rin(:,:) = rinv(:,:)
  work(:,:) = rinv(:,:)
  call inverse(rin,rinv,nn-2)


  !! check pre rinv 
  write(*,*) 'rinv'
  write(*,'(3f18.5)') rinv(:,:)
  rin(:,:) = 0.0
  do j = 1,nn-2
    do i = 1,nn-2
      do k = 1,nn-2
        rin(i,j) = rin(i,j) + rinv(i,k)*work(k,j)
      end do 
      if (i == j ) then 
        if (abs(rin(i,j) -1.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)-1.0
        endif 
      else
        if (abs(rin(i,j) -0.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)
        endif 
      endif 
    end do 
  end do 
  !write(*,*) 'work'
  !write(*,*) rin(:,:)

end subroutine 
subroutine exppre_diag(cc,ns,ew,ne,se,rinv,n,m)
  implicit none
  integer,intent(in) :: n,m
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,se

  real*16,dimension(n +m -5,n +m -5),intent(inout) :: rinv

  !local 
  integer :: i,j,k,ii,info
  integer :: nm 
  real*16,dimension(n,m) :: y
  real*16,dimension(n +m -5,n +m -5) :: work,ipiv
  real*16,dimension(n +m -5,n +m -5) :: rin
  real*16,dimension(n +m -5) :: r

  nm = n +m -5
  y(:,:) = 0.0
  ! left 
  do ii = 1,m-2
    y(2,m-ii) = 1.0
    do j = 2, m-1
      do i = 2, n-1
        y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  se(i,j-1)   * y(i+1,j-1)  &
               -  se(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
      end do
    end do
    ! get F error
    
    do i = 1,n-2
      rinv(ii,i) = -y(i+2,m)
    end do 

    do j = 1,m-3
      rinv(ii,n-2+j) = -y(n,m-j)
    end do 

    y(2,m-ii) = 0.0
  end do 
  ! buttom
  do ii = 1,n-3
    y(ii+2,2) = 1.0
    do j = 2, m-1
      do i = 2, n-1
        y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  se(i,j-1)   * y(i+1,j-1)  &
               -  se(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
      end do
    end do
    ! get F error
    
    do i = 1,n-2
      rinv(m-2+ii,i) = -y(i+2,m)
    end do 

    do j = 1,m-3
      rinv(m-2+ii,n-2+j) = -y(n,m-j)
    end do 

    y(ii+2,2) = 0.0
  end do 

  write(*,*) 'rin'
  write(*,'(5e18.5)') rinv(:,:)
  rin(:,:) = rinv(:,:)
  work(:,:) = rinv(:,:)
  call inverse(rin,rinv,nm)


  !! check pre rinv 
  write(*,*) 'rinv'
  write(*,'(5f18.5)') rinv(:,:)
  rin(:,:) = 0.0
  do j = 1,nm
    do i = 1,nm
      do k = 1,nm
        rin(i,j) = rin(i,j) + rinv(i,k)*work(k,j)
      end do 
      if (i == j ) then 
        if (abs(rin(i,j) -1.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)-1.0
        endif 
      else
        if (abs(rin(i,j) -0.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)
        endif 
      endif 
    end do 
  end do 
  !write(*,*) 'work'
  !write(*,*) rin(:,:)

end subroutine 
subroutine exppre(cc,ns,ew,ne,rinv,n,m)
  implicit none
  integer,intent(in) :: n,m
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne

  real*16,dimension(n +m -5,n +m -5),intent(inout) :: rinv

  !local 
  integer :: i,j,k,ii,info
  integer :: nm 
  real*16,dimension(n,m) :: y
  real*16,dimension(n +m -5,n +m -5) :: work,ipiv
  real*16,dimension(n +m -5,n +m -5) :: rin
  real*16,dimension(n +m -5) :: r

  nm = n +m -5
  y(:,:) = 0.0
  ! left 
  do ii = 1,m-2
    y(2,m-ii) = 1.0
    do j = 2, m-1
      do i = 2, n-1
        y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  ne(i,j-1)   * y(i+1,j-1)  &
               -  ne(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
      end do
    end do
    ! get F error
    
    do i = 1,n-2
      rinv(ii,i) = -y(i+2,m)
    end do 

    do j = 1,m-3
      rinv(ii,n-2+j) = -y(n,m-j)
    end do 

    y(2,m-ii) = 0.0
  end do 
  ! buttom
  do ii = 1,n-3
    y(ii+2,2) = 1.0
    do j = 2, m-1
      do i = 2, n-1
        y(i+1,j+1)  = (- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  ne(i,j-1)   * y(i+1,j-1)  &
               -  ne(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
      end do
    end do
    ! get F error
    
    do i = 1,n-2
      rinv(m-2+ii,i) = -y(i+2,m)
    end do 

    do j = 1,m-3
      rinv(m-2+ii,n-2+j) = -y(n,m-j)
    end do 

    y(ii+2,2) = 0.0
  end do 

  write(*,*) 'rin'
  write(*,'(5e18.5)') rinv(:,:)
  rin(:,:) = rinv(:,:)
  work(:,:) = rinv(:,:)
  call inverse(rin,rinv,nm)


  !! check pre rinv 
  write(*,*) 'rinv'
  write(*,'(5f18.5)') rinv(:,:)
  rin(:,:) = 0.0
  do j = 1,nm
    do i = 1,nm
      do k = 1,nm
        rin(i,j) = rin(i,j) + rinv(i,k)*work(k,j)
      end do 
      if (i == j ) then 
        if (abs(rin(i,j) -1.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)-1.0
        endif 
      else
        if (abs(rin(i,j) -0.0) > 1.0e-5 ) then 
          write(*,*) 'fail in pre',i,j,rin(i,j)
        endif 
      endif 
    end do 
  end do 
  !write(*,*) 'work'
  !write(*,*) rin(:,:)

end subroutine 
subroutine expevp_diag(cc,ns,ew,ne,se,rinv,u,f,n,m)
  implicit none
  integer:: n,m
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,se
  real*16,dimension(n,m),intent(in) :: f
  real*16,dimension(n,m),intent(inout) :: u

  real*16,dimension(n+m-5,n+m-5),intent(in) :: rinv

  !local 
  integer :: i,j,k,nm
  real*16,dimension(n,m) :: y,ry,ff
  real*16,dimension(n+m-5) :: r
  real*16 :: rr

  nm = n+m-5
  write(*,*) 'inital u'
  write(*,'(5f18.5)')  u(:,:)


  y(:,:) = u(:,:) 
  y(2:n-1,2) = y(2:n-1,1)
  y(2,3:m-1) = y(1,3:m-1)
  do j = 2, m-1
    do i = 2, n-1
        y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  se(i,j-1)   * y(i+1,j-1)  &
               -  se(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
    end do
  end do
  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)

  do i = 1,n-2
    r(i) = y(i+2,m)-u(i+2,m)
  end do 

  do j = 1,m-3
    r(n-2+j) = y(n,m-j) -u(n,m-j)
  end do 

  do j = 1,m-2
      do k = 1,nm
        y(2,m-j)  = y(2,m-j) + rinv(k,j)*r(k)
      end do 
  end do 
  do i = 1,n-3
      do k = 1,nm
        y(i+2,2)  = y(i+2,2) + rinv(k,m-2+i)*r(k)
      end do 
  end do 

  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)

  do j = 2, m-1
    do i = 2, n-1
        y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  se(i,j-1)   * y(i+1,j-1)  &
               -  se(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
    end do
  end do
  u(2:n-1,2:m-1) = y(2:n-1,2:m-1) 
  !
  !u(:,:) = u(:,:)*cc(:,:)
  write(*,*) 'final u'
  write(*,'(5f18.5)')  u(:,:)
  write(*,*) 'final f'
  write(*,'(5f18.5)')  f(:,:)
  y(:,:) = 0.
  call calc_rhs_diag(cc,ns,ew,ne,se,u,f,y,n,m,rr)
  write(*,*) 'rr in evp  :', rr
  write(*,'(5e18.5)') y

end subroutine 
subroutine expevp(cc,ns,ew,ne,rinv,u,f,n,m)
  implicit none
  integer:: n,m
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*16,dimension(n,m),intent(in) :: f
  real*16,dimension(n,m),intent(inout) :: u

  real*16,dimension(n+m-5,n+m-5),intent(in) :: rinv

  !local 
  integer :: i,j,k,nm
  real*16,dimension(n,m) :: y,ry,ff
  real*16,dimension(n+m-5) :: r
  real*16 :: rr

  nm = n+m-5
  write(*,*) 'inital u'
  write(*,'(5f18.5)')  u(:,:)



  y(:,:) = u(:,:) 
  y(2:n-1,2) = y(2:n-1,1)
  y(2,3:m-1) = y(1,3:m-1)
  do j = 2, m-1
    do i = 2, n-1
        y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  ne(i,j-1)   * y(i+1,j-1)  &
               -  ne(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
    end do
  end do
  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)

  do i = 1,n-2
    r(i) = y(i+2,m)-u(i+2,m)
  end do 

  do j = 1,m-3
    r(n-2+j) = y(n,m-j) -u(n,m-j)
  end do 

  do j = 1,m-2
      do k = 1,nm
        y(2,m-j)  = y(2,m-j) + rinv(k,j)*r(k)
      end do 
  end do 
  do i = 1,n-3
      do k = 1,nm
        y(i+2,2)  = y(i+2,2) + rinv(k,m-2+i)*r(k)
      end do 
  end do 

  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)

  do j = 2, m-1
    do i = 2, n-1
        y(i+1,j+1)  = (f(i,j)- cc(i,j)     * y(i,j )     & 
               -  ns(i,j)     * y(i,j+1)    &
               -  ns(i,j-1)   * y(i,j-1)    &
               -  ew(i,j)     * y(i+1,j)    &
               -  ew(i-1,j)   * y(i-1,j)    &
               -  ne(i,j-1)   * y(i+1,j-1)  &
               -  ne(i-1,j)   * y(i-1,j+1)  & 
               -  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
    end do
  end do
  u(2:n-1,2:m-1) = y(2:n-1,2:m-1) 
  !
  !u(:,:) = u(:,:)*cc(:,:)
  write(*,*) 'final u'
  write(*,'(5f18.5)')  u(:,:)
  write(*,*) 'final f'
  write(*,'(5f18.5)')  f(:,:)
  y(:,:) = 0.
  call calc_rhs(cc,ns,ew,ne,u,f,y,n,m,rr)
  write(*,*) 'rr in evp  :', rr
  write(*,'(5e18.5)') y

end subroutine 


! 9 points EVP 
subroutine evp(cc,ns,ew,ne,rinv,u,f,n,m)
  implicit none
  integer:: n,m
  real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*16,dimension(n,m),intent(in) :: f
  real*16,dimension(n,m),intent(inout) :: u

  real*16,dimension(n-2,n-2),intent(in) :: rinv

  !local 
  integer :: i,j,k
  real*16,dimension(n,m) :: y,ry,ff
  real*16,dimension(n-2) :: r
  real*16 :: rr
  write(*,*) 'inital u'
  write(*,'(5f18.5)')  u(:,:)

  u(:,:) = u(:,:) !/cc(:,:)
  ff(:,:) = f(:,:) !/cc(:,:)
  write(*,*) 'inital u1'
  write(*,'(5f18.5)')  u(:,:)
  write(*,*) 'rinv'
  write(*,'(3f18.5)')  rinv(:,:)


  y(:,:) = u(:,:) 
  y(2:n-1,2) = y(2:n-1,1)
  do j = 2, m-2
    do i = 2, n-1

      r(i-1) = ff(i,j)-cc(i,j)     * y(i,j )     & 
      -  ns(i,j-1)   * y(i,j-1)    &
      -  ew(i,j)     * y(i+1,j)    &
      -  ew(i-1,j)   * y(i-1,j)    &
      -  ne(i,j-1)   * y(i+1,j-1)  &
      -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    call tridiagsolver(ns(:,j),ne(:,j),y(:,j+1),r(:),n)
  end do
  j = m-1
    do i = 2, n-1
      r(i-1) = ff(i,j)-cc(i,j)     * y(i,j )     & 
           -  ns(i,j)     * y(i,j+1)    &
           -  ns(i,j-1)   * y(i,j-1)    &
           -  ew(i,j)     * y(i+1,j)    &
           -  ew(i-1,j)   * y(i-1,j)    &
           -  ne(i,j)     * y(i+1,j+1)  &
           -  ne(i,j-1)   * y(i+1,j-1)  &
           -  ne(i-1,j)   * y(i-1,j+1)  & 
           -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    write(*,*) 'r1 in j',j
    write(*,*) r(:)
  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)
  do i = 1,n-2 
      do k = 1,n-2
        y(i+1,2)  = y(i+1,2) + rinv(k,i)*r(k)
      end do 
  end do 
  write(*,*) 'y(:,2)'
  write(*,'(5f18.5)') y(:,:)
  do j = 2, m-2
    do i = 2, n-1

      r(i-1) = ff(i,j) -cc(i,j)     * y(i,j )     & 
      -  ns(i,j-1)   * y(i,j-1)    &
      -  ew(i,j)     * y(i+1,j)    &
      -  ew(i-1,j)   * y(i-1,j)    &
      -  ne(i,j-1)   * y(i+1,j-1)  &
      -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    call tridiagsolver(ns(:,j),ne(:,j),y(:,j+1),r(:),n)
    do i = 2, n-1
      r(i-1) = ff(i,j)-cc(i,j)     * y(i,j )     & 
           -  ns(i,j)     * y(i,j+1)    &
           -  ns(i,j-1)   * y(i,j-1)    &
           -  ew(i,j)     * y(i+1,j)    &
           -  ew(i-1,j)   * y(i-1,j)    &
           -  ne(i,j)     * y(i+1,j+1)  &
           -  ne(i,j-1)   * y(i+1,j-1)  &
           -  ne(i-1,j)   * y(i-1,j+1)  & 
           -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    write(*,*) 'r2 in j',j
    write(*,*) r(:)
  end do
  !last line 
  y(:,m) = u(:,m)
  j = m-1
    do i = 2, n-1
      r(i-1) = ff(i,j)-cc(i,j)     * y(i,j )     & 
           -  ns(i,j)     * y(i,j+1)    &
           -  ns(i,j-1)   * y(i,j-1)    &
           -  ew(i,j)     * y(i+1,j)    &
           -  ew(i-1,j)   * y(i-1,j)    &
           -  ne(i,j)     * y(i+1,j+1)  &
           -  ne(i,j-1)   * y(i+1,j-1)  &
           -  ne(i-1,j)   * y(i-1,j+1)  & 
           -  ne(i-1,j-1) * y(i-1,j-1) 
    end do
    write(*,*) 'r2 in j',j
    write(*,*) r(:)
  u(2:n-1,2:m-1) = y(2:n-1,2:m-1) 
  !
  !u(:,:) = u(:,:)*cc(:,:)
  write(*,*) 'final u'
  write(*,'(5f18.5)')  u(:,:)
  write(*,*) 'final f'
  write(*,'(5f18.5)')  f(:,:)
  y(:,:) = 0.
  call calc_rhs(cc,ns,ew,ne,u,f,y,n,m,rr)
  write(*,*) 'rr in evp  :', rr
  write(*,'(5e18.5)') y

end subroutine 

subroutine matrixmultiple_diag(cc,ns,ew,ne,se,u,s,n,m)
implicit none 
integer :: n,m
real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,se,u
real*16,dimension(n,m),intent(inout) :: s

! local 
integer :: i,j
s(:,:) = 0.
do i = 2, n-1 
  do j = 2, m-1
    s(i,j) =  cc(i,j)     * u(i,j )     & 
           +  ns(i,j)     * u(i,j+1)    &
           +  ns(i,j-1)   * u(i,j-1)    &
           +  ew(i,j)     * u(i+1,j)    &
           +  ew(i-1,j)   * u(i-1,j)    &
           +  ne(i,j)     * u(i+1,j+1)  &
           +  se(i,j-1)   * u(i+1,j-1)  &
           +  se(i-1,j)   * u(i-1,j+1)  & 
           +  ne(i-1,j-1) * u(i-1,j-1) 
  end do 
end do 
end subroutine 

subroutine matrixmultiple(cc,ns,ew,ne,u,s,n,m)
implicit none 
integer :: n,m
real*16,dimension(n,m),intent(in) :: cc,ns,ew,ne,u
real*16,dimension(n,m),intent(inout) :: s

! local 
integer :: i,j
s(:,:) = 0.
do i = 2, n-1 
  do j = 2, m-1
    s(i,j) =  cc(i,j)     * u(i,j )     & 
           +  ns(i,j)     * u(i,j+1)    &
           +  ns(i,j-1)   * u(i,j-1)    &
           +  ew(i,j)     * u(i+1,j)    &
           +  ew(i-1,j)   * u(i-1,j)    &
           +  ne(i,j)     * u(i+1,j+1)  &
           +  ne(i,j-1)   * u(i+1,j-1)  &
           +  ne(i-1,j)   * u(i-1,j+1)  & 
           +  ne(i-1,j-1) * u(i-1,j-1) 
  end do 
end do 
end subroutine 

subroutine geoplot(cc,n,m,cnt)
implicit none 
integer :: n,m,cnt
real*16,dimension(n,m),intent(in) :: cc

! local 
integer :: i,j
integer,dimension(n,m) :: geo
  where(cc == 0.0) 
    geo = 1
  elsewhere 
    geo = 0
  end where
  cnt = sum(geo)
  write(*,*) 'land ',cnt, 'ratio',dble(cnt)/(m*n)
  write(*,'(44I1.1)') geo
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
  real*16,intent(inout) :: a(n,n)
  real*16,intent(inout) :: c(n,n)
  real*16 :: L(n,n), U(n,n), b(n), d(n), x(n)
  real*16 :: coeff
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


SUBROUTINE LU0(a,b,c,r,n)

       implicit none

       integer               :: i,n
       real*16                :: ai,sn,rn

       real*16,dimension(1:n) :: a,b,c,r
       real*16,dimension(1:n) :: s,t

       s(1)=a(1)
       t(1)=c(n)

       sn=0.0
       rn=0.0

       do i=2,n-1
          ai=a(i)/b(i-1)
          b(i)=b(i)-ai*c(i-1)
          r(i)=r(i)-ai*r(i-1)
          s(i)=-ai*s(i-1)

          ai=t(i-1)/b(i-1)
          t(i)=-ai*c(i-1)
          sn  =sn-ai*s(i-1)
          rn  =rn-ai*r(i-1)
       enddo

       a(n)=a(n)+t(n-1)
       b(n)=b(n)+sn
       c(n-1)=c(n-1)+s(n-1)
       r(n)=r(n)+rn

       ai=a(n)/b(n-1)
       b(n)=b(n)-ai*c(n-1)
       r(n)=r(n)-ai*r(n-1)

       r(n)=r(n)/b(n)
       r(n-1)=(r(n-1)-c(n-1)*r(n))/b(n-1)

       do i=n-2,1,-1
          ai=r(i)-s(i)*r(n)-c(i)*r(i+1)
          r(i)=ai/b(i)
       enddo

       return

END SUBROUTINE LU0
subroutine tridiagsolver(ns,ne,x,r,n)
implicit none 
  integer, intent(in) :: n
  real*16, dimension(n), intent(in) :: ns,ne
  real*16, dimension(n-2), intent(in) :: r
  real*16, dimension(n), intent(inout) :: x
 
  !local
  integer :: i  
  real*16, dimension(n-2) :: tr
  real*16, dimension(n-2) :: a,b,c
  tr(:)  = r(:)
  tr(1) = tr(1) - ne(1)*x(1)
  tr(n-2) = tr(n-2) - ne(n-1)*x(n)
  a(:) = 0.0
  b(:) = 0.0
  c(:) = 0.0
  a(2:n-2) = ne(2:n-2)
  b(:) = ns(2:n-1)
  c(1:n-3) = ne(2:n-2)
  write(*,*) 'a',a
  write(*,*) 'b',b
  write(*,*) 'c',c
  call LU0(a,b,c,tr,n-2)
  x(2:n-1) = tr(:)
  write(*,*) 'x',x
  write(*,*) 'r',r
  ! check 
  do i = 2, n-1
    tr(i-1) = r(i-1) -(ne(i-1)*x(i-1) + ns(i) *x(i)+ ne(i)*x(i+1))
    if (abs(tr(i-1))  > 1.0e-8 ) then
      write(*,*) 'tri-solver fail',i, tr(i-1)
    end if 
  end do 
end subroutine 
subroutine diagprecond(cc,ns,ew,ne,se,n,m)
implicit none 
  integer, intent(in) :: n,m
  real*16, dimension(n,m),intent(inout) :: cc,ns,ew,ne
  real*16, dimension(n,m),intent(inout) :: se
  
  integer :: i,j
  real*16 :: maxcc,maxns,maxew,maxne

  maxcc  = maxval(-cc)
  maxns  = maxval(-ns)
  maxew  = maxval(-ew)
  maxne  = maxval(-ne)
  
  where(cc(:,:) == 0.0) 
    cc = maxcc
  end where
  where(ns(:,:) == 0.0) 
    ns = maxns
  end where
  where(ew(:,:) == 0.0) 
    ew = maxew
  end where
  where(cc(:,:) == 0.0) 
    ne = maxne
  end where
  
  se(:,:) = ne(:,:)
  do i = 1,n-1
    do j = 1,m-1
      ns(i,j) = ns(i,j) /sqrt(cc(i,j) *cc(i,j+1))
      ew(i,j) = ew(i,j) /sqrt(cc(i,j) *cc(i+1,j))
      ne(i,j) = ne(i,j) /sqrt(cc(i,j) *cc(i+1,j+1))
      se(i,j) = se(i,j) /sqrt(cc(i,j+1) *cc(i+1,j+1))
    end do 
  end do 
  cc(:,:) = -1.0

end subroutine
