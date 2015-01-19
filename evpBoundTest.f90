! 2D multigrid example from MUDPACK
!   solve f''(t) = -1, f(0)=f(1)=0
program main 
  implicit none
  real*8,parameter :: rtol = 1.0e-3
  integer, parameter :: n = 44
  integer, parameter :: m = 36
  integer, parameter :: nn =10
  integer, parameter :: mm =10
  integer, parameter :: miter = 200
  integer, parameter :: choice = 2
  integer, parameter :: solver = 1
  character*2 :: grd  = '00'
  real*8  :: tol,rr,mtmp,ntmp
  real*8,dimension(n,m) :: f,u0,u,r,tu
  real*8,dimension(n,m) :: cc,ns,ew,ne
  real*8,dimension(n,m) :: fcc,fns,few,fne,ff,fu0
  integer :: nb, mb  ! blocks on x and y direction 
  integer,dimension(:),allocatable :: ndi,mdi ! bound index 

  ! function 
  !integer :: grid_num

  ! local variable
  integer :: i,j,block,iter,cnt

  print *, '2D EVP preconditioned PCG solver'

 !  set problem 
 do block = 9,9
   !do i = 2,2
   write(GRD,'(I2.2)') block
   write(*,*) 'block ',block
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
   !========================= pop test case ====== 
   !!u0(10,:) = 1.0 
   !!u0(9+nn,:) = 1.0 
   !!u0(:,10) = 1.0 
   !!u0(:,9+mm) = 1.0 
   !!u0(:,:) = -u0(:,:)/100000.
   f(:,:) = 0.
   call matrixmultiple(cc,ns,ew,ne,u0,f,n,m)
   u(:,:) = u0
   u(2:n-1,2:m-1) = 0.
   !!call calc_rhs(tcc,tns,tew,tne,tu0,tf,f,nn,mm,rr)
   !!write(*,*) 'rr =', rr
   !!call calc_rhs(tcc,tns,tew,tne,tu,tf,f,nn,mm,rr)
   !!write(*,*) 'rr =', rr
   !=========================pop test case ====== 

   !=========================ideal test case ====== 
   !cc(:,:) = 1
   !ns(:,:) = -0.05
   !ew(:,:) = -0.05
   !ne(:,:) = -0.25
   !u0 = 1.0
   !f(:,:) = 0.
   !call matrixmultiple(cc,ns,ew,ne,u0,f,n,m)
   !write(*,*) 'f'
   !write(*,'(5f18.5)') f
   !u(:,:) = u0(1:n,1:m)
   !u(2:n-1,2:m-1) = 0.
   !write(*,*) 'inital u'
   !write(*,'(5e18.5)') u(1:n,1:m)
   !=========================ideal test case ====== 

   !call testninevp(tcc,tns,tew,tne,tu,tf,nn,mm)
   !call testevp(tcc,tns,tew,tne,tu,tf,nn,mm)

   write(*,*) 'tcc'
   write(*,'(5f15.5)') cc
   write(*,*) 'tns'
   write(*,'(5f15.5)') ns
   write(*,*) 'tew'
   write(*,'(5f15.5)') ew
   write(*,*) 'tne'
   write(*,'(5f15.5)') ne
   write(*,*) 'tu'
   write(*,'(5f15.5)') u
   write(*,*) 'tf'
   write(*,'(5f15.5)') f
   tol = rtol

   ! partition
   nb = (n-3)/nn + 1
   allocate(ndi(nb+1))
   ndi(1) = 2 
   if ( nb == 1 ) then 
     ndi(nb+1) = n
   else 
     do i = 1, nb-2
       ndi(i+1) = 2+i*nn
     end do 
     ndi(nb) = (ndi(nb-1) + n)/2 
     ndi(nb+1) = n
   end if 

   mb = (m-3)/mm + 1
   allocate(mdi(mb+1))
   mdi(1) = 2 
   if ( mb == 1 ) then 
     mdi(mb+1) = m
   else 
     do i = 1, mb-2
       mdi(i+1) = 2+i*mm
     end do 
     mdi(mb) = (mdi(mb-1) +m )/2 
     mdi(mb+1) = m
   end if 

   write(*,*) 'bound index --x '
   write(*,*) ndi(:)
   write(*,*) 'bound index --y '
   write(*,*) mdi(:)


   tu(:,:) = u(:,:) 
   call evppcg(cc,ns,ew,ne,tu,f,n,m,nn,mm,ndi,mdi,nb,mb,tol,miter)
   tu(:,:) = u(:,:) 
   call diagpcg(cc,ns,ew,ne,tu,f,n,m,tol,miter)
   tu(:,:) = u(:,:) 
   call pcg(cc,ns,ew,ne,tu,f,n,m,tol,miter)
   !call testevp(cc,ns,ew,ne,u,f,n,m)

   write(*,*) 'tu0'
   write(*,'(5e18.5)') u0(1:n,1:m)
   write(*,*) 'tu'
   write(*,'(5e18.5)') u(1:n,1:m)
   deallocate(ndi)
   deallocate(mdi)

 end do  
  end program

  subroutine testevp(cc,ns,ew,ne,u,r,n,m)
  integer,intent(in) :: n,m

  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,r
  real*8,dimension(n,m),intent(inout) :: u
  real*8,dimension(n-2,m-2) :: tu

  real*8,dimension(n+m-5,n+m-5) :: rinv
  write(*,*) 'r'
  write(*,'(5f18.5)') r
  write(*,*) 'pre'
  call exppre(cc,ns,ew,ne,rinv,n,m)
  write(*,*) 'evp'
  call expevp(cc,ns,ew,ne,rinv,u,tu,r,n,m)
  end subroutine

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


  subroutine evppre(cc,ns,ew,ne,rinv,landindx,n,m,nn,mm,ndi,mdi,nb,mb)
  implicit none
  integer,intent(in) :: n,m  ! total block size
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*8,dimension(nb*mb,nn+mm-1,nn+mm-1),intent(inout):: rinv
  integer,dimension(nb,mb),intent(inout):: landindx
  integer,intent(in) :: nn,mm  ! small block ideal size
  integer,intent(in) :: nb, mb  ! blocks on x and y direction 
  integer,intent(in) :: ndi(nb+1),mdi(mb+1) ! bound index 

  ! local 
  integer :: i,j,is,ie,js,je,ln,lm,l
  do i = 1, nb
    is = ndi(i)-1
    ie = ndi(i+1)
    ln = (ie-is) +1
    do j = 1, mb
      js = mdi(j)-1
      je = mdi(j+1) 
      if(minval(abs(ne(is+1:ie-1,js+1:je-1))) == 0.0 ) then 
        landindx(i,j) = 1
        rinv((i-1)*mb+j,:,:) = 0.
      else 
        landindx(i,j) = 0
        lm = (je-js) +1
        l  = ln + lm -5
        write(*,*) 'i : ', is, ie
        write(*,*) 'j : ', js, je
        write(*,*) 'l : ', l,'nn+mm-3',nn+mm-3
        call exppre(cc(is:ie,js:je),ns(is:ie,js:je),ew(is:ie,js:je),ne(is:ie,js:je),rinv((i-1)*mb+j,1:l,1:l),ln,lm)
      end if 
    end do 
  end do 
  write(*,*) 'total block ', nb*mb, ' land block ',sum(landindx)

  end subroutine 
  subroutine evpprecond(cc,ns,ew,ne,rinv,landindx,u,f,n,m,nn,mm,ndi,mdi,nb,mb)
  implicit none
  integer,intent(in) :: n,m  ! total block size
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,f
  real*8,dimension(n,m),intent(inout) :: u
  real*8,dimension(nb*mb,nn+mm-1,nn+mm-1),intent(inout):: rinv
  integer,dimension(nb,mb),intent(in):: landindx
  integer,intent(in) :: nn,mm  ! small block ideal size
  integer,intent(in) :: nb, mb  ! blocks on x and y direction 
  integer,intent(in) :: ndi(nb+1),mdi(mb+1) ! bound index 

  ! local 
  integer :: i,j,is,ie,js,je,ln,lm,l
  real*8,dimension(n,m) :: tu

  write(*,'(a30)') 'EVPPRECOND input'
  write(*,'(5e15.5)') f(:,:)
 
  tu = 0.0 
  do i = 1, nb
    is = ndi(i) -1
    ie = ndi(i+1) 
    ln = (ie-is) +1
    do j = 1, mb
      js = mdi(j) -1
      je = mdi(j+1) 
      lm = (je-js) +1
      l  = ln + lm -5
      write(*,*) 'i : ', is, ie,n
      write(*,*) 'j : ', js, je,m
      write(*,*) 'l : ', l,(i-1)*mb+j
      if (landindx(i,j) == 1 ) then 
        ! diagonal preconditioning for land blocks
        where(cc(is+1:ie-1,js+1:je-1) /= 0.0) 
          tu(is+1:ie-1,js+1:je-1) = & 
                    f(is+1:ie-1,js+1:je-1)/cc(is+1:ie-1,js+1:je-1)
        end where

      else 
        call expevp(cc(is:ie,js:je),ns(is:ie,js:je),ew(is:ie,js:je),ne(is:ie,js:je),rinv((i-1)*mb+j,1:l,1:l),u(is:ie,js:je),tu(is+1:ie-1,js+1:je-1),f(is:ie,js:je),ln,lm)
      end if 
    end do 
  end do 
  do i = 1, nb
    is = ndi(i) -1
    ie = ndi(i+1) 
    do j = 1, mb
      js = mdi(j)-1 
      je = mdi(j+1) 
      write(*,*) 'is+1,ie-1,js+1,je-1'
      write(*,*) is+1,ie-1,js+1,je-1
      u(is+1:ie-1,js+1:je-1) = tu(is+1:ie-1,js+1:je-1)
    end do 
  end do 
      write(*,'(a30)') 'EVPPRECOND output'
      write(*,'(5e15.5)') u(:,:)


  end subroutine 

subroutine exppre(cc,ns,ew,ne,rinv,n,m)
  implicit none
  integer,intent(in) :: n,m
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne

  real*8,dimension(n +m -5,n +m -5),intent(inout) :: rinv

  !local 
  integer :: i,j,k,ii,info
  integer :: nm 
  real*8,dimension(n,m) :: y
  real*8,dimension(n +m -5,n +m -5) :: work,ipiv
  real*8,dimension(n +m -5,n +m -5) :: rin
  real*8,dimension(n +m -5) :: r

  nm = n +m -5
  y(:,:) = 0.0
  ! left 
  do ii = 1,m-2
    y(2,m-ii) = 1.0
    i = 2
    j = 2
    ! left-bottom corner 
    y(i+1,j+1)  = -((cc(i,j)+ne(i-1,j-1)+ew(i-1,j)+ns(i,j-1))     * y(i,j )     & 
               +  (ns(i,j)+ne(i-1,j))     * y(i,j+1)    &
               +  (ew(i,j)+ne(i,j-1))     * y(i+1,j)    &
               +  ne(i-1,j-1) * y(i-1,j-1) ) /ne(i,j) 
    ! left boundary
    do j = 3, m-2
      y(i+1,j+1)  = -( (cc(i,j)+ew(i-1,j))     * y(i,j )     & 
      +  (ns(i,j)+ne(i-1,j))     * y(i,j+1)    &
      +  ns(i,j-1)   * y(i,j-1)    &
      +  ew(i,j)     * y(i+1,j)    &
      +  (ne(i,j-1)+ne(i-1,j-1))   * y(i+1,j-1)  &
      ) /ne(i,j) 
    end do
    ! bottom boundary
    do i = 2, n-1
      y(i+1,j+1)  = -((cc(i,j)+ns(i,j-1))     * y(i,j )     & 
      +  ns(i,j)                 * y(i,j+1)    &
      +  (ew(i,j)+ne(i,j-1))     * y(i+1,j)    &
      +  (ew(i-1,j)+ne(i-1,j-1)) * y(i-1,j)    &
      +  ne(i-1,j)               * y(i-1,j+1)  & 
      ) /ne(i,j) 
    end do

    do j = 3, m-2
      do i = 3, n-2
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
    ! left-bottom corner 
    rinv(ii,1)  = -((cc(i,j)+ns(i,j) +ew(i-1,j)+ne(i-1,j))     * y(i,j )     & 
               +  (ns(i,j-1) +ne(i-1,j-1))  * y(i,j-1)    &
               +  (ew(i,j)  +ne(i,j))   * y(i+1,j)    &
               +  ne(i,j-1)   * y(i+1,j-1) )
    ! top boundary
    do i = 3, n-2
      rinv(ii,i)  = -((cc(i,j)+ns(i,j))     * y(i,j )     & 
               +  ns(i,j-1)   * y(i,j-1)    &
               +  (ew(i,j)   +ne(i,j))      * y(i+1,j)    &
               +  (ew(i-1,j) +ne(i-1,j))    * y(i-1,j)    &
               +  ne(i,j-1)   * y(i+1,j-1)  &
               +  ne(i-1,j-1) * y(i-1,j-1) )
    end do
    ! right-top corner
    rinv(ii,1)  = -((cc(i,j)+ns(i,j) +ew(i-1,j)+ne(i-1,j))     * y(i,j )     & 
               +  (ns(i,j-1) +ne(i-1,j-1))  * y(i,j-1)    &
               +  (ew(i,j)  +ne(i,j))   * y(i+1,j)    &
               +  ne(i,j-1)   * y(i+1,j-1) )
    ! left boundary
    do j = 3, m-2
      y(i+1,j+1)  = -( (cc(i,j)+ew(i-1,j))     * y(i,j )     & 
      +  (ns(i,j)+ne(i-1,j))     * y(i,j+1)    &
      +  ns(i,j-1)   * y(i,j-1)    &
      +  ew(i,j)     * y(i+1,j)    &
      +  (ne(i,j-1)+ne(i-1,j-1))   * y(i+1,j-1)  &
      ) /ne(i,j) 
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

subroutine expevp(cc,ns,ew,ne,rinv,u,tu,f,n,m)
  implicit none
  integer:: n,m
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*8,dimension(n,m),intent(in) :: f
  real*8,dimension(n,m),intent(in) :: u
  real*8,dimension(n-2,m-2),intent(inout) :: tu

  real*8,dimension(n+m-5,n+m-5),intent(in) :: rinv

  !local 
  integer :: i,j,k,nm
  real*8,dimension(n,m) :: y,ry,ff
  real*8,dimension(n+m-5) :: r
  real*8 :: rr

  nm = n+m-5
  write(*,*) 'inital u'
  write(*,'(5f18.5)')  u(:,:)
  write(*,*) 'inital f'
  write(*,'(5f18.5)')  f(:,:)

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

  write(*,*) 'r'
  write(*,'(5f18.5)') r(:)
  write(*,*) 'rinv'
  write(*,'(5f18.5)') rinv(:,:)
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

  do j = 2, m-2
    do i = 2, n-2
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
  tu(1:n-2,1:m-2) = y(2:n-1,2:m-1) 
  !
  !u(:,:) = u(:,:)*cc(:,:)
  !write(*,*) 'final u'
  !write(*,'(5f18.5)')  u(:,:)
  !write(*,*) 'final f'
  !write(*,'(5f18.5)')  f(:,:)
  !y(:,:) = 0.
  ry(:,:) = u(:,:)
  ry(2:n-1,2:m-1) = tu
  y(:,:) = 0.0
  call calc_rhs(cc,ns,ew,ne,ry,f,y,n,m,rr)
  write(*,*) 'rr in evp  :', rr
  write(*,'(5e18.5)') y

end subroutine 


  subroutine calc_rhs(cc,ns,ew,ne,u,f,r,n,m,rr)
  implicit none
  integer,intent(in) :: n,m
  real*8, intent(inout) :: rr
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne
  real*8,dimension(n,m),intent(in) :: u,f
  real*8,dimension(n,m),intent(inout) :: r

  !local 
  integer :: i,j
  rr = 0.0
  r(:,:) = 0.0
  call matrixmultiple(cc,ns,ew,ne,u,r,n,m)
  r(:,:) = f(:,:) -r(:,:)
  !write(*,*) 'r'
  !write(*,'(5e18.5)') r
  do i = 2, n-1
    do j = 2, m-1
     rr = rr + r(i,j)*r(i,j)
    end do
  end do
end subroutine 

subroutine pcg(cc,ns,ew,ne,u,f,n,m,tol,miter)

implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,f
real*8,dimension(n,m),intent(inout) :: u
integer,intent(in):: miter
real*8,intent(in) :: tol
! local 
integer:: iter
real*8,dimension(n,m) :: p,ap,s,au,r
real*8 :: mu,nu,rho,alpha

p(:,:) = 0.
s(:,:) = 0.
r(:,:) = 0.
au(:,:) = 0.
ap(:,:) = 0.
call matrixmultiple(cc,ns,ew,ne,u,au,n,m)
r(:,:) = f(:,:) - au(:,:)
p(:,:) = r(:,:) 
call vectornorm(r,r,mu,n,m)

iter = 0
write(*,'(A5,I5.3,3e15.5)') 'CG ',iter,mu,tol
do while ((mu > tol) .and. (iter < miter))
  call matrixmultiple(cc,ns,ew,ne,p,ap,n,m)
  !write(*,*) ap
  call vectornorm(p,ap,rho,n,m)
  alpha = mu/rho
  if(rho ==0.)  exit
  write(*,*) 'mu,rho,alpha',mu,rho,alpha
  u(:,:) = u(:,:) + 1.*alpha*p(:,:)
  r(:,:) = r(:,:) - alpha*ap(:,:)
  call  vectornorm(r,r,nu,n,m)
  p(:,:) = r(:,:)+nu/mu*p(:,:)
  mu = nu
  iter = iter +1
  write(*,'(A5,I5.3,3e15.5)') 'CG ',iter,mu,tol
end do 
  write(*,'(A15,I5.3)') 'CG ITER ',iter
end
subroutine evppcg(cc,ns,ew,ne,u,f,n,m,nn,mm,ndi,mdi,nb,mb,tol,miter)
implicit none 
  integer,intent(in)  :: n,m
  integer,intent(in)  :: nn,mm
  integer,intent(in)  :: nb,mb
  integer,intent(in)  :: ndi(nb+1),mdi(mb+1) ! bound index 
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,f
  real*8,dimension(n,m),intent(inout) :: u
  integer,intent(in):: miter
  real*8,intent(in) :: tol
  ! local 
  real*8,dimension(nb*mb,nn+mm-1,nn+mm-1):: rinv
  integer,dimension(nb,mb) :: landindx
  integer:: iter
  real*8,dimension(n,m) :: p,ap,s,au,r,z
  real*8,dimension(n,m) :: sw
  real*8 :: mu,nu,rho,alpha,rr
  
  sw(:,:) = 0.
  call evppre(cc,sw,sw,ne,rinv,landindx,n,m,nn,mm,ndi,mdi,nb,mb)
  !u(2:n-1,2:m-1) = 0.
  
  !call evpprecond(cc,ns,ew,ne,rinv,landindx,u,f,n,m,nn,mm,ndi,mdi,nb,mb)
  p(:,:) = 0.
  s(:,:) = 0.
  r(:,:) = 0.
  z(:,:) = 0.
  au(:,:) = 0.
  ap(:,:) = 0.
  call matrixmultiple(cc,ns,ew,ne,u,au,n,m)
  r(:,:) = f(:,:) - au(:,:)
  !z(:,:) = u(:,:)!/cc(:,:) 
  z(2:n-1,2:m-1) = 0.
  !where(cc(:,:) /= 0.0) 
  !  z(:,:) = r(:,:)/cc(:,:)
  !else where 
  !  z(:,:) = 0.0
  !end where

  !write(*,*) 'initial r '
  !write(*,'(5e18.5)') r
  call evpprecond(cc,sw,sw,ne,rinv,landindx,z,r,n,m,nn,mm,ndi,mdi,nb,mb)
  !call matrixmultiple(cc,ns,ew,ne,z,au,n,m)
  !write(*,*) 'initial z '
  !write(*,'(5e18.5)') z(:,:)
  !write(*,*) 'initial az-r'
  !write(*,'(5e18.5)') au(:,:)-r(:,:)
  p(:,:) = z(:,:) 
  call vectornorm(r,z,mu,n,m)
  iter = 0
  call  vectornorm(r,r,rr,n,m)
  write(*,'(A15,I5.3,3e15.5)') 'EVPPCG ',iter,mu,rr,tol
  do while ((rr > tol) .and. (iter < miter))
    call matrixmultiple(cc,ns,ew,ne,p,ap,n,m)
    !write(*,*) ap
    call vectornorm(p,ap,rho,n,m)
    alpha = mu/rho
    write(*,*) 'mu,rho,alpha',mu,rho,alpha
    u(:,:) = u(:,:) + alpha*p(:,:)
    r(:,:) = r(:,:) - alpha*ap(:,:)
    !z(:,:) = u(:,:) !/cc(:,:)
    z(2:n-1,2:m-1) = 0.
    !where(cc(:,:) /= 0.0) 
    !  z(:,:) = r(:,:)/cc(:,:)
    !else where 
    !  z(:,:) = 0.0
    !end where
    call evpprecond(cc,sw,sw,ne,rinv,landindx,z,r,n,m,nn,mm,ndi,mdi,nb,mb)
    call  vectornorm(r,r,rr,n,m)
    !call lev_evp(ax,bx,cc,ay,by,rinv,z,r,n,nb,m,mb,nn)
    if(rr==0.)  exit
    call  vectornorm(z,r,nu,n,m)
    p(:,:) = z(:,:)+nu/mu*p(:,:)
    iter = iter +1
    mu = nu
    write(*,'(A15,I5.3,3e15.5)') 'EVPPCG ',iter,mu,rr,tol
  end do 
  write(*,'(A15,I5.3)') 'EVPPCG ITER ',iter
 end subroutine

subroutine diagpcg(cc,ns,ew,ne,u,f,n,m,tol,miter)

implicit none 
  integer :: n,m
  real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,f
  real*8,dimension(n,m),intent(inout) :: u
  integer,intent(in):: miter
  real*8,intent(in) :: tol
  ! local 
  integer:: iter
  real*8,dimension(n,m) :: p,ap,s,au,r,z
  real*8 :: mu,nu,rho,alpha,rr

  p(:,:) = 0.
  s(:,:) = 0.
  r(:,:) = 0.
  z(:,:) = 0.
  au(:,:) = 0.
  ap(:,:) = 0.
  call matrixmultiple(cc,ns,ew,ne,u,au,n,m)
  r(:,:) = f(:,:) - au(:,:)
  where(cc(:,:) /= 0.0 ) 
    z(:,:) = r(:,:)/cc(:,:) 
  end where
  p(:,:) = z(:,:) 
  call vectornorm(r,z,mu,n,m)
  iter = 0
  call  vectornorm(r,r,rr,n,m)
  write(*,'(A15,I5.3,3e15.5)') 'DIAGPCG ',iter,mu,rr,tol
  do while ((rr > tol) .and. (iter < miter))
    call matrixmultiple(cc,ns,ew,ne,p,ap,n,m)
    !write(*,*) ap
    call vectornorm(p,ap,rho,n,m)
    alpha = mu/rho
    write(*,*) 'mu,rho,alpha',mu,rho,alpha
    u(:,:) = u(:,:) + alpha*p(:,:)
    r(:,:) = r(:,:) - alpha*ap(:,:)
    where(cc(:,:) /= 0.0 ) 
      z(:,:) = r(:,:)/cc(:,:) 
    end where
    call  vectornorm(r,r,rr,n,m)
    !call lev_evp(ax,bx,cc,ay,by,rinv,z,r,n,nb,m,mb,nn)
    if(rr==0.)  exit
    call  vectornorm(z,r,nu,n,m)
    p(:,:) = z(:,:)+nu/mu*p(:,:)
    iter = iter +1
    mu = nu
    write(*,'(A15,I5.3,3e15.5)') 'DIAGPCG ',iter,mu,rr,tol
  end do 
  write(*,'(A15,I5.3)') 'DIAGPCG ITER ',iter
  end subroutine

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
subroutine matrixmultiple(cc,ns,ew,ne,u,s,n,m)
implicit none 
integer :: n,m
real*8,dimension(n,m),intent(in) :: cc,ns,ew,ne,u
real*8,dimension(n,m),intent(inout) :: s

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
real*8,dimension(n,m),intent(in) :: cc

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


subroutine lu0(a,b,c,r,n)

implicit none
  integer               :: i,n
  real*8                :: ai,sn,rn

  real*8,dimension(1:n) :: a,b,c,r
  real*8,dimension(1:n) :: s,t

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

end subroutine lu0
