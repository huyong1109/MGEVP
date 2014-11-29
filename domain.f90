module domain
  implicit none 
  integer, parameter      :: nx = 8    ! grid numbers along x-direction
  !integer, parameter      :: ny = 80     ! grid numbers along y-direction
  integer, parameter      :: ny = 8     ! grid numbers along y-direction
  real(8), parameter      :: r1 = 1.0d0 ! double 1.
  real(8), parameter      :: r0 = 0.0d0 ! double 0.
  real(8), parameter      :: PI = 4.d0*datan(r1) !pi
  real(8), parameter      :: lx = PI    ! domain size along x-direction
  real(8), parameter      :: ly = r1    ! domain size along y-direction
  real(8), parameter      :: hx = lx/nx ! grid size along x-direction
  real(8), parameter      :: hy = ly/ny ! grid size along y-direction
  real(8), dimension(1:nx,1:ny)    :: ax, ay, bb, cx, cy ! A coefficient
  real(8), dimension(0:nx+1,0:ny+1)    :: u ! variable
  real(8), dimension(0:nx+1,0:ny+1)    :: b ! right hand value


  ! evp variables
  integer, parameter                            :: nblock = 4
  integer, parameter                            :: nblk = ny/nblock
  integer, dimension(0:nblk)                    :: ie

  real(8), dimension(nx, nx, nblk) :: rinv,rinv1


contains

 subroutine init()
    implicit none
   
    ! local variables 
    integer     :: i, j, n
    real(8)     :: x, y
    

    ax(:,:) = r0
    ay(:,:) = r0
    bb(:,:) = r0
    cx(:,:) = r0
    cy(:,:) = r0
    u (:,:) = r0
    b (:,:) = r0

    do j = 1, ny
        y = hy*j
        do i = 1, nx
        ax(i,j) = r1
        ay(i,j) = r1 
        cx(i,j) = r1 
        cy(i,j) = r1 
        bb(i,j) = -ax(i,j) - ay(i,j) - cx(i,j) - cy(i,j)
        x = hx*i
        u(i,j) = r1/(9.0D0+PI**2)*cos(3.0d0*x)*sin(PI*y)
        !u(i,j) = r1 ! r1/(9.0D0+PI**2)*cos(3.0d0*x)*sin(PI*y)
        b(i,j) = -hx*hy*cos(3.0D0*x)*sin(PI*y)
        !b(i,j) = 0.0 !-hx*hy*cos(3.0D0*x)*sin(PI*y)
        end do 
    end do 

    ! extending east-west boundary
    !do j = 1, ny
    !    !y = r1*j/ny 
    !    bb(1,j) =  bb(1,j) +ax(1,j)
    !    ax(1,j) =  0.
    !    u(0,j) = u(1,j)
    !    !cc(0,j) =  - xr(i,j) 
    !    !xl(nx,j) = r1 
    !    !cc(nx,j) = -xl(i,j) 
    !    bb(nx,j) =  bb(nx,j) +cx(nx,j)
    !    cx(nx,j) =  0.
    !    u(nx+1,j) = u(nx,j)
    !end do 
    !do i = 1, nx
        !y = r1*j/ny 
        !bb(i,1) =  bb(i,1) +ay(i,1)
        !ay(i,1) =  0.
        !u(i,1) = 0.
        !cc(0,j) =  - xr(i,j) 
        !xl(nx,j) = r1 
        !cc(nx,j) = -xl(i,j) 
        !bb(i,ny) =  bb(i, ny) +cy(i, ny)
        !cy(i, ny) =  0.
        !u(nx+1,j) =  0.
    !end do 

    call write_grads(u, nx+2, ny+2, 'u.dat', 11) 

    call init_test(nx, ny, ax, ay, bb, cx, cy, b, u)
   
    ! ie start from 0, in consideration of convinence in DDEVP
    do n = 0, nblk
        ie(n) = n*nblock+1
    end do 
    
    write(*,*) "===================================="
    write(*,*) "domain info : "
    write(*,*) nx, ny, PI, hx, hy
    write(*,*) ie
    write(*,*) "===================================="

    call pre(ax, ay, bb, cx, cy, rinv, rinv1,nx, ny, ie, nblk)
    call check_nan(ax,nx, ny, 'ax')
    call check_nan(ay,nx, ny, 'ay')
    call check_nan(bb,nx, ny, 'bb')
    call check_nan(cx,nx, ny, 'cx')
    call check_nan(cy,nx, ny, 'cy')
    call check_nan(rinv, nx, nx, 'rinv')
    call check_nan(rinv1, nx, nx, 'rinv1')


    write(*,*) "===================================="
    write(*,*) "Preprocessing end "
    write(*,*) "===================================="

 end subroutine init

 subroutine init_test(nx, ny, ax, ay, bb, cx, cy, b, u)
   implicit none 
   integer, intent(in)          :: nx, ny
   real(8), dimension(1:nx, 1:ny),intent(in)       ::  ax, ay, bb, cx, cy 
   real(8), dimension(0:nx+1, 0:ny+1), intent(in)           ::  b, u
   ! local variable
   integer      :: i, j, flag
   real(8)      :: tmp
  
   flag = 0 
   do j = 1, ny
   do i = 1, nx

        tmp = b(i,j)  -ax(i,j)*u(i-1,j) - cx(i,j)*u(i+1,j) - ay(i,j)*u(i,j-1) - cy(i,j)*u(i,j+1) -bb(i,j)*u(i,j)

        if(abs(tmp )> 1.0E-3) then
            flag  = flag + 1
            print *, '(', i, ',', j, ')'
            write(*,'(6f8.4)') b(i,j) , ax(i,j), ay(i,j), cx(i,j),  cy(i,j), bb(i,j)
            write(*,'(6f8.4)') tmp, u(i-1,j), u(i+1,j), u(i,j-1), u(i,j+1), u(i,j)
        endif 
   end do 
   end do 
   
   write(*,*) "===================================="
   if(flag .eq. 0) then
       print *, "Initial value test pass !"
   else
       print *, "Initial value test fail !"
   endif 
   write(*,*) "===================================="



 end subroutine init_test
end module domain
