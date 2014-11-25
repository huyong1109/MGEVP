!==============================================================
!   Domain Decomposition based Error Vector Propogation Method 
!==============================================================

program main
use domain, only: nx, ny, ax, ay, bb, cx, cy, rinv, rinv1, u, b, ie, nblk
use domain, only: init

    implicit none
    integer :: i, j , flag
    real(8)	:: tmp
    real(8), dimension(0:nx+1,0:ny+1)    :: x ! variable
    real(8), dimension(1:nx,1:ny)	 ::  f ! variable
    character(LEN=80)	:: fmt
    
    
    fmt = "(80f8.3)"

    

    call  init
    !write(*,*) nx, ny
    !write(*,*) xl, xr, cc, yb, yt
    do j = 1, ny
    do i = 1, nx

         f(i,j) = ax(i,j)*u(i-1,j) + cx(i,j)*u(i+1,j) + ay(i,j)*u(i,j-1) + cy(i,j)*u(i,j+1) +bb(i,j)*u(i,j)

    end do 
    end do 
    
    x(:,:) = u(:,:)
    x(1:nx,1:ny) =0.0d0

    !call rep(ax,ay,bb,cx,cy,rinv,rinv1,f,x,ie,nx, ny, nblk)
    
	    fmt = "(80f8.3)"
    	    !write(*,*) "coefficient ax "
    	    !write(*,fmt) ax
    	    !write(*,*) "coefficient ay "
    	    !write(*,fmt) ay
    	    !write(*,*) "coefficient bb "
    	    !write(*,fmt) bb
    	    !write(*,*) "coefficient cx "
    	    !write(*,fmt) cx
    	    !write(*,*) "coefficient cy "
    	    !write(*,fmt) cy
    	    !write(*,*) "rinv"
    	    !write(*,fmt) rinv(:,:,:), rinv1(:,:,:)
    	    !write(*,*) "fffg"
    	    !write(*,fmt) f
    	    !write(*,*) "ie "
    	    !write(*,*) ie(:)
    
   ! write(*,'(8f8.3)') rinv1(:,:,:)
    write(*,*) "Before ddevp"
    call ddevp(ax,ay,bb,cx,cy,rinv,rinv1,f,x,ie,nx, ny, nblk)
    write(*,*) "After ddevp"
    
    ! extending east-west boundary
    do j = 1, ny
        x(0,j) = x(1,j)
        x(nx+1,j) = x(nx,j)
    end do 

    write(*,*) 'check rep		err		    u'
    flag = 0
    do i = 1, nx
	do j = 1, ny
	    tmp = abs(x(i,j) - u(i, j))
	    !write(*,*) i, j, tmp, u(i,j)
	    if(tmp > 1.0D-8) then 
		write(*,*) i, j, tmp, u(i,j)
		flag = flag + 1
    	    endif 
	    end do 
    end do 

    if (flag == 0 ) then 
	print *, "Case pass "
    else
	print *, "Case fail "
    endif 

    call write_grads(x, nx+2, ny+2, 'x.dat', 12) 

end program
