!
SUBROUTINE euler(dt,iter)

   use module_para
   use module_array
   use module_mpi
   use module_io
!$   use omp_lib
   implicit none 
   integer i,i1,i2,j,k,iter,kni,bk_cir
   integer::stat(MPI_STATUS_SIZE)

   real*8,dimension(1:n1,1:n)  :: tu,tv,th,h 
   real*8,dimension(1:n1,1:n)  :: du,dv,dh ,rh
   real*8,dimension(1:n1,1:n)      :: a1,a2
   real*8,dimension(1:n1)	:: ru 
   real*8,dimension(1:kn)      :: fm,fp,f0,gm,gp,g0,rf,rg 
   real*8,dimension(1:kn*core_x_num) :: mfm,mfp,mf0,mgm,mgp,mg0,mrf,mrg 

   real*8                      :: ai,aj,dt,dt2,en,en0,den

   real*8,external             :: inner
   real*8			    :: sum1,sum2
   integer	:: ierr

!@   call omp_set_num_threads(numthreads)
   dt2=dt*0.5d0

   !$OMP PARALLEL DO PRIVATE(j,i)
   do j=1,n
    do i=1,n1
     tu(i,j)=wu(i,j)
     tv(i,j)=wv(i,j)
     th(i,j)=wh(i,j)
    end do
   end do
   !$OMP END PARALLEL DO    
       
   en0=inner(wu,wv,wh,wu,wv,wh)
   call MPI_Reduce(en0,sum1,1,MPI_DBL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   
   do k=1,1000

   !$OMP PARALLEL DO PRIVATE(j,i)
      do j=1,n
         do i=1,n1
            h(i,j)=dsqrt(th(i,j))
         end do
      end do
   !$OMP END PARALLEL DO    
 
   !print *, k
   call DIFUH(tu,tv,du,dh,h)
   call bound_trans(du)
   call bound_trans(dh)
  
   !$OMP PARALLEL DO PRIVATE(j,i)
   do j=1,n
      do i=1,n1
         tu(i,j)=wu(i,j)-dt2*du(i,j)
         th(i,j)=wh(i,j)-dt2*dh(i,j)
      end do
   end do
   !$OMP END PARALLEL DO


   !$OMP PARALLEL DO PRIVATE(j,i,ai,ru)
   do j=s_y,net_y+1
        ai=dt2*c11(j)

         do i=1,n1
	     a1(i,j)=ai*h(i,j)
	     ru(i)=tu(i,j)*ai*h(i,j)
	     a2(i,j)=ai*h(i,j)*ai*h(i,j)
	 end do

       	do i=2,n_x
       	     rh(i,j)=th(i,j)-ru(i+1)+ru(i-1)
        enddo
   enddo
   !$OMP END PARALLEL DO

   call bound_trans(rh)

   !!$OMP PARALLEL DO PRIVATE(j,i,i1,i2,kni,fm,fp,f0,gm,gp,g0,rf,rg,mfm,mfp,mf0,mgm,mgp,mg0,mrf,mrg)
   do j=s_y,net_y+1
       	do i=1,kn
       		 i1=i*2
       		 i2=i1+1

       		 fp(i)=-a2(i2,j)
       		 rf(i)=rh(i1,j)

       		 gm(i)=-a2(i1,j)
       		 rg(i)=rh(i2,j)

       		 mfp(i)=fp(i)
       		 mgm(i)=gm(i)
       		 mrf(i)=rf(i)
       		 mrg(i)=rg(i)
       	  enddo

       if(core_i /= 0) then 
           call MPI_Send(fp,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_RECV(mfp(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,stat,ierr)
           enddo
       endif


       if(core_i /= 0) then 
           call MPI_Send(gm,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_RECV(mgm(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,stat,ierr)
           enddo
       endif
       

       if(core_i /= 0) then 
           call MPI_Send(rf,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_RECV(mrf(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,stat,ierr)
           enddo
       endif


       if(core_i /= 0) then 
           call MPI_Send(rg,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_RECV(mrg(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,stat,ierr)
           enddo
       endif


       if( core_i .eq. 0) then
       	  kni = kn*core_x_num
       	  do i=2,kni
       	     mfm(i)=mfp(i-1)
       	  enddo
       	  mfm(1)=mfp(kni)
       	  
       	  do i=1,kni-1
       		 mgp(i)=mgm(i+1)
       	  enddo
       	  mgp(kni)=mgm(1)

       	  do i=1,kni
       		mf0(i)=1.0-mfm(i)-mfp(i)
       		mg0(i)=1.0-mgm(i)-mgp(i)
       	  enddo
       	  
       	  call LU0(mgm,mg0,mgp,mrg,kni)
       	  call LU0(mfm,mf0,mfp,mrf,kni)
       
       	    do i = 1,kn
       		 fp(i)=mfp(i)
       		 gm(i)=mgm(i)
       		 rf(i)=mrf(i)
       		 rg(i)=mrg(i)
       	    enddo
       endif


       if(core_i /= 0) then 
           call MPI_RECV(rf,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,stat,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_SEND(mrf(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,ierr)
           enddo
       endif


       if(core_i /= 0) then 
           call MPI_RECV(rg,kn,MPI_DBL,rank-MOD(rank,core_x_num),rank,MPI_COMM_WORLD,stat,ierr)
       else 
           do i = 1,core_x_num-1
       	call MPI_SEND(mrg(kn*i+1),kn,MPI_DBL,rank+i,rank+i,MPI_COMM_WORLD,ierr)
           enddo
       endif

       do i=1,kn
            i1=i*2
            i2=i1+1
	    th(i1,j)=rf(i)
            th(i2,j)=rg(i)
       enddo
   enddo
   !!$OMP END PARALLEL DO

    call bound_trans(th)

   !$OMP PARALLEL DO PRIVATE(i,j)
    do j=s_y,net_y+1
	do i=2,n_x
            tu(i,j)=tu(i,j)-dt2*c11(j)*h(i,j)*(th(i+1,j)-th(i-1,j))
       enddo
    enddo
   !$OMP END PARALLEL DO
   
    call bound_trans(tu)
    call DIFV(tu,tv,th,dv,h)
    call bound_trans(dv)
    
   !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,n
        do i=1,n1
           tv(i,j)=wv(i,j)-dt2*dv(i,j)
        end do
    end do
   !$OMP END PARALLEL DO

    en=inner(tu,tv,th,tu,tv,th)
    bk_cir = 0 
    call MPI_Reduce(en,sum2,1,MPI_DBL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   
    if(rank == 0) then
  	den=dabs(sum2-sum1)*2.0/(sum2+sum1)
  	!write(*,*),k,sum1,sum2
  	sum1=sum2
  	if (den.lt.1d-15) bk_cir = 1
    endif
    
    call MPI_Bcast(bk_cir,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
    if(bk_cir .eq. 1)  goto 10 
  
 end do

10   continue
   
       iter=k

   !$OMP PARALLEL DO PRIVATE(i,j)
       do j=1,n
       do i=1,n1
        wu(i,j)=tu(i,j)*2.0d0-wu(i,j)
        wv(i,j)=tv(i,j)*2.0d0-wv(i,j)
        wh(i,j)=th(i,j)*2.0d0-wh(i,j)
       end do
       end do
   !$OMP END PARALLEL DO
       return
END SUBROUTINE euler

SUBROUTINE LU0(a,b,c,r,n)

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

END SUBROUTINE LU0

