subroutine write_grads(array, lx, ly, fn, nfn)

        integer, intent(in)            :: lx, ly, nfn ! array dimension
        character(*), intent(in)       :: fn    ! file index name
        real(8), dimension(lx,ly), intent(in)     :: array

        ! local variable 

        integer           :: i, j
        real(4), dimension(lx,ly)     :: array_io

array_io(:,:) = array(:,:)


        open(unit= nfn, file = fn, form='unformatted', access='direct', recl=lx*ly)
        !open(unit= nfn, file = fn, form='binary')
        write(nfn, rec=1) ((array_io(i,j), i=1,lx), j=1,ly)
close(nfn)

        end subroutine write_grads

subroutine check_nan(mat, nx, ny, str)
        integer :: i,j, flag
real*8 :: mat(1:nx,1:ny)
        character (*) :: str	
        flag = 0
        do j = 1, nx
        do i = 1, ny
        if (abs(mat(i,j)) > 1.0E+10) then 
        flag = flag+1
        endif 
        enddo
        enddo

        if (flag > 0) then 
        write(*,*) "NaN detacted in ", str
        endif
        end subroutine
