subroutine printmatrix ( mat, n, m, iounit)
 
 use myconstants
 implicit none
 
 ! arguments
 real(dp), intent(in) :: mat(n,m)
 integer, intent(in) :: n, m, iounit

 ! local constants and variables
 integer, parameter :: w = 15, w2 = 6
 integer :: i, j
 
 write ( iounit, '( A, i12, A, i12, A )') "size: (", n, " * ", m, " )"

101 format( 3x, e10.3 )
102 format( "  ...  " )
103 format( 6x, "  ... ... ...  " )

 if ( n <= w .and. m <= w ) then
   ! print all the matrix 
   do i = 1, n
      do j = 1, m
        write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * ) 
   enddo
 endif

 if ( n <= w .and. m > w ) then
   ! print part of the matrix
   do i = 1, n
      do j = 1, w2
        write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, 102, advance='no' )
      do j = m-w2+1, m
        write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

 if ( n > w .and. m <= w ) then
   ! print part of the matrix
   do i = 1, w2
      do j = 1, m
        write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
   write ( iounit, 103 )
   do i = n-w2+1, n
      do j = 1, m
        write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

 if ( n > w .and. m > w ) then
   ! print part of the matrix
   do i = 1, w2
      do j = 1, w2
         write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, 102, advance='no' )
      do j = m-w2+1, m
         write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
   write ( iounit, 103 )
   do i = n-w2+1, n
      do j = 1, w2
         write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, 102, advance='no' )
      do j = m-w2+1, m
         write ( iounit, 101, advance='no' ) mat(i,j)
      enddo
      write ( iounit, * )
   enddo
 endif

end subroutine printmatrix
