module mod_csr_matrix

  implicit none

!*******************************************************************!
! define csr_matrix type                                            !
!*******************************************************************!
  type csr_matrix

  ! dimension of the matrix
  integer n
  ! position of first value of each row in col_ind and val
  integer, pointer, dimension(:) :: row_ptr
  ! total number of non zero coefficients, size of col_ind and val
  integer nz
  ! column number of each non zero value
  integer, pointer, dimension(:) :: col_ind
  ! value of non zero coefficients
  real(kind=kind(0.d0)), allocatable, dimension(:) :: val

  end type csr_matrix

contains
!*******************************************************************!
! create csr matrix                                                 !
!*******************************************************************!
  subroutine create_csr_matrix(n, nz, a)
    integer :: n, nz
    type(csr_matrix) :: a

    a%n = n
    a%nz = nz
    allocate(a%row_ptr(n+1))
    allocate(a%val(nz))
    allocate(a%col_ind(nz))
  end subroutine create_csr_matrix

!*******************************************************************!
! destroy csr matrix                                                !
!*******************************************************************!
  subroutine destroy_csr_matrix(a)
    type(csr_matrix) :: a

    deallocate(a%val)
    deallocate(a%col_ind)
    deallocate(a%row_ptr)
  end subroutine destroy_csr_matrix
  
!*******************************************************************!
! compute matrix vector product                                     !
!*******************************************************************!
  subroutine csr_mat_vec_prod(a, x, y)
    type(csr_matrix) :: a
    real(kind=kind(0.d0)), dimension(a%n) :: x, y
    integer :: irow, icol

    if (size(x) /= a%n) then
      print *, 'csr_mat_vec_prod: incompatible size'
      stop
    end if
 
    do irow = 1, a%n
      y(irow) = 0.d0
      do icol = a%row_ptr(irow), a%row_ptr(irow+1)-1
        y(irow) = y(irow) + a%val(icol)*x(a%col_ind(icol))
      end do
    end do

  end subroutine csr_mat_vec_prod

!*******************************************************************!
! print matrix                                                      !
!*******************************************************************!
  subroutine print_csr_matrix(a)
    type(csr_matrix) :: a
    integer :: irow, iptr
    
    do irow = 1, a%n
       print *, 'line :', irow
       do iptr = a%row_ptr(irow), A%row_ptr(irow+1)-1
          print*, '  col :', a%col_ind(iptr), ':', a%val(iptr)
       end do
    end do
  end subroutine print_csr_matrix

end module mod_csr_matrix
