module mod_bz_2eq_1d

use mod_precision
use mod_cartesian_grid
use mod_csr_matrix

integer, private, parameter :: neq = 2
real(kind=dp), private, parameter :: eps = 1.d-02
real(kind=dp), private, parameter :: f = 3.d0
real(kind=dp), private, parameter :: q = 2.d-04
real(kind=dp), private, parameter :: db = 1.d0
real(kind=dp), private, parameter :: dc = 6.d-01
real(kind=dp), private :: dboverdxdx
real(kind=dp), private :: dcoverdxdx
real(kind=dp), private :: oneovereps
integer, private :: nx

contains

subroutine init_bz_2eq(grid)
  implicit none
  type(cartesian_grid_type) :: grid
  real(kind=dp) :: dx

  dx = grid%dx

  dboverdxdx = db / (dx*dx)
  dcoverdxdx = dc / (dx*dx)

  oneovereps = 1.d0 / eps

  nx = grid%nx

end subroutine init_bz_2eq

function f_bz_2eq(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq
  ! local 
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  bi   = u(1)
  ci   = u(2)
  bip1 = u(3)
  cip1 = u(4)

  f_bz_2eq(1) = (dboverdxdx * (bip1 - bi)) + oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_bz_2eq(2) = (dcoverdxdx * (cip1 - ci)) + (bi - ci)

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    bim1 = u(irow-2)
    cim1 = u(irow-1)
    bi   = u(irow)
    ci   = u(irow+1)
    bip1 = u(irow+2)
    cip1 = u(irow+3)

    f_bz_2eq(irow)   = (dboverdxdx * (bim1 - 2.d0*bi + bip1)) + oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
    f_bz_2eq(irow+1) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1)) + (bi - ci)

  end do 

  ! right boundary
  bim1 = u(n-3)
  cim1 = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_2eq(n-1) = (dboverdxdx * (bim1 - bi)) +  oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_bz_2eq(n)   = (dcoverdxdx * (cim1 - ci)) + (bi - ci)

end function f_bz_2eq

function f_loc_bz_2eq_reac(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_loc_bz_2eq_reac
  ! local
  real(kind=dp) :: bi
  real(kind=dp) :: ci

  bi   = u(1)
  ci   = u(2)

  f_loc_bz_2eq_reac(1) = oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_loc_bz_2eq_reac(2) = (bi - ci)

end function f_loc_bz_2eq_reac

function f_bz_2eq_reac(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq_reac
  ! local 
  real(kind=dp) :: bi
  real(kind=dp) :: ci
  integer :: inx, irow

  do inx = 1, nx

    irow = 1 + (inx-1)*neq

    bi   = u(irow)
    ci   = u(irow+1)

    f_bz_2eq_reac(irow)   = oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
    f_bz_2eq_reac(irow+1) = (bi - ci)

  end do 

end function f_bz_2eq_reac

function f_bz_2eq_diff(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq_diff
  ! local 
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  bi   = u(1)
  ci   = u(2)
  bip1 = u(3)
  cip1 = u(4)

  f_bz_2eq_diff(1) = (dboverdxdx * (bip1 - bi))
  f_bz_2eq_diff(2) = (dcoverdxdx * (cip1 - ci))

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    bim1 = u(irow-2)
    cim1 = u(irow-1)
    bi   = u(irow)
    ci   = u(irow+1)
    bip1 = u(irow+2)
    cip1 = u(irow+3)

    f_bz_2eq_diff(irow)   = (dboverdxdx * (bim1 - 2.d0*bi + bip1))
    f_bz_2eq_diff(irow+1) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1))

  end do 

  ! right boundary
  bim1 = u(n-3)
  cim1 = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_2eq_diff(n-1) = (dboverdxdx * (bim1 - bi))
  f_bz_2eq_diff(n)   = (dcoverdxdx * (cim1 - ci))

end function f_bz_2eq_diff

subroutine build_bz_2eq_diff_matrix(n, a)
  implicit none
  integer :: n
  type(csr_matrix) :: a
  ! local variables
  integer :: nz
  integer :: inx, irow, inz

  nz = 2*2 + 2*3*(nx-2) + 2*2
  call create_csr_matrix(n, nz, a)

  a%row_ptr(1) = 1
  inz = 1

  ! left boundary
  a%col_ind(inz) = 1
  a%val(inz) = -dboverdxdx
  inz = inz + 1
  a%col_ind(inz) = 3
  a%val(inz) = dboverdxdx
  inz = inz + 1
  a%row_ptr(2) = inz

  a%col_ind(inz) = 2
  a%val(inz) = -dcoverdxdx
  inz = inz + 1
  a%col_ind(inz) = 4
  a%val(inz) = dcoverdxdx
  inz = inz + 1
  a%row_ptr(3) = inz

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    a%col_ind(inz) = irow-2
    a%val(inz) = dboverdxdx
    inz = inz + 1
    a%col_ind(inz) = irow
    a%val(inz) = -2.d0 * dboverdxdx
    inz = inz + 1
    a%col_ind(inz) = irow+2
    a%val(inz) = dboverdxdx
    inz = inz + 1
    a%row_ptr(irow+1) = inz

    a%col_ind(inz) = irow-1
    a%val(inz) = dcoverdxdx
    inz = inz + 1
    a%col_ind(inz) = irow+1
    a%val(inz) = -2.d0 * dcoverdxdx
    inz = inz + 1
    a%col_ind(inz) = irow+3
    a%val(inz) = dcoverdxdx
    inz = inz + 1
    a%row_ptr(irow+2) = inz

  end do

  ! right boundary
  a%col_ind(inz) = n-3
  a%val(inz) = dboverdxdx
  inz = inz + 1
  a%col_ind(inz) = n-1
  a%val(inz) = -dboverdxdx
  inz = inz + 1
  a%row_ptr(n) = inz

  a%col_ind(inz) = n-2
  a%val(inz) = dcoverdxdx
  inz = inz + 1
  a%col_ind(inz) = n
  a%val(inz) = -dcoverdxdx
  inz = inz + 1
  a%row_ptr(n+1) = inz

  print *, "LS: inz-1 = ", inz-1
  print *, "LS:  nz   = ", nz 

end subroutine build_bz_2eq_diff_matrix

subroutine bz_2eq_init_sol(np, neq, u)
  implicit none
  integer :: np
  integer :: neq
  real(kind=dp), dimension(:) :: u
  ! local
  real(kind=dp), allocatable, dimension(:) :: b, c
  real(kind=dp) :: xcoor, ycoor
  integer :: jy
  real(kind=dp), parameter :: pi = 4.d0 * atan(1.d0)
  integer :: i, j

  allocate(b(np))
  allocate(c(np))

  do  jy = 1, int(np/20)
    xcoor = 0.5d0
    ycoor = 1.d0*(dfloat(jy))/dfloat(np/20)-.05d0

    if (.3d0*xcoor.ge.ycoor .and.  ycoor.ge.0.d0 .and. xcoor.ge.0.d0) then
      b(jy) = .8d0
    else
      b(jy) = q*(f+1)/(f-1)
    endif

    if (xcoor.gt.0.d0) then

      if (ycoor.ge.0.d0) then
        c(jy) = q*(f+1.d0)/(f-1.d0) + (datan(ycoor/xcoor))/(8.d0*pi*f)
      else
        c(jy) = q*(f+1.d0)/(f-1.d0) + (datan(ycoor/xcoor)+2.d0*pi)/(8.d0*pi*f)
      endif

    else

      if (xcoor.lt.0.d0) then
        c(jy) = q*(f+1.d0)/(f-1.d0) + (datan((ycoor)/xcoor)+pi)/(8.d0*pi*f)
      else

        if (ycoor.ge.0.d0) then
          c(jy) = q*(f+1)/(f-1)+1.d0/16.d0/f
        else
          c(jy) = q*(f+1)/(f-1)+3.d0/16.d0/f
        endif

      endif

    endif

  enddo

  do jy=int(np/20)+1,np
    b(jy)=b(int(np/20))
    c(jy)=c(int(np/20))
  enddo

 do i = 1, np
   j = neq*(i-1) + 1
   u(j)   = b(i)
   u(j+1) = c(i)
 end do

 deallocate(b, c)

end subroutine bz_2eq_init_sol

end module mod_bz_2eq_1d
