module mod_brusselator

use mod_mesh

!! Brusselator with diffusion
!!
!! du/dt = a - (b+1).u + a.u.u.v + alpha * d2u/dxdx
!! dv/dt = b.u - u.u.v           + alpha * d2v/dxdx
!!
!! dy/dt = R.y + alpha * D.y with y = | u
!!                                    | v  

integer, parameter, private :: dp = kind(1.d0)
real(kind=dp) :: a
real(kind=dp) :: b
real(kind=dp) :: alpha
type(mesh1d) :: m
real(kind=dp), allocatable, dimension(:) :: y, yini

contains

subroutine init_brusselator()
  implicit none
    
  a = 1.d0
  b = 3.d0

  allocate(y(2), yini(2))

end subroutine init_brusselator

subroutine init_brusselator_diff()
  implicit none
  real(kind=dp) :: xmin = 0.d0
  real(kind=dp) :: xmax = 1.d0
  integer :: nx = 101
    
  a = 1.d0
  b = 3.d0
  alpha = 1.d0 / 50.d0

  call init_mesh1d(xmin, xmax, nx, m)

  allocate(y(2*m%nx), yini(2*m%nx))

  call init_sol_diff()

end subroutine init_brusselator_diff

subroutine init_sol()
  implicit none
  ! local

  yini(1) = 1.5d0
  yini(2) = 3.d0

end subroutine init_sol

subroutine init_sol_diff()
  implicit none
  real(kind=dp), parameter :: pi = 4.d0 * atan(1.d0)
  ! local
  real(kind=dp) :: ui, vi, xi
  integer :: i, j

  open(unit=10, file='sol_ini.dat')

  do i = 1, m%nx

    xi = m%xmin + (i-1)*m%dx
    ui = 1.d0 + sin(2.d0*pi*xi)
    vi = 3.d0

    j = 2*i - 1
    yini(j)   = ui
    yini(j+1) = vi

   
    write(10,*) xi, ui, vi

  end do

  close(10)

end subroutine init_sol_diff

function f_brusselator(neq, y)
  implicit none
  integer :: neq
  real(kind=dp), dimension(neq) :: y
  real(kind=dp), dimension(neq) :: f_brusselator

  f_brusselator(1) = a + y(1)*y(1)*y(2) - (b + 1.d0)*y(1)
  f_brusselator(2) = b*y(1) - y(1)*y(1)*y(2)
end function


function f_brusselator_diff(neq, y)
  implicit none
  integer :: neq
  real(kind=dp), dimension(neq) :: y
  real(kind=dp), dimension(neq) :: f_brusselator_diff
  ! local variables
  real(kind=dp) :: uim1, ui, uip1
  real(kind=dp) :: vim1, vi, vip1
  real(kind=dp) :: dxdx
  integer :: i, j

  
  dxdx = m%dx * m%dx

  ! left boundary
  ui   = y(1)
  vi   = y(2)
  uip1 = y(3)
  vip1 = y(4)

  f_brusselator_diff(1) = a + ui*ui*vi - (b+1.d0)*ui + alpha*((uip1 - ui)/dxdx)
  f_brusselator_diff(2) = b*ui - ui*ui*vi + alpha*((vip1 - vi)/dxdx)

  do i = 2, m%nx-1
 
    j = 2*i - 1

    uim1 = y(j-2)
    vim1 = y(j-1)

    ui   = y(j)
    vi   = y(j+1)
  
    uip1 = y(j+2)
    vip1 = y(j+3)

    f_brusselator_diff(j)   = a + ui*ui*vi - (b+1.d0)*ui + alpha*((uip1 - 2.d0*ui + uim1)/dxdx)
    f_brusselator_diff(j+1) = b*ui - ui*ui*vi + alpha*((vip1 - 2.d0*vi + vim1)/dxdx)

  end do
 
  ! left boundary
  uim1 = y(neq-3)
  vim1 = y(neq-2)
  ui   = y(neq-1)
  vi   = y(neq)

  f_brusselator_diff(neq-1) = a + ui*ui*vi - (b+1.d0)*ui + alpha*((uim1 - ui)/dxdx)
  f_brusselator_diff(neq)   = b*ui - ui*ui*vi + alpha*((vim1 - vi)/dxdx)
 
end function f_brusselator_diff

end module mod_brusselator
