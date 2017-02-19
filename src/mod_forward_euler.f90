module mod_forward_euler

use mod_precision

contains

subroutine forward_euler_integration(f, n, u, nt, tini, tend)
  implicit none
  ! arguments
  interface
    function f(n, u)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: u, f
    end function f
  end interface
  integer :: n
  real(kind=dp), dimension(n) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  real(dp) :: dt
  integer :: it
 
  ! init 
  dt = (tend-tini)/(nt-1)
 
  ! euler
  do it = 1, nt-1
    u = u + (f(n,u)*dt)
  end do

end subroutine forward_euler_integration

end module mod_forward_euler
