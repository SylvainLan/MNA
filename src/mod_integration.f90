module mod_integration

use mpi
use mod_precision
use mod_radau
use mod_strang
use mod_imex
use mod_bz_2eq_1d

contains

subroutine integrate(method, tol, neq, nx, u, nt, tini, tend)
  implicit none
  ! arguments
  character(len=20) :: method
  real(kind=dp) :: tol
  integer :: neq
  integer :: nx
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  integer :: ntot
  integer :: proc_nb
  integer :: ierr

  ntot = neq * nx

  call mpi_comm_rank(mpi_comm_world, proc_nb, ierr) 
 
  select case (method) 
    case ("radau5")
      if (proc_nb == 0) then
        print *
        print *, "Radau5 integration"
        call radau5_integration(f_bz_2eq, tol, ntot, u, nt, tini, tend, .true.)
      endif
    case ("strang")
      if (proc_nb == 0) then
        print *
        print *, "Strang splitting integration"
        print *, "  RADAU5 for reaction"
        print *, "  ROCK4  for diffusion"
        call strang_integration(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend)
      end if
    case ("imex")
      if (proc_nb == 0) then
        print *
        print *, "Second-order IMEX integration"
      endif
      call imex_integration(build_bz_2eq_diff_matrix, f_bz_2eq_reac, ntot, u, nt, tini, tend)
    case default
      print *
      print *, method, "Unknown integration method"
      call mpi_finalize(ierr)
      call exit(0)
  end select

end subroutine integrate

end module mod_integration
