module mod_backward_euler

use mpi
use mod_precision
use mod_csr_matrix

contains

subroutine backward_euler_integration(build_matrix, n, u, nt, tini, tend)
  implicit none
  include 'dmumps_struc.h'
  ! arguments
  interface
    subroutine build_matrix(n, a)
       use mod_csr_matrix
       integer :: n
       type(csr_matrix) :: a
    end subroutine build_matrix
  end interface
  integer :: n
  real(kind=dp), dimension(n) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  type (dmumps_struc) :: mumps
  type(csr_matrix) :: a
  real(dp) :: dt, dx
  integer :: irow, iptr, inz, it
  real(dp) :: t_mumps, tg_mumps
  integer :: nb_procs, ierr

  mumps%comm = mpi_comm_world
  mumps%par  = 1
  mumps%sym  = 0
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)

  !! initialize an instance of mumps
  mumps%job  = -1
  call dmumps(mumps)

  if (mumps%myid == 0) then

    print *, "  Parallel version of mumps"
    print *, "  Number of process used :", nb_procs

    dt = (tend-tini)/(nt-1)

    call build_matrix(n, a)

    mumps%n  = a%n
    mumps%nz = a%nz
 
    allocate(mumps%irn(mumps%nz))
    allocate(mumps%jcn(mumps%nz))
    allocate(mumps%a(mumps%nz))
    allocate(mumps%rhs(mumps%n))

    inz = 1
    do irow = 1, a%n
      do iptr = a%row_ptr(irow), a%row_ptr(irow+1)-1
        mumps%irn(inz) = irow
        mumps%jcn(inz) = a%col_ind(iptr)
        if ((a%col_ind(iptr)) == irow) then
          mumps%a(inz) = 1.d0 - dt * a%val(iptr) 
        else 
          mumps%a(inz) = - dt * a%val(iptr) 
        end if
        inz = inz + 1 
      end do
    end do

    mumps%rhs = u 

  end if

  !! no outputs
  mumps%icntl(4) = 1

  !! parallel computation of analysis step
  !!mumps%icntl(28) = 2
  !!mumps%icntl(29) = 2
  !! use Scotch
  !!mumps%icntl(7) = 7

  ! call mumps for analysis step
  mumps%job = 1
  t_mumps = mpi_wtime() 
  call dmumps(mumps)
  t_mumps = mpi_wtime() - t_mumps
  call mpi_reduce(t_mumps, tg_mumps, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr) 
  if (mumps%myid == 0) then
    print *, "  Time (s) for mumps to perform analysis :", tg_mumps
    print *, "  Type of analysis actually done :", mumps%infog(32)
    print *, "  ordering method actually used  :", mumps%infog(7)
  end if

  ! call mumps for factorisation step
  mumps%job = 2
  t_mumps = mpi_wtime() 
  call dmumps(mumps)
  t_mumps = mpi_wtime() - t_mumps
  call mpi_reduce(t_mumps, tg_mumps, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr) 
  if (mumps%myid == 0) then
    print *, "  Time (s) for mumps to perform factorisation :", tg_mumps
  end if

  ! euler iterations
  t_mumps = mpi_wtime() 
  do it = 1, nt-1 
    mumps%job = 3
    call dmumps(mumps)
  end do 
  t_mumps = mpi_wtime() - t_mumps
  call mpi_reduce(t_mumps, tg_mumps, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr) 

  if (mumps%myid == 0) then
    print *, "  Total   time (s) for mumps to compute solutions :", tg_mumps
    print *, "  Average time (s) for mumps to compute solutions :", tg_mumps/(nt-1)
    u = mumps%rhs
    ! deallocate   
    deallocate(mumps%irn)
    deallocate(mumps%jcn)
    deallocate(mumps%a)
    deallocate(mumps%rhs)
  end if

  !! destroy the instance of mumps
  mumps%job = -2
  call dmumps(mumps)

end subroutine backward_euler_integration

end module mod_backward_euler
