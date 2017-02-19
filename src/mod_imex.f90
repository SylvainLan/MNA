module mod_imex

use mpi
use mod_precision
use mod_csr_matrix

contains

subroutine imex_integration(build_matrix, f, n, u, nt, tini, tend)
  implicit none
  include 'dmumps_struc.h'
  ! arguments
  interface
    subroutine build_matrix(n, a)
       use mod_csr_matrix
       integer :: n
       type(csr_matrix) :: a
    end subroutine build_matrix
    function f(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f
    end function f
  end interface
  integer :: n
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  real(kind=dp), dimension(n) :: um1, um2
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
    print *, "  dt :", dt

    call build_matrix(n, a)

    mumps%n  = a%n
    mumps%nz = a%nz
    print *, "  a%n :", a%n
    print *, "  a%nz :", a%nz
 
    allocate(mumps%irn(mumps%nz))
    allocate(mumps%jcn(mumps%nz))
    allocate(mumps%a(mumps%nz))
    allocate(mumps%rhs(mumps%n))

    ! matrix for the first iteration : first order IMEX scheme
    inz = 1
    do irow = 1, a%n
      !!print *, "  irow :", irow
      do iptr = a%row_ptr(irow), a%row_ptr(irow+1)-1
        !!print *, "      icol :", a%col_ind(iptr)
        mumps%irn(inz) = irow
        mumps%jcn(inz) = a%col_ind(iptr)
        if ((a%col_ind(iptr)) == irow) then
          mumps%a(inz) = 1.d0 - dt * a%val(iptr)
        else
          mumps%a(inz) = - dt * a%val(iptr)
        end if
        inz = inz + 1
        !!print *, "      inz :", inz
      end do
    end do

    mumps%rhs = u + dt*f(n,u)

  end if

  !! no outputs
  mumps%icntl(4) = 1

  ! compute solution of first iteration
  mumps%job = 6
  call dmumps(mumps)

  if (mumps%myid == 0) then

    !do irow = 1, a%n
    !  norm = u[
    !end do  

    um1 = mumps%rhs

    inz = 1
    do irow = 1, a%n
      do iptr = a%row_ptr(irow), a%row_ptr(irow+1)-1
        mumps%irn(inz) = irow
        mumps%jcn(inz) = a%col_ind(iptr)
        if ((a%col_ind(iptr)) == irow) then
          mumps%a(inz) = 1.5d0 - dt * a%val(iptr) 
        else 
          mumps%a(inz) = - dt * a%val(iptr) 
        end if
        inz = inz + 1 
      end do
    end do

  end if  

  !!!! no outputs
  !!!!mumps%icntl(4) = 1
  !!
  ! call mumps for analysis step
  mumps%job = 1
  t_mumps = mpi_wtime() 
  call dmumps(mumps)
  t_mumps = mpi_wtime() - t_mumps
  call mpi_reduce(t_mumps, tg_mumps, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr) 
  if (mumps%myid == 0) then
    print *, "  Time (s) for mumps to perform analysis :", tg_mumps
    print *, "  Type of analysis actually done :", mumps%infog(32)
    print *, "  Ordering method actually used  :", mumps%infog(7)
  end if

  ! call mumps to factorisation step
  mumps%job = 2
  t_mumps = mpi_wtime() 
  call dmumps(mumps)
  t_mumps = mpi_wtime() - t_mumps
  call mpi_reduce(t_mumps, tg_mumps, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr) 
  if (mumps%myid == 0) then
    print *, "  Time (s) for mumps to perform factorisation :", tg_mumps
  end if

  ! second-order imex iterations
  t_mumps = mpi_wtime() 
  do it = 2, nt-1 
    if (mumps%myid == 0) then
      um2 = um1
      um1 = mumps%rhs
      mumps%rhs = 2.d0*dt*f(n,um1) - dt*f(n, um2) + 2.d0*um1 - 0.5*um2
    end if
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

end subroutine imex_integration

end module mod_imex
