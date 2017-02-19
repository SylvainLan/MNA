!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program bz_2eq_1d_main

  use mpi
  use mod_precision
  use mod_cartesian_grid
  use mod_utils
  use mod_bz_2eq_1d
  use mod_integration

  implicit none

  integer :: neq = 2
  real(kind=dp) :: tini, tend
  integer :: nt
  real(kind=dp) :: xmin, xmax
  integer :: nxib
  type(cartesian_grid_type) :: grid
  real(kind=dp), allocatable, dimension(:) :: unum
  integer :: ntot = 2
  character(len=20) :: method
  real(kind=dp) :: tol
  real(kind=dp) :: norm_err
  real(kind=dp) :: loc_elapsed, elapsed
  integer :: proc_nb, nb_procs
  integer :: ierr

  ! init mpi environnement for mumps
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nb_procs, ierr)
  call mpi_comm_rank(mpi_comm_world, proc_nb, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! initialisation

  call read_data(tini, tend, nt, xmin, xmax, nxib, method, tol)

  if (proc_nb == 0) then

    print *, "Resolution of 1d heat equation" 
 
    print "(a, f12.3)", "   tini = ", tini
    print "(a, f12.3)", "   tend = ", tend
    print *, "  nt   =", nt
    print "(a, f12.3)", "   xmin =", xmin
    print "(a, f12.3)", "   xmax =", xmax
    print *, "  nxib =", nxib
    print *, "  integration method : ", method
    print "(a, es12.3)", "   tolerance (for Radau5 and Rock4) =", tol
 
    call init_cartesian_grid(xmin, xmax, nxib, grid)
  
    ! init heat 
    call init_bz_2eq(grid)
 
    ! allocation
    ntot = neq*grid%nx
    allocate(unum(ntot)) 
  
    ! compute and save initial solution 
    call bz_2eq_init_sol(grid%nx, neq, unum)
    call save_sol("sol_ini.dat", neq, grid, unum)

  endif
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! integration
  loc_elapsed = mpi_wtime()
  call integrate(method, tol, neq, grid%nx, unum, nt, tini, tend)
  loc_elapsed = mpi_wtime() - loc_elapsed
  call mpi_reduce(loc_elapsed, elapsed, 1, mpi_double, mpi_max, 0, mpi_comm_world, ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! print sol and compute error
  if (proc_nb == 0) then

    print *
    print *, "Time (s) to integrate :", elapsed
 
    ! write numeric solution
    call save_sol("sol_num.dat", neq, grid, unum)
 
    ! compute error
    call compute_error(neq, grid, unum)
   
    deallocate(unum) 

  end if
 
  ! terminate mpi environment
  call mpi_finalize(ierr)

end program bz_2eq_1d_main
