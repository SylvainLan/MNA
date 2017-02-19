
FC = mpif90
LD = mpif90

INCLUDE_DIR = /applications/libraries/mumps/5.0.2/include

SCOTCH_LIBS = -L/applications/libraries/mumps/5.0.2/ptscotch/6.0.4/lib -lptesmumps -lptscotch -lptscotcherr -lscotch
METIS_LIBS  = -L/applications/libraries/mumps/5.0.2/parmetis/4.0.3/lib -lparmetis -lmetis
MKL_LIBS    = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIBS = -L/applications/libraries/mumps/5.0.2/lib -ldmumps -lmumps_common -lpord $(SCOTCH_LIBS) $(METIS_LIBS) $(MKL_LIBS)

FCFLAGS = -qopenmp -module ./obj -I$(INCLUDE_DIR) 
LDFLAGS = -qopenmp
