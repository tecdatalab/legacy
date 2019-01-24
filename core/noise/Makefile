UNAME=$(shell uname)
CPP=g++
MPICPP=mpiCC.mpich
WARNINGSOKFLAGS=-O3 -Wall
CFLAGS=-O3 -Wall -Werror
MPICFLAGS=-O3 -Wall
GLUTLINUXFLAGS=-lglut -lGL -lGLU
GLUTMACFLAGS=-framework GLUT -framework OpenGL -framework Cocoa
ifeq ($(UNAME),Linux)
  GLUTFLAGS=$(GLUTLINUXFLAGS)
else
  GLUTFLAGS=$(GLUTMACFLAGS)
endif

ANN_LIBS=-L$(HOME)/libraries/ann/lib
ANN_INCLUDE=-I$(HOME)/libraries/ann/include
ANN_LDFLAGS=-lANN

GSL_LIBS=-L$(HOME)/libraries/gsl/lib
GSL_INCLUDE=-I$(HOME)/libraries/gsl/include
GSL_LDFLAGS=-lgsl -lgslcblas

FINALFLAGS=-static
LD=$(CPP)

EM_CORRELATION_OBJS=mrc.o correlation_results.o density_lattice.o point_lattice.o point_transformation.o main_em_correlation.o
EM_MPI_CORRELATION_OBJS=mrc.o correlation_results.o density_lattice_mpi.o point_lattice.o point_transformation.o main_em_correlation_mpi.o
EM_CORRELATION_CLUSTERING_OBJS=mrc.o correlation_results.o density_lattice.o point_lattice.o point_transformation.o main_em_correlation_clustering.o
EM_VIEWER_OBJS=mrc.o correlation_results.o density_lattice.o point_lattice.o point_transformation.o marching_cubes_surface.o main_em_viewer.o opengl_utils.o pdb_painter.o pdb.o rmsd.o contact.o align.o utils.o transformable_pdb.o
EM_VOLUME_OBJS=mrc.o correlation_results.o density_lattice.o point_lattice.o point_transformation.o main_em_volume.o
EM_RESET_ORIGIN_OBJS=mrc.o correlation_results.o density_lattice.o point_lattice.o pdb.o point_transformation.o utils.o main_em_reset_origin.o

EM_REMOVE_NOISE_OBJS=mrc.o density_lattice.o point_lattice.o main_test_noise_removal.o point_transformation.o correlation_results.o

ifeq ($(UNAME),Linux)
BINARIES=em_correlation em_correlation_clustering em_correlation_mpi em_viewer em_volume em_reset_origin
else
BINARIES=em_correlation em_correlation_clustering em_viewer em_volume em_reset_origin
endif

all: $(BINARIES)

em_correlation: $(EM_CORRELATION_OBJS)
	$(LD) $(CFLAGS) $(EM_CORRELATION_OBJS) -o $@

em_correlation_mpi: $(EM_MPI_CORRELATION_OBJS)
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(EM_MPI_CORRELATION_OBJS) -o $@

em_correlation_clustering: $(EM_CORRELATION_CLUSTERING_OBJS)
	$(LD) $(CFLAGS) $(EM_CORRELATION_CLUSTERING_OBJS) -o $@

ifeq ($(UNAME),Linux)
em_viewer:
else
em_viewer: $(EM_VIEWER_OBJS)
	$(LD) $(WARNINGSOKFLAGS) $(GLUTFLAGS) $(ANN_LIBS) $(ANN_LDFLAGS) $(GSL_LIBS) $(GSL_LDFLAGS) $(EM_VIEWER_OBJS) -o $@
endif

em_volume: $(EM_VOLUME_OBJS)
	$(LD) $(CFLAGS) $(EM_VOLUME_OBJS) -o $@

em_reset_origin: $(EM_RESET_ORIGIN_OBJS)
	$(LD) $(CFLAGS) $(ANN_LIBS) $(EM_RESET_ORIGIN_OBJS) $(ANN_LDFLAGS) -o $@

em_remove_noise: $(EM_REMOVE_NOISE_OBJS)
	$(LD) $(CFLAGS) $(EM_REMOVE_NOISE_OBJS) -o $@

main_em_correlation_mpi.o: main_em_correlation.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) -c $(subst _mpi.o,.cc,$@) -o $@

density_lattice_mpi.o: density_lattice.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) -c $(subst _mpi.o,.cc,$@) -o $@

install:
	cp $(BINARIES) ~/bin

clean:
	rm -f $(EM_CORRELATION_OBJS) $(EM_MPI_CORRELATION_OBJS) $(EM_CORRELATION_CLUSTERING_OBJS) $(EM_VIEWER_OBJS) $(EM_VOLUME_OBJS) $(EM_RESET_ORIGIN_OBJS) deps $(BINARIES)

.SUFFIXES: .cc .o

ifneq ($(UNAME),Linux)
opengl_utils.o:
	$(CPP) $(WARNINGSOKFLAGS) -c $(subst .o,.cc,$@) -o $@
endif

pdb_painter.o: ../proteinstructure/pdb_painter.h ../proteinstructure/pdb_painter.cc
	$(CPP) $(CFLAGS) -c -o $@ ../proteinstructure/pdb_painter.cc

point_transformation.o: ../proteinstructure/point_transformation.h ../proteinstructure/point_transformation.cc
	$(CPP) $(CFLAGS) -c -o $@ ../proteinstructure/point_transformation.cc

transformable_pdb.o: ../proteinstructure/transformable_pdb.h ../proteinstructure/transformable_pdb.cc
	$(CPP) $(CFLAGS) -c -o $@ ../proteinstructure/transformable_pdb.cc

pdb.o: ../proteinstructure/pdb.h ../proteinstructure/atom.h ../proteinstructure/utils.h
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) -c -o $@ ../proteinstructure/pdb.cc

rmsd.o: ../proteinstructure/rmsd.h ../proteinstructure/atom.h ../proteinstructure/rmsd.cc
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c -o $@ ../proteinstructure/rmsd.cc

align.o: ../proteinstructure/align.h ../proteinstructure/atom.h ../proteinstructure/align.cc
	$(CPP) $(CFLAGS) $(GSL_INCLUDE) -c -o $@ ../proteinstructure/align.cc

contact.o: ../proteinstructure/contact.h ../proteinstructure/atom.h ../proteinstructure/contact.cc
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c -o $@ ../proteinstructure/contact.cc

utils.o: ../proteinstructure/utils.h
	$(CPP) $(CFLAGS) -c -o $@ ../proteinstructure/utils.cc

main_em_viewer.o:
	$(CPP) $(WARNINGSOKFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst .o,.cc,$@)

.cc.o:
	$(CPP) $(CFLAGS) -c $< -o $@

deps:
	$(CPP) -MM *.cc $(GSL_INCLUDE) $(ANN_INCLUDE) > deps

include deps
