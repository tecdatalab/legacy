CPP=g++
MPICPP=mpiCC.mpich
CFLAGS=-O3 -Wall -Werror
#CFLAGS=-Wall -Werror -ggdb
#CFLAGS=-O3 -Wall -Werror
# mpi.h throws warnings so we can't treat them as errors
MPICFLAGS=-O3 -Wall
#CFLAGS=-ggdb
#FINALFLAGS=-static
FINALFLAGS=
LD=$(CPP)
UNAME=$(shell uname -s)

ANN_LIBS=-L$(HOME)/libraries/ann/lib
ANN_INCLUDE=-I$(HOME)/libraries/ann/include
ANN_LDFLAGS=-lANN

GSL_LIBS=-L$(HOME)/libraries/gsl/lib
GSL_INCLUDE=-I$(HOME)/libraries/gsl/include
GSL_LDFLAGS=-lgsl -lgslcblas

MDGA_OBJS=md_main.o mdstate.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph.o rmsd.o align.o cluster_manager.o point_transformation.o
MDGA_MPI_OBJS=md_main_mpi.o mpimessaging.o mpicoordinator.o mdstate.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph_mpi.o rmsd.o align.o point_transformation.o cluster_manager_mpi.o
MDGA_MPI_NODISTCLUST_OBJS=md_main_mpi.o mpimessaging.o mpicoordinator.o mdstate.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph_mpi.o rmsd.o align.o point_transformation.o cluster_manager.o
MDGA_LIBS=$(ANN_LIBS) $(GSL_LIBS)
MDGA_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
MDGA_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

MDGA_CREATEPDB_OBJS=create_pdb.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph.o
MDGA_CREATEPDB_LIBS=$(ANN_LIBS) $(GSL_LIBS)
MDGA_CREATEPDB_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
MDGA_CREATEPDB_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

MDGA_STATS_OBJS=stats.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph.o rmsd.o align.o point_transformation.o
MDGA_STATS_MPI_OBJS=stats_mpi.o mpimessaging.o mpicoordinator.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph_mpi.o cluster_manager_mpi.o rmsd.o align.o point_transformation.o
MDGA_STATS_LIBS=$(ANN_LIBS) $(GSL_LIBS)
MDGA_STATS_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
MDGA_STATS_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

MDGA_SUBUNIT_STATS_OBJS=multilzerd_subunit_stats.o statssubunits.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph.o rmsd.o align.o point_transformation.o

READER_OBJS=reader.o pdb.o utils.o
READER_LIBS=$(ANN_LIBS)
READER_INCLUDE=$(ANN_INCLUDE)
READER_LDFLAGS=$(ANN_LDFLAGS)

SPLITTER_OBJS=splitter.o pdb.o utils.o
SPLITTER_LIBS=$(ANN_LIBS)
SPLITTER_INCLUDE=$(ANN_INCLUDE)
SPLITTER_LDFLAGS=$(ANN_LDFLAGS)

CLASHES_OBJS=clashing_atoms.o pdb.o utils.o
CLASHES_LIBS=$(ANN_LIBS)
CLASHES_INCLUDE=$(ANN_INCLUDE)
CLASHES_LDFLAGS=$(ANN_LDFLAGS)

SCORE_OBJS=calculate_score.o pdb.o utils.o contact.o DisjointSets.o scoring.o md.o mdgraph.o rmsd.o align.o point_transformation.o
SCORE_LIBS=$(MDGA_STATS_LIBS)
SCORE_INCLUDE=$(MDGA_STATS_INCLUDE)
SCORE_LDFLAGS=$(MDGA_STATS_LDFLAGS)

CLUSTER_OBJS=pairwiseclustering.o pdb.o utils.o contact.o align.o rmsd.o cluster_manager.o DisjointSets.o mdgraph.o scoring.o point_transformation.o 
CLUSTER_LIBS=$(MDGA_STATS_LIBS)
CLUSTER_INCLUDE=$(MDGA_STATS_INCLUDE)
CLUSTER_LDFLAGS=$(MDGA_STATS_LDFLAGS)

MULTICLUSTER_OBJS=multiclustering.o mdstate.o md.o contact.o DisjointSets.o pdb.o scoring.o utils.o mdgraph.o rmsd.o align.o cluster_manager.o point_transformation.o
MULTICLUSTER_MPI_OBJS=multiclustering_mpi.o mdstate.o md.o contact.o DisjointSets.o pdb.o scoring.o utils.o mpimessaging.o mpicoordinator.o mdgraph_mpi.o rmsd.o align.o cluster_manager_mpi.o point_transformation.o
MULTICLUSTER_LIBS=$(ANN_LIBS) $(GSL_LIBS)
MULTICLUSTER_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
MULTICLUSTER_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

REFINE_OBJS=lnsrch.o dfpmin.o montecarlo.o pdb.o utils.o md.o contact.o scoring.o mdgraph.o DisjointSets.o
REFINE_LIBS=$(ANN_LIBS)
REFINE_INCLUDE=$(ANN_INCLUDE)
REFINE_LDFLAGS=$(ANN_LDFLAGS) -lm

PAIRWISE_STATS_OBJS=pairwisestats.o pdb.o utils.o contact.o align.o rmsd.o point_transformation.o
PAIRWISE_STATS_LIBS=$(MDGA_STATS_LIBS)
PAIRWISE_STATS_INCLUDE=$(MDGA_STATS_INCLUDE)
PAIRWISE_STATS_LDFLAGS=$(MDGA_STATS_LDFLAGS)

PDB_ALIGN_OBJS=main_pdb_align.o pdb.o utils.o align.o point_transformation.o
PDB_ALIGN_LIBS=$(ANN_LIBS) $(GSL_LIBS)
PDB_ALIGN_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
PDB_ALIGN_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

PDBGEN_OBJS=main_pdbgen.o pdb.o utils.o
PDBGEN_LIBS=$(ANN_LIBS)
PDBGEN_INCLUDE=$(ANN_INCLUDE)
PDBGEN_LDFLAGS=$(ANN_LDFLAGS)

MRF_SAMPLER_OBJS=main_mrf_sampler.o mrf_sampler.o mrf_features.o pdb.o utils.o align.o transformable_pdb.o point_transformation.o rmsd.o contact.o scoring.o density_lattice.o point_lattice.o mrc.o correlation_results.o
MRF_SAMPLER_LIBS=$(ANN_LIBS) $(GSL_LIBS)
MRF_SAMPLER_INCLUDE=$(ANN_INCLUDE) $(GSL_INCLUDE)
MRF_SAMPLER_LDFLAGS=$(ANN_LDFLAGS) $(GSL_LDFLAGS)

ifeq ($(UNAME),Linux)
BINARIES=multilzerd multilzerd_mpi multilzerd_mpi_nodistclust multilzerd_stats_mpi multilzerd_stats multilzerd_create_pdb multilzerd_refine multilzerd_cluster multilzerd_cluster_mpi reader splitter clashes score lzerd_cluster lzerd_stats multilzerd_subunit_stats pdb_align PDBGEN mrf_sampler
else
BINARIES=multilzerd multilzerd_stats multilzerd_create_pdb multilzerd_refine multilzerd_cluster reader splitter clashes score lzerd_cluster lzerd_stats multilzerd_subunit_stats pdb_align PDBGEN mrf_sampler
endif

all: $(BINARIES)

multilzerd: $(MDGA_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MDGA_INCLUDE) $(MDGA_LIBS) $(MDGA_OBJS) $(MDGA_LDFLAGS) -o $@

ifeq ($(UNAME),Linux)
multilzerd_mpi: $(MDGA_MPI_OBJS)
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(FINALFLAGS) $(MDGA_INCLUDE) $(MDGA_LIBS) $(MDGA_MPI_OBJS) $(MDGA_LDFLAGS) -o $@

multilzerd_mpi_nodistclust: $(MDGA_MPI_NODISTCLUST_OBJS)
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(FINALFLAGS) $(MDGA_INCLUDE) $(MDGA_LIBS) $(MDGA_MPI_NODISTCLUST_OBJS) $(MDGA_LDFLAGS) -o $@
endif

multilzerd_stats: $(MDGA_STATS_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MDGA_STATS_INCLUDE) $(MDGA_STATS_LIBS) $(MDGA_STATS_OBJS) $(MDGA_STATS_LDFLAGS) -o $@

multilzerd_subunit_stats: $(MDGA_SUBUNIT_STATS_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MDGA_STATS_INCLUDE) $(MDGA_STATS_LIBS) $(MDGA_SUBUNIT_STATS_OBJS) $(MDGA_STATS_LDFLAGS) -o $@

ifeq ($(UNAME),Linux)
multilzerd_stats_mpi: $(MDGA_STATS_MPI_OBJS)
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(FINALFLAGS) $(MDGA_STATS_INCLUDE) $(MDGA_STATS_LIBS) $(MDGA_STATS_MPI_OBJS) $(MDGA_STATS_LDFLAGS) -o $@
endif

multilzerd_create_pdb: $(MDGA_CREATEPDB_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MDGA_CREATEPDB_INCLUDE) $(MDGA_CREATEPDB_LIBS) $(MDGA_CREATEPDB_OBJS) $(MDGA_CREATEPDB_LDFLAGS) -o $@

multilzerd_refine: $(REFINE_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(REFINE_INCLUDE) $(REFINE_LIBS) $(REFINE_OBJS) $(REFINE_LDFLAGS) -o $@

reader: $(READER_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(READER_INCLUDE) $(READER_LIBS) $(READER_OBJS) $(READER_LDFLAGS) -o $@

splitter: $(SPLITTER_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(SPLITTER_INCLUDE) $(SPLITTER_LIBS) $(SPLITTER_OBJS) $(SPLITTER_LDFLAGS) -o $@

clashes: $(CLASHES_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(CLASHES_INCLUDE) $(CLASHES_LIBS) $(CLASHES_OBJS) $(CLASHES_LDFLAGS) -o $@

score: $(SCORE_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(SCORE_INCLUDE) $(SCORE_LIBS) $(SCORE_OBJS) $(SCORE_LDFLAGS) -o $@

lzerd_cluster: $(CLUSTER_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(CLUSTER_INCLUDE) $(CLUSTER_LIBS) $(CLUSTER_OBJS) $(CLUSTER_LDFLAGS) -o $@

multilzerd_cluster: $(MULTICLUSTER_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MULTICLUSTER_INCLUDE) $(MULTICLUSTER_LIBS) $(MULTICLUSTER_OBJS) $(MULTICLUSTER_LDFLAGS) -o $@

ifeq ($(UNAME),Linux)
multilzerd_cluster_mpi: $(MULTICLUSTER_MPI_OBJS)
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(FINALFLAGS) $(MULTICLUSTER_INCLUDE) $(MULTICLUSTER_LIBS) $(MULTICLUSTER_MPI_OBJS) $(MULTICLUSTER_LDFLAGS) -o $@
endif

lzerd_stats: $(PAIRWISE_STATS_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(PAIRWISE_STATS_INCLUDE) $(PAIRWISE_STATS_LIBS) $(PAIRWISE_STATS_OBJS) $(PAIRWISE_STATS_LDFLAGS) -o $@

pdb_align: $(PDB_ALIGN_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(PDB_ALIGN_INCLUDE) $(PDB_ALIGN_LIBS) $(PDB_ALIGN_OBJS) $(PDB_ALIGN_LDFLAGS) -o $@ 

PDBGEN: $(PDBGEN_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(PDBGEN_INCLUDE) $(PDBGEN_LIBS) $(PDBGEN_OBJS) $(PDBGEN_LDFLAGS) -o $@ 

pdb.o multilzerd_subunit_stats.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) -c $(subst .o,.cc,$@)

mrf_sampler: $(MRF_SAMPLER_OBJS)
	$(LD) $(CFLAGS) $(FINALFLAGS) $(MRF_SAMPLER_INCLUDE) $(MRF_SAMPLER_LIBS) $(MRF_SAMPLER_OBJS) $(MRF_SAMPLER_LDFLAGS) -o $@ 

create_pdb.o contact.o scoring.o mdgraph.o md_main.o stats.o rmsd.o calculate_score.o montecarlo.o pairwiseclustering.o pairwisestats.o mdstate.o multiclustering.o cluster_manager.o statssubunits.o main_pdb_align.o mrf_sampler.o:
	$(CPP) $(CFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst .o,.cc,$@)

ifeq ($(UNAME),Linux)
mdgraph_mpi.o: mdgraph.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

cluster_manager_mpi.o: cluster_manager.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

stats_mpi.o: stats.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

md_main_mpi.o: md_main.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@

multiclustering_mpi.o: multiclustering.cc
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) $(ANN_INCLUDE) $(GSL_INCLUDE) -c $(subst _mpi.o,.cc,$@) -o $@
	
mpimessaging.o mpicoordinator.o:
	$(MPICPP) -DWITH_MPI $(MPICFLAGS) -c $(subst .o,.cc,$@)
endif

align.o:
	$(CPP) $(CFLAGS) $(GSL_INCLUDE) -c $(subst .o,.cc,$@)

density_lattice.o: ../em/density_lattice.h ../em/point_lattice.h ../em/mrc.h ../em/correlation_results.h
	$(CPP) $(CFLAGS) -c -o $@ ../em/density_lattice.cc

point_lattice.o: ../em/point_lattice.h
	$(CPP) $(CFLAGS) -c -o $@ ../em/point_lattice.cc

mrc.o: ../em/mrc.h ../em/density_lattice.h ../em/point_lattice.h
	$(CPP) $(CFLAGS) -c -o $@ ../em/mrc.cc

correlation_results.o: ../em/correlation_results.h
	$(CPP) $(CFLAGS) -c -o $@ ../em/correlation_results.cc

install:
	cp $(BINARIES) ~/bin

clean:
	rm -f $(MDGA_OBJS) $(MDGA_MPI_OBJS) $(MDGA_STATS_OBJS) $(MDGA_STATS_MPI_OBJS) $(MDGA_CREATEPDB_OBJS) $(READER_OBJS) $(SPLITTER_OBJS) $(CLASHES_OBJS) $(SCORE_OBJS) $(REFINE_OBJS) $(MULTICLUSTER_OBJS) $(MULTICLUSTER_MPI_OBJS) $(CLUSTER_OBJS) $(PAIRWISE_STATS_OBJS) $(MDGA_SUBUNIT_STATS_OBJS) $(PDB_ALIGN_OBJS) $(PDBGEN_OBJS) $(MRF_SAMPLER_OBJS) deps $(BINARIES)

.SUFFIXES: .cc .cpp .o

.cc.o:
	$(CPP) $(CFLAGS) -c $< -o $@

.cpp.o:
	$(CPP) $(CFLAGS) -c $< -o $@


deps:
	$(CPP) -MM *.cc $(GSL_INCLUDE) $(ANN_INCLUDE) > deps

include deps
