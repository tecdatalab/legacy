#####
# Makefile for Situs 2.7
#####

SRC_DIR	= .
BIN_DIR = bin
SHELL = /bin/sh

#####
# name (directory path) of compiler to be used
#####

CC         = gcc
CXX        = g++

#####
# any extra compiler and linker options (see the compiler man page)
#####

CFLAGS     = -O3
XFLAGS     = -O3
LFLAGS     = -lm 

#####
# enable POSIX thread support (comment out if not supported)
#####

SMP        = -D_SMP_ -DMAX_NUM_THREADS=16
SMPLIB     = -lpthread

#####
# enable OPENMP thread support (comment out if not supported)
#####

OMP        = -fopenmp
OMPLIB     = -lgomp 

#####
# Situs files
#####

BINARIES = em2zernike


#####
# misc commands
#####

CD = cd
RM = rm -f
MV = mv
OBJ_FILES = $(wildcard *.o)
MKDIR = mkdir -p

#####
# main makefile rules
#####

default: $(BINARIES)

all: $(BINARIES)

em2zernike:	lib_err.o lib_vwk.o lib_vec.o lib_std.o lib_vio.o 
	$(CXX) main.cpp lib_vio.o lib_vwk.o lib_vec.o lib_std.o lib_err.o $(LFLAGS) -o em2zernike $(XFLAGS) $(SMP) $(OMP)

clean:
	$(RM) $(OBJ_FILES)
	$(MKDIR) $(BIN_DIR)
	$(MV) $(BINARIES) $(BIN_DIR)/.

veryclean:
	$(RM) $(OBJ_FILES)
	$(RM) -r $(BIN_DIR)

.c.o:
	$(CC) $(CFLAGS) $(SMP) $(OMP) -c $(SRC_DIR)/$<
.cpp.o:
	$(CXX) $(XFLAGS) $(SMP) $(OMP) -c $(SRC_DIR)/$<
