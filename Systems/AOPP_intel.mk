# This make file is for the AOPP Oxford Physics department cluster in 2023.
# Before compilation run the following in the terminal to load the relevant modules and then build the executable:
# module load intel-compilers/2022
# module load openmpi/4.1.4-intel
# module load hdf5/1.12.2-intel-parallel
# module load netcdf/netcdf-c-4.9.0-parallel

# Specify compilers
CXX     ?= g++
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -fopenmp -lm -lsz -lbz2 -lxml2

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3

# Modules are automatically on lib dir
NETCDF_LIBS=-L/network/software/ubuntu_bionic/netcdf/netcdf-c-4.9.0-parallel/lib -L/network/software/ubuntu_bionic/hdf5/1.12.2-intel-parallel/lib
NETCDF_INCS=-I/network/software/ubuntu_bionic/netcdf/netcdf-c-4.9.0-parallel/include -I/network/software/ubuntu_bionic/hdf5/1.12.2-intel-parallel/include


LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
