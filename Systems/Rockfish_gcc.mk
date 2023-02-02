# The following modules were last used
#
#  Modules: 
#  gcc/9.3.0
#  openmpi/3.1.6
#  netcdf-c/4.7.4
#
#  Note: load python/anaconda modules *before* 
#  		 loading those listed above, otherwise 
#  		 netcdf libraries may have conflicts
#
#
#

# Specify compilers
CXX     ?= g++
MPICXX  ?= mpic++

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -fopenmp -lm

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g3 -ggdb
DEBUG_LDFLAGS:=-g3 -ggdb

# Basic optimization flags
OPT_FLAGS:=-O3

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3

# Modules are automatically on lib dir
NETCDF_LIBS="-L/data/apps/linux-centos8-cascadelake/gcc-9.3.0/netcdf-c-4.7.4-cblqp3w3sfun4bnlysfs2paclfihrugz/lib"
NETCDF_INCS=`nc-config --cflags`

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
