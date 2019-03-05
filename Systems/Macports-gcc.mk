# Specify compilers
CXX     ?= gcc-mp-8
MPICXX  ?= mpicxx-openmpi-gcc8

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl

# Default compiler flags
CFLAGS:=-Wall 

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3

# Extra optimization flags
EXTRA_OPT_FLAGS:=

# Add Macports-installed libraries to the path
LIB_DIRS:=-L /opt/local/lib
INC_DIRS:=-I /opt/local/include
