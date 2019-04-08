# Specify compilers
CXX     ?= g++ #icpc
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -fopenmp

# Default compiler flags
CFLAGS:=-Wall 

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3 #-fp-model fast=2

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=-ip -ipo

# Modules are automatically on lib dir
LIB_DIRS:=
INC_DIRS:=
