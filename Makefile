CXX     ?= gcc-mp-8
MPICXX  ?= mpicxx-openmpi-gcc8
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
CFLAGS:=-O3 -Wall #-qopenmp -fp-model fast=2

DEBUG_FLAGS:=-g -DDEBUG
DEBUG_LDFLAGS:=-g

EXTRA_OPT_FLAGS:=#-ip -ipo

DEBUG:=true
EXTRA_OPT:=false

LIB_DIRS:=-L /opt/local/lib
INC_DIRS:=-I /opt/local/include

ifeq ($(DEBUG),true)
	CFLAGS:=$(CFLAGS) $(DEBUG_FLAGS)
	LINKS:=$(LINKS) $(DEBUG_LDFLAGS)
endif

ifeq ($(EXTRA_OPT),true)
	CFLAGS:=$(CFLAGS) $(EXTRA_OPT_FLAGS)
endif

CFLAGS:=$(CFLAGS) $(LIB_DIRS)
LDFLAGS:=$(LDFLAGS) $(INC_DIRS)

NETCDF_IO_CPPS := $(wildcard NETCDF_IO/*.cpp)
NETCDF_IO_OBJS := $(addprefix NETCDF_IO/,$(notdir $(NETCDF_IO_CPPS:.cpp=.o)))

FUNCTIONS_CPPS := $(wildcard Functions/*.cpp)
FUNCTIONS_OBJS := $(addprefix FUNCTIONS/,$(notdir $(FUNCTIONS_CPPS:.cpp=.o)))

.PHONY: clean
clean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o

all: filter.x

%.o: %.cpp
	$(MPICXX) $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

filter.x: ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} filter.o
	$(MPICXX) $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^
