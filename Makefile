CXX     ?= gcc-mp-8
MPICXX  ?= mpicxx-openmpi-gcc8
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
CFLAGS:=-O3 -Wall #-qopenmp -fp-model fast=2

# Debug output level
CFLAGS:=-DDEBUG=1 $(CFLAGS)

# Do you want vorticity computed?
CFLAGS:=-DCOMP_VORT=true $(CFLAGS)

# Do you want energy transfers computed?
CFLAGS:=-DCOMP_TRANSFERS=true $(CFLAGS)

# Turn on/off debug flags or additional optimization flags
DEBUG:=true
EXTRA_OPT:=false

##
## Shouldn't need to modify anything beyond this point
##

DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

EXTRA_OPT_FLAGS:=#-ip -ipo

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
FUNCTIONS_OBJS := $(addprefix Functions/,$(notdir $(FUNCTIONS_CPPS:.cpp=.o)))

.PHONY: clean hardclean docs
clean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o
hardclean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o coarse_grain.x
	rm -r coarse_grain.x.dSYM
	rm -r Documentation/html
	rm -r Documentation/latex

all: coarse_grain.x

%.o: %.cpp
	$(MPICXX) $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

coarse_grain.x: ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} coarse_grain.o
	$(MPICXX) $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^

docs:
	doxygen Doxyfile
