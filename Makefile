# The "system.mk" file in the current directory will contain
# system-specific make variables, most notably the
# C/C++ compiler/linker,
include system.mk
include VERSION

# Debug output level
CFLAGS:=-DDEBUG=2 $(CFLAGS)

# Do you want vorticity computed?
CFLAGS:=-DCOMP_VORT=true $(CFLAGS)

# Do you want energy transfers computed?
CFLAGS:=-DCOMP_TRANSFERS=true $(CFLAGS)

# Do you want baroclinic energy transfers computed?
CFLAGS:=-DCOMP_BC_TRANSFERS=true $(CFLAGS)

# Turn on/off debug flags or additional optimization flags
OPT:=true
DEBUG:=false
EXTRA_OPT:=false

##
## Shouldn't need to modify anything beyond this point
##

# Flag to pass version info to code
VERSION:= -DMAJOR_VERSION=${MAJOR_VERSION} -DMINOR_VERSION=${MINOR_VERSION} -DPATCH_VERSION=${PATCH_VERSION}
DOXY_VERSION:="${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}"

ifeq ($(DEBUG),true)
	CFLAGS:=$(CFLAGS) $(DEBUG_FLAGS)
	LINKS:=$(LINKS) $(DEBUG_LDFLAGS)
endif

ifeq ($(OPT),true)
    CFLAGS:=$(CFLAGS) $(OPT_FLAGS)
endif

ifeq ($(EXTRA_OPT),true)
	CFLAGS:=$(CFLAGS) $(EXTRA_OPT_FLAGS)
endif

CFLAGS:=$(CFLAGS) $(LIB_DIRS)
LDFLAGS:=$(LDFLAGS) $(INC_DIRS)

# Get list of netcdf IO cpp files
NETCDF_IO_CPPS := $(wildcard NETCDF_IO/*.cpp)
NETCDF_IO_OBJS := $(addprefix NETCDF_IO/,$(notdir $(NETCDF_IO_CPPS:.cpp=.o)))

# Get list of function cpp files
FUNCTIONS_CPPS := $(wildcard Functions/*.cpp)
FUNCTIONS_OBJS := $(addprefix Functions/,$(notdir $(FUNCTIONS_CPPS:.cpp=.o)))

# Get list of differentation cpp files
DIFF_TOOL_CPPS := $(wildcard Functions/Differentiation_Tools/*.cpp)
DIFF_TOOL_OBJS := $(addprefix Functions/Differentiation_Tools/,$(notdir $(DIFF_TOOL_CPPS:.cpp=.o)))

# Get list of ALGLIB object files
ALGLIB_CPPS := $(wildcard ALGLIB/*.cpp)
ALGLIB_OBJS := $(addprefix ALGLIB/,$(notdir $(ALGLIB_CPPS:.cpp=.o)))

# Get the list of preprocessing  cpp files
PREPROCESS_CPPS := $(wildcard Preprocess/*.cpp)
PREPROCESS_OBJS := $(addprefix Preprocess/,$(notdir $(PREPROCESS_CPPS:.cpp=.o)))

# Get list of test executables
TEST_CPPS := $(wildcard Tests/*.cpp)
TEST_EXES := $(addprefix Tests/,$(notdir $(TEST_CPPS:.cpp=.x)))

.PHONY: clean hardclean docs cleandocs tests all ALGLIB
clean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o Functions/Differentiation_Tools/*.o Tests/*.o Preprocess/*.o
hardclean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o Functions/Differentiation_Tools/*.o coarse_grain.x Tests/*.o Tests/*.x
	rm -r coarse_grain.x.dSYM
	rm ALGLIB/*.o
	rm -r Documentation/html
	rm -r Documentation/latex
cleandocs:
	rm -r Documentation/html
	rm -r Documentation/latex

all: coarse_grain.x

tests: ${TEST_EXES}

ALGLIB: ${ALGLIB_OBJS}

# ALGLIB requires a different set of compile flags
ALGLIB/%.o: ALGLIB/%.cpp
	$(CXX) -I ./ALGLIB -c -O3 -DAE_CPU=AE_INTEL -o $@ $<

# The core coarse_graining functions use the same compile flags
Functions/Differentiation_Tools/%.o: Functions/Differentiation_Tools/%.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

Functions/%.o: Functions/%.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

NETCDF_IO/%.o: NETCDF_IO/%.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

Preprocess/%.o: Preprocess/%.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $<

# Building test scripts
Tests/%.o: Tests/%.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

Tests/%.x: Tests/%.o ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS}
	$(MPICXX) ${VERSION} $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^

# Building coarse_grain executable
coarse_grain.o: coarse_grain.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

coarse_grain.x: ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} coarse_grain.o
	$(MPICXX) ${VERSION} $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^

# Interpolator needs to link in ALGLIB
interpolator.o: interpolator.cpp
	$(MPICXX) $(LINKS) $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $<

interpolator.x: ${NETCDF_IO_OBJS} ${ALGLIB_OBJS} ${PREPROCESS_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} interpolator.o
	$(MPICXX) $(LINKS) $(CFLAGS) -I ./ALGLIB $(LDFLAGS) -o $@ $^

# Documentation build with Doxygen
docs:
	DOXY_VERSION=${DOXY_VERSION} doxygen Doxyfile
