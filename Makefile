# The "system.mk" file in the current directory will contain
# system-specific make variables, most notably the
# C/C++ compiler/linker,
include system.mk
include VERSION

# Debug output level
CFLAGS:=-DDEBUG=1 $(CFLAGS)

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

# Get list of test executables
TEST_CPPS := $(wildcard Tests/*.cpp)
TEST_EXES := $(addprefix Tests/,$(notdir $(TEST_CPPS:.cpp=.x)))

.PHONY: clean hardclean docs cleandocs tests all
clean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o Functions/Differentiation_Tools/*.o Tests/*.o
hardclean:
	rm -f *.o NETCDF_IO/*.o Functions/*.o Functions/Differentiation_Tools/*.o coarse_grain.x Tests/*.o Tests/*.x
	rm -r coarse_grain.x.dSYM
	rm -r Documentation/html
	rm -r Documentation/latex
cleandocs:
	rm -r Documentation/html
	rm -r Documentation/latex

all: coarse_grain.x

tests: ${TEST_EXES}

%.o: %.cpp
	$(MPICXX) ${VERSION} $(LINKS) $(LDFLAGS) -c $(CFLAGS) -o $@ $<

coarse_grain.x: ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} coarse_grain.o
	$(MPICXX) ${VERSION} $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^

Tests/%.x: Tests/%.o ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS}
	$(MPICXX) ${VERSION} $(LINKS) $(CFLAGS) $(LDFLAGS) -o $@ $^

docs:
	DOXY_VERSION=${DOXY_VERSION} doxygen Doxyfile
