# The "system.mk" file in the current directory will contain
# system-specific make variables, most notably the
# C/C++ compiler/linker,
include system.mk
include VERSION

# Debug output level
CFLAGS:=-DDEBUG=1 $(CFLAGS)

# Turn on/off debug flags or additional optimization flags
OPT:=true
DEBUG:=false
EXTRA_OPT:=true

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

# Get list of interface cpp files
INTERFACE_CPPS := $(wildcard Functions/Interface_Tools/*.cpp)
INTERFACE_OBJS := $(addprefix Functions/Interface_Tools/,$(notdir $(INTERFACE_CPPS:.cpp=.o)))

# Get list of fftw-based files
FFT_BASED_CPPS := $(wildcard Functions/FFTW_versions/*.cpp)
FFT_BASED_OBJS := $(addprefix Functions/FFTW_versions/,$(notdir $(FFT_BASED_CPPS:.cpp=.o)))

# Get list of ALGLIB object files
ALGLIB_CPPS := $(wildcard ALGLIB/*.cpp)
ALGLIB_OBJS := $(addprefix ALGLIB/,$(notdir $(ALGLIB_CPPS:.cpp=.o)))

# Get the list of preprocessing  cpp files
PREPROCESS_CPPS := $(wildcard Preprocess/*.cpp)
PREPROCESS_OBJS := $(addprefix Preprocess/,$(notdir $(PREPROCESS_CPPS:.cpp=.o)))

# Get the list of post-processing  cpp files
POSTPROCESS_CPPS := $(wildcard Postprocess/*.cpp)
POSTPROCESS_OBJS := $(addprefix Postprocess/,$(notdir $(POSTPROCESS_CPPS:.cpp=.o)))

# Get list of test executables
TEST_CPPS := $(wildcard Tests/*.cpp)
TEST_EXES := $(addprefix Tests/,$(notdir $(TEST_CPPS:.cpp=.x)))


.PHONY: clean hardclean docs cleandocs tests all ALGLIB
clean:
	rm -f *.o 
	rm -f NETCDF_IO/*.o 
	rm -f Functions/*.o 
	rm -f Functions/Differentiation_Tools/*.o 
	rm -f Functions/Interface_Tools/*.o 
	rm -f Functions/FFTW_versions/*.o 
	rm -f Tests/*.o 
	rm -f Preprocess/*.o
	rm -f Postprocess/*.o
hardclean:
	rm -f *.o 
	rm -f *.x
	rm -f NETCDF_IO/*.o 
	rm -f Functions/*.o 
	rm -f Functions/Differentiation_Tools/*.o 
	rm -f Tests/*.[o,x] 
	rm -f Functions/Interface_Tools/*.o 
	rm -f Functions/FFTW_versions/*.o 
	rm -f coarse_grain.x 
	rm -r coarse_grain.x.dSYM
	rm ALGLIB/*.o
	rm -r docs/html
	rm -r docs/latex
	rm -f Preprocess/*.[o,x]
	rm -f Postprocess/*.[o,x]
cleandocs:
	rm -r docs/html
	rm -r docs/latex

all: coarse_grain.x integrator.x ${TEST_EXES}

tests: ${TEST_EXES}

ALGLIB: ${ALGLIB_OBJS}


#
# Commands for building object files (does not include those for executables)
#

# ALGLIB requires a different set of compile flags
ALGLIB/%.o: ALGLIB/%.cpp
	$(CXX) -I ./ALGLIB -c -O3 -DAE_CPU=AE_INTEL -o $@ $<

# Group together object files with similar compilation requirements
CORE_OBJS := ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} ${POSTPROCESS_OBJS}

$(CORE_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

$(INTERFACE_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

FFT_BASED_OBJS:= ${FFT_BASED_OBJS} coarse_grain_fftw.o
$(FFT_BASED_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< -lfftw3 -lm $(LINKS) 

$(PREPROCESS_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $< $(LINKS) 


# 
# Commands for building executables (and related object files)
#

# Group together executables with similar compilations
CORE_TARGET_EXES := coarse_grain.x integrator.x coarse_grain_sw.x coarse_grain_subset.x
CORE_TARGET_OBJS := coarse_grain.o integrator.o coarse_grain_sw.o coarse_grain_subset.o

$(CORE_TARGET_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

$(CORE_TARGET_EXES): %.x : ${CORE_OBJS} ${INTERFACE_OBJS} %.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 

# Building test scripts
Tests/%.o: Tests/%.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

Tests/%.x: Tests/%.o ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS}
	$(MPICXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 

# Building fftw-based coarse_grain executable
coarse_grain_fftw.x: ${CORE_OBJS} ${INTERFACE_OBJS} ${FFT_BASED_OBJS} coarse_grain_fftw.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ -lfftw3 -lm $(LINKS) 

# Interpolator needs to link in ALGLIB
interpolator.o: interpolator.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $< $(LINKS) 

interpolator.x: ${NETCDF_IO_OBJS} ${ALGLIB_OBJS} ${PREPROCESS_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} interpolator.o
	$(MPICXX) $(CFLAGS) -I ./ALGLIB $(LDFLAGS) -o $@ $^ $(LINKS) 


#
# Documentation
#


# Documentation build with Doxygen
docs:
	DOXY_VERSION=${DOXY_VERSION} doxygen Doxyfile
