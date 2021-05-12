# The "system.mk" file in the current directory will contain
# system-specific make variables, most notably the
# C/C++ compiler/linker,
include system.mk
include VERSION

# Debug output level
CFLAGS:=-DDEBUG=0 $(CFLAGS)

# Turn on/off debug flags or additional optimization flags
OPT:=true
DEBUG:=false
EXTRA_OPT:=true
USE_GPROF:=false

##
## Shouldn't need to modify anything beyond this point
##

# Flag to pass version info to code
VERSION:= -DMAJOR_VERSION=${MAJOR_VERSION} -DMINOR_VERSION=${MINOR_VERSION} -DPATCH_VERSION=${PATCH_VERSION}
DOXY_VERSION:="${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}"
GIT_VERSION := "$(shell git describe --dirty --always --all --long)"
VERSION:=${VERSION} -DGIT_VERSION=\"$(GIT_VERSION)\"

ifeq ($(DEBUG),true)
    CFLAGS:=$(CFLAGS) $(DEBUG_FLAGS)
    LINKS:=$(LINKS) $(DEBUG_LDFLAGS)
endif

ifeq ($(OPT),true)
    CFLAGS:=$(CFLAGS) $(OPT_FLAGS)
endif

ifeq ($(EXTRA_OPT),true)
    CFLAGS:=$(CFLAGS) $(EXTRA_OPT_FLAGS)
    ALGLIB_OPT_FLAGS:=$(ALGLIB_OPT_FLAGS) $(EXTRA_OPT_FLAGS)
endif

ifeq ($(USE_GPROF),true)
    LINKS:=$(LINKS) -pg
endif

CFLAGS:=$(CFLAGS) $(LIB_DIRS)
LDFLAGS:=$(LDFLAGS) $(INC_DIRS)

##
## Group the object files into handy names.
##    Also specify the compile command for them.
##


# Get list of netcdf IO cpp files
NETCDF_IO_CPPS := $(wildcard  NETCDF_IO/*.cpp)
NETCDF_IO_OBJS := $(addprefix NETCDF_IO/,$(notdir $(NETCDF_IO_CPPS:.cpp=.o)))

# Get list of function cpp files
FUNCTIONS_CPPS := $(wildcard  Functions/*.cpp)
FUNCTIONS_OBJS := $(addprefix Functions/,$(notdir $(FUNCTIONS_CPPS:.cpp=.o)))

# Get list of differentation cpp files
DIFF_TOOL_CPPS := $(wildcard  Functions/Differentiation_Tools/*.cpp)
DIFF_TOOL_OBJS := $(addprefix Functions/Differentiation_Tools/,$(notdir $(DIFF_TOOL_CPPS:.cpp=.o)))

# Get the list of post-processing  cpp files
POSTPROCESS_CPPS := $(wildcard  Postprocess/*.cpp)
POSTPROCESS_OBJS := $(addprefix Postprocess/,$(notdir $(POSTPROCESS_CPPS:.cpp=.o)))

# Get the list of particle cpp files
PARTICLE_CPPS := $(wildcard  Functions/Particles/*.cpp)
PARTICLE_OBJS := $(addprefix Functions/Particles/,$(notdir $(PARTICLE_CPPS:.cpp=.o)))

# Group together object files with similar compilation requirements
CORE_OBJS := ${NETCDF_IO_OBJS} ${FUNCTIONS_OBJS} \
	${DIFF_TOOL_OBJS} ${POSTPROCESS_OBJS} ${PARTICLE_OBJS}

$(CORE_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 


# Get list of SW cpp files
SW_TOOL_CPPS := $(wildcard  Functions/SW_Tools/*.cpp)
SW_TOOL_OBJS := $(addprefix Functions/SW_Tools/,$(notdir $(SW_TOOL_CPPS:.cpp=.o)))

$(SW_TOOL_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 


# Get list of interface cpp files
INTERFACE_CPPS := $(wildcard  Functions/Interface_Tools/*.cpp)
INTERFACE_OBJS := $(addprefix Functions/Interface_Tools/,$(notdir $(INTERFACE_CPPS:.cpp=.o)))

$(INTERFACE_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 


# Get list of fftw-based files
FFT_BASED_CPPS := $(wildcard  Functions/FFTW_versions/*.cpp)
FFT_BASED_OBJS := $(addprefix Functions/FFTW_versions/,$(notdir $(FFT_BASED_CPPS:.cpp=.o)))

FFT_BASED_OBJS:= ${FFT_BASED_OBJS} Case_Files/coarse_grain_fftw.o
$(FFT_BASED_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< -lfftw3 -lm $(LINKS) 


# Get list of toroidal projection files
TOROIDAL_CPPS := $(wildcard  Functions/Toroidal_Projection/*.cpp)
TOROIDAL_OBJS := $(addprefix Functions/Toroidal_Projection/,$(notdir $(TOROIDAL_CPPS:.cpp=.o)))

$(TOROIDAL_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 


# Get list of toroidal projection files
HELMHOLTZ_CPPS := $(wildcard  Functions/Helmholtz/*.cpp)
HELMHOLTZ_OBJS := $(addprefix Functions/Helmholtz/,$(notdir $(HELMHOLTZ_CPPS:.cpp=.o)))

$(HELMHOLTZ_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS)


# Get list of ALGLIB object files
ALGLIB_CPPS := $(wildcard  ALGLIB/*.cpp)
ALGLIB_OBJS := $(addprefix ALGLIB/,$(notdir $(ALGLIB_CPPS:.cpp=.o)))

ALGLIB/%.o: ALGLIB/%.cpp
	$(CXX) -I ./ALGLIB -c $(ALGLIB_OPT_FLAGS) -o $@ $<


# Get the list of preprocessing  cpp files
PREPROCESS_CPPS := $(wildcard  Preprocess/*.cpp)
PREPROCESS_OBJS := $(addprefix Preprocess/,$(notdir $(PREPROCESS_CPPS:.cpp=.o)))

$(PREPROCESS_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $< $(LINKS) 


# Get list of test executables
TEST_CPPS := $(wildcard Tests/*.cpp)
TEST_EXES := $(addprefix Tests/,$(notdir $(TEST_CPPS:.cpp=.x)))


.PHONY: clean hardclean docs cleandocs tests all ALGLIB
clean:
	rm -f *.o 
	rm -f NETCDF_IO/*.o 
	rm -f Functions/*.o 
	rm -f Functions/Differentiation_Tools/*.o 
	rm -f Functions/Helmholtz/*.o 
	rm -f Functions/SW_Tools/*.o 
	rm -f Functions/Interface_Tools/*.o 
	rm -f Functions/FFTW_versions/*.o 
	rm -f Functions/Particles/*.o 
	rm -f Tests/*.o 
	rm -f Preprocess/*.o
	rm -f Postprocess/*.o
	rm -f Case_Files/*.o
hardclean:
	rm -f *.[o,x] 
	rm -f NETCDF_IO/*.o 
	rm -f Functions/*.o 
	rm -f Functions/Differentiation_Tools/*.o 
	rm -f Functions/Helmholtz/*.o 
	rm -f Functions/SW_Tools/*.o 
	rm -f Functions/Interface_Tools/*.o 
	rm -f Functions/FFTW_versions/*.o 
	rm -f Functions/Particles/*.o 
	rm -f Tests/*.[o,x] 
	rm -f Preprocess/*.o
	rm -f Postprocess/*.o
	rm -f Case_Files/*.[o,x]
	rm ALGLIB/*.o
cleandocs:
	rm -r docs/html
	rm -r docs/latex

all: coarse_grain.x integrator.x ${TEST_EXES}

tests: ${TEST_EXES}

ALGLIB: ${ALGLIB_OBJS}


# 
# Commands for building executables (and related object files)
#

# Group together executables with similar compilations
CORE_TARGET_EXES := Case_Files/coarse_grain.x \
					Case_Files/particles.x \
					Case_Files/compare_particles.x \
					Case_Files/project_onto_particles.x
CORE_TARGET_OBJS := Case_Files/coarse_grain.o \
					Case_Files/particles.o \
					Case_Files/compare_particles.o \
					Case_Files/project_onto_particles.o

$(CORE_TARGET_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

$(CORE_TARGET_EXES): %.x : ${CORE_OBJS} ${INTERFACE_OBJS} %.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 


# Helmholtz
HELM_TARGET_EXES := Case_Files/coarse_grain_helmholtz.x 
HELM_TARGET_OBJS := Case_Files/coarse_grain_helmholtz.o 

$(HELM_TARGET_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

$(HELM_TARGET_EXES): %.x : ${CORE_OBJS} ${INTERFACE_OBJS} ${PREPROCESS_OBJS} ${ALGLIB_OBJS} ${HELMHOLTZ_OBJS} %.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 

# Building shallow water
SW_TARGET_EXES := Case_Files/coarse_grain_sw.x
SW_TARGET_OBJS := Case_Files/coarse_grain_sw.o 

$(SW_TARGET_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

$(SW_TARGET_EXES): %.x : ${SW_TOOL_OBJS} ${CORE_OBJS} ${INTERFACE_OBJS} %.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 

# Building toroidal projection
TOROID_TARGET_EXES := 	Case_Files/toroidal_projection.x \
						Case_Files/potential_projection.x \
						Case_Files/interpolator.x \
						Case_Files/geostrophic_vel.x \
						Case_Files/coarsen_grid.x \
						Case_Files/refine_Helmholtz_seed.x
TOROID_TARGET_OBJS := 	Case_Files/toroidal_projection.o \
						Case_Files/potential_projection.o \
						Case_Files/interpolator.o \
						Case_Files/geostrophic_vel.o \
						Case_Files/coarsen_grid.o \
						Case_Files/refine_Helmholtz_seed.o

$(TOROID_TARGET_OBJS): %.o : %.cpp constants.hpp
	$(MPICXX) ${VERSION} $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $< $(LINKS) 

$(TOROID_TARGET_EXES): %.x : ${CORE_OBJS} ${INTERFACE_OBJS} ${PREPROCESS_OBJS} ${ALGLIB_OBJS} %.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -I ./ALGLIB -o $@ $^ $(LINKS) 

# Building test scripts
Tests/%.o: Tests/%.cpp constants.hpp
	$(MPICXX) $(LDFLAGS) -c $(CFLAGS) -o $@ $< $(LINKS) 

Tests/%.x: Tests/%.o ${DIFF_TOOL_OBJS} ${CORE_OBJS} ${INTERFACE_OBJS}
	$(MPICXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LINKS) 

# Building fftw-based coarse_grain executable
Case_Files/coarse_grain_fftw.x: ${CORE_OBJS} ${INTERFACE_OBJS} ${FFT_BASED_OBJS} Case_Files/coarse_grain_fftw.o
	$(MPICXX) ${VERSION} $(CFLAGS) $(LDFLAGS) -o $@ $^ -lfftw3 -lm $(LINKS) 

# Interpolator needs to link in ALGLIB
#Case_Files/interpolator.o: Case_Files/interpolator.cpp constants.hpp
#	$(MPICXX) $(LDFLAGS) -I ./ALGLIB -c $(CFLAGS) -o $@ $< $(LINKS) 
#
#Case_Files/interpolator.x: ${NETCDF_IO_OBJS} ${ALGLIB_OBJS} ${PREPROCESS_OBJS} ${FUNCTIONS_OBJS} ${DIFF_TOOL_OBJS} Case_Files/interpolator.o
#	$(MPICXX) $(CFLAGS) -I ./ALGLIB $(LDFLAGS) -o $@ $^ $(LINKS) 


#
# Documentation
#


# Documentation build with Doxygen
docs:
	DOXY_VERSION=${DOXY_VERSION} doxygen Doxyfile
