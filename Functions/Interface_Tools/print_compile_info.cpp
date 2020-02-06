#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "../../constants.hpp"
#include "../../functions.hpp"

void print_compile_info(
        const std::vector<double> &scales   /**< [in] Scales used for filtering */
        ) {

    fprintf(stdout, "\n");
    fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
    fprintf(stdout, "  Code version: %d.%d.%d \n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    #if defined(__ICC) || defined(__INTEL_COMPILER)
    fprintf(stdout, "  compiler: intel (%s)\n", __VERSION__);
    #elif defined(__GNUC__) || defined(__GNUG__)
    fprintf(stdout, "  compiler: gnu (%s)\n", __VERSION__);
    #else
    fprintf(stdout, "  compiler: unknown\n");
    #endif
    /*
    C++ pre-C++98: __cplusplus is 1.
    C++98: __cplusplus is 199711L.
    C++11: __cplusplus is 201103L.
    C++14: __cplusplus is 201402L.
    C++17: __cplusplus is 201703L.
    */
    switch(__cplusplus) {
        case 1 :
            fprintf(stdout, "  c++ version: pre C++98\n");
            break;
        case 199711L :
            fprintf(stdout, "  c++ version: C++98\n");
            break;
        case 201103L :
            fprintf(stdout, "  c++ version: C++11\n");
            break;
        case 201402L :
            fprintf(stdout, "  c++ version: C++14\n");
            break;
        case 201703L :
            fprintf(stdout, "  c++ version: C++17\n");
            break;
        default :
            fprintf(stdout, "  c++ version: unknown (%zu)\n", __cplusplus);
            break;
    }
    fprintf(stdout, "\n");

    fprintf(stdout, "Flags\n");
    fprintf(stdout, "  DEBUG             = %d\n", DEBUG);
    fprintf(stdout, "  KERNEL_OPT        = %d\n", constants::KERNEL_OPT);
    fprintf(stdout, "  CARTESIAN         = %s\n", constants::CARTESIAN         ? "true" : "false");
    fprintf(stdout, "\n");
    fprintf(stdout, "  PERIODIC_X        = %s\n", constants::PERIODIC_X        ? "true" : "false");
    fprintf(stdout, "  PERIODIC_Y        = %s\n", constants::PERIODIC_Y        ? "true" : "false");
    fprintf(stdout, "  UNIFORM_LAT_GRID  = %s\n", constants::UNIFORM_LAT_GRID  ? "true" : "false");
    fprintf(stdout, "\n");
    fprintf(stdout, "  COMP_VORT         = %s\n", constants::COMP_VORT         ? "true" : "false");
    fprintf(stdout, "  COMP_TRANSFERS    = %s\n", constants::COMP_TRANSFERS    ? "true" : "false");
    fprintf(stdout, "  COMP_BC_TRANSFERS = %s\n", constants::COMP_BC_TRANSFERS ? "true" : "false");
    fprintf(stdout, "\n");
    fprintf(stdout, "  MINIMAL_OUTPUT    = %s\n", constants::MINIMAL_OUTPUT    ? "true" : "false");
    fprintf(stdout, "  CAST_TO_INT       = %s\n", constants::CAST_TO_INT       ? "true" : "false");
    fprintf(stdout, "\n");
    fprintf(stdout, "  DO_TIMING         = %s\n", constants::DO_TIMING         ? "true" : "false");
    fprintf(stdout, "\n");

    fprintf(stdout, "Constants\n");
    fprintf(stdout, "  R_earth = %g\n", constants::R_earth);
    fprintf(stdout, "  rho0    = %g\n", constants::rho0);
    fprintf(stdout, "  g       = %g\n", constants::g);
    fprintf(stdout, "  DiffOrd = %d\n", constants::DiffOrd);
    fprintf(stdout, "  KernPad = %g\n", constants::KernPad);
    fprintf(stdout, "\n");

    fprintf(stdout, "%zu Filter Scales (km)\n", scales.size());
    fprintf(stdout, "  ");
    for (size_t II = 0; II < scales.size(); ++II) {
        fprintf(stdout, "%g", scales.at(II)/1e3);
        if (II < scales.size() - 1) { fprintf(stdout, ",  "); }
    }
    fprintf(stdout, "\n\n");
}
