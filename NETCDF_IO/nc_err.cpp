#include "../netcdf_io.hpp"
#include <string.h>
#include "../constants.hpp"
#include <cassert>
#include <stdexcept>

void NC_ERR(
        const int e,
        const int line_num,
        const char* file_name,
        const bool is_fatal
        ) {

    #if DEBUG >= -2
    // Print error statement
    fprintf(stderr, "Error: [%s] at line %d in %s\n", nc_strerror(e), line_num, file_name);

    // Halt
    if ( is_fatal ) { throw std::runtime_error("Fatal netCDF error. Halting."); }
    #endif

}

