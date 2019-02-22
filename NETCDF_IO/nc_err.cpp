/*
 *
 * Short-hand function to handle netcdf errors
 *
 * Prints: the netcdf error
 *         the line number at which the error occurred
 *         the file in which the error occurred
 *
 */

#include "../netcdf_io.hpp"
#include <string.h>

#ifndef DEBUG
    #define DEBUG 0
#endif

void NC_ERR(const int e, const int line_num, const char* file_name) {
    #if DEBUG <= -2
    fprintf(stderr, "Error: [%s] at line %d in %s\n", nc_strerror(e), line_num, file_name);
    #endif
}

