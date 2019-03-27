
#include "../netcdf_io.hpp"
#include <string.h>
#include "../constants.hpp"

void NC_ERR(
        const int e,          /**< [in] Error code provided by netcdf function */
        const int line_num,   /**< [in] The line number at which the error occurred */
        const char* file_name /**< [in] The name of the file in which the error occurred. */
        ) {
    #if DEBUG <= -2
    fprintf(stderr, "Error: [%s] at line %d in %s\n", nc_strerror(e), line_num, file_name);
    #endif
}

