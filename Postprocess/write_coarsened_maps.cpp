#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void write_coarsened_maps(
        const std::vector< std::vector< double > > & coarsened_maps,
        const std::vector<std::string> & vars_to_process,
        const char * filename,
        const int Stime,
        const int Sdepth,
        const int Ntime,
        const int Ndepth,
        const int coarse_Nlat,
        const int coarse_Nlon,
        const int num_fields
        ){

    // Dimension order: time - depth - lat - lon
    size_t start[4], count[4];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    start[2] = 0;
    count[2] = coarse_Nlat;

    start[3] = 0;
    count[3] = coarse_Nlon;

    for ( int Ifield = 0; Ifield < num_fields; ++Ifield ) {
        write_field_to_output( coarsened_maps.at(Ifield), vars_to_process.at(Ifield) + "_coarsened_map", start, count, filename );
    }

}
