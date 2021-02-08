#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void write_region_avg_and_std(
        const std::vector< std::vector< double > > & field_averages,
        const std::vector< std::vector< double > > & field_std_devs,
        const std::vector<std::string> & vars_to_process,
        const char * filename,
        const int Stime,
        const int Sdepth,
        const int Ntime,
        const int Ndepth,
        const int num_regions,
        const int num_fields
        ){

    // Dimension order: time - depth - region
    size_t start[3], count[3];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    start[2] = 0;
    count[2] = num_regions;

    for ( int Ifield = 0; Ifield < num_fields; ++Ifield ) {
        write_field_to_output( field_averages.at(Ifield), vars_to_process.at(Ifield) + "_avg", start, count, filename );
    }

}
