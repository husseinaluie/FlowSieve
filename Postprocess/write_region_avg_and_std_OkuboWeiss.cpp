#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void write_region_avg_and_std_OkuboWeiss(
        const std::vector< std::vector< double > > & field_averages_OW,
        const std::vector< std::vector< double > > & field_std_devs_OW,
        const std::vector< double > & OkuboWeiss_areas,
        const std::vector<std::string> & vars_to_process,
        const char * filename,
        const int Stime,
        const int Sdepth,
        const int Ntime,
        const int Ndepth,
        const int Nokubo,
        const int num_regions,
        const int num_fields
        ){

    // Dimension order: time - depth - region
    size_t start[4], count[4];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    start[2] = 0;
    count[2] = Nokubo;

    start[3] = 0;
    count[3] = num_regions;

    write_field_to_output(OkuboWeiss_areas, "area_OkuboWeiss", start, count, filename);
    for ( int Ifield = 0; Ifield < num_fields; ++Ifield ) {
        write_field_to_output( field_averages_OW.at(Ifield), vars_to_process.at(Ifield) + "_OkuboWeiss_average", start, count, filename );
    }

}
