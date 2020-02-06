#ifndef POSTPROCESS_HPP
#define POSTPROCESS_HPP 1

#include <vector>
#include <string>
#include <math.h>
#include <mpi.h>

void Apply_Postprocess_Routines(
        const std::vector<double> & coarse_u_r, 
        const std::vector<double> & coarse_u_lon, 
        const std::vector<double> & coarse_u_lat, 
        const std::vector<double> & coarse_vort_r, 
        const std::vector<double> & coarse_vort_lon, 
        const std::vector<double> & coarse_vort_lat,
        const std::vector<double> & energy_transfer, 
        const std::vector<double> & div_J, 
        const std::vector<double> & lambda_m, 
        const std::vector<double> & PEtoKE,
        const std::vector<double> & div,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const double filter_scale,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void write_regions(
        const char * filename,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

namespace RegionTest 
{
    const double D2R = M_PI / 180.;
    const double R2D = 180. / M_PI;

    extern bool Global(double latitude, double longitude);
    extern bool GulfofMexico(double latitude, double longitude);
    extern bool GulfStream(double latitude, double longitude);
    extern bool Equator(double latitude, double longitude);
    extern bool NorthofEquator(double latitude, double longitude);
    extern bool SouthofEquator(double latitude, double longitude);
    extern bool Kuroshio(double latitude, double longitude);
    extern bool AntarcticCircumpolar(double latitude, double longitude);

    // Add a variable which points to all of the ( currently defined ) 
    //    regions. This will make it easier for post-process routines
    //    to simply loop through the whole set.
    // First one points to functions, second one gives (human-readable) names
    //    currently, printing a dimension with string values doesn't seem to 
    //    work in parallel netcdf. Until I find a way around this, I'm just going
    //    to flub it and add a bunch of attributes (hopefully)
    extern const std::vector< bool(*)(double, double) > all_regions;
    //extern char const ** region_names;
    extern const std::vector< std::string > region_names;

}

#endif
