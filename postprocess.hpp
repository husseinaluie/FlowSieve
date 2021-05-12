#ifndef POSTPROCESS_HPP
#define POSTPROCESS_HPP 1

#include <vector>
#include <string>
#include <math.h>
#include <mpi.h>
#include "functions.hpp" // provides classes

void Apply_Postprocess_Routines(
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<std::string> & vars_to_process,
        const std::vector<double> & OkuboWeiss,
        const double filter_scale,
        const std::string filename_base = "postprocess",
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void write_regions(
        const char * filename,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<bool> & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void compute_region_areas(
        std::vector<double> & region_areas,
        const std::vector<double> & areas,
        const std::vector<bool> & mask,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int num_regions,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        );

void compute_region_avg_and_std(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void compute_region_avg_and_std_OkuboWeiss(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        std::vector< double > & OkuboWeiss_areas,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<double> & OkuboWeiss,
        const std::vector<double> OkuboWeiss_bounds,
        const int NOkubo
        );

void compute_time_avg_std(
        std::vector<std::vector<double> > & time_average,
        std::vector<std::vector<double> > & time_std_dev,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<int> & mask_count,
        const std::vector<bool> & always_masked,
        const int full_Ntime,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

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
        );

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
        );

namespace RegionTest 
{
    const double D2R = M_PI / 180.;
    const double R2D = 180. / M_PI;

    extern bool Global(double latitude, double longitude);
    extern bool GulfofMexico(double latitude, double longitude);
    extern bool GulfStream(double latitude, double longitude);
    extern bool Tropics(double latitude, double longitude);
    extern bool NorthofTropics(double latitude, double longitude);
    extern bool SouthofTropics(double latitude, double longitude);
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
