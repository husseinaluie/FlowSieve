#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <mpi.h>

#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../constants.hpp"
#include "../functions.hpp"

int main(int argc, char **argv)
{
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string &zonal_vel_name    = input.getCmdOption("--zonal_vel",   "uo");
    const std::string &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo");
    const std::string &salinity_name     = input.getCmdOption("--salinity",    "so");
    const std::string &temperature_name  = input.getCmdOption("--temperature", "thetao");
    const std::string &sla_name          = input.getCmdOption("--sla",         "zos");
    const std::string &input_fname       = input.getCmdOption("--input_file",  "input.nc");
    const std::string &output_fname      = input.getCmdOption("--output_file", "interp.nc");

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth,
                        u_lat, thetao, sal, ssh, 
                        mask_thetao, mask_sal, mask_ssh, mask;
    std::vector<int> myCounts, myStarts;

    read_var_from_file(longitude, "longitude", input_fname, NULL);
    read_var_from_file(latitude,  "latitude",  input_fname, NULL);
    read_var_from_file(time,      "time",      input_fname, NULL);
    read_var_from_file(depth,     "depth",     input_fname, NULL);

    convert_coordinates(longitude, latitude);

    read_var_from_file(thetao, temperature_name,    input_fname, &mask_thetao, &myCounts, &myStarts );
    read_var_from_file(sal ,   salinity_name,       input_fname, &mask_sal );
    read_var_from_file(ssh,    sla_name,            input_fname, &mask_ssh );

    // Combine masks (only keep points where none are masked)
    mask.resize( mask_thetao.size() );
    for (size_t II = 0; II < mask_thetao.size(); ++II) {
        mask.at(II) = mask_thetao.at(II) * mask_sal.at(II) * mask_ssh.at(II);
    }

    // Free up memory from the old masks
    mask_thetao.clear();  mask_thetao.shrink_to_fit();
    mask_sal.clear();     mask_sal.shrink_to_fit();
    mask_ssh.clear();     mask_ssh.shrink_to_fit();

    // Compute the density and pressure at each water point
    //    by converting first, it means that we only need
    //    to interpolate two fields instead of three, which
    //    saves us a bit.
    std::vector<double> pressure(        sal.size() ),
                        pressure_interp( sal.size() ),
                        density(         sal.size() ),
                        density_interp(  sal.size() );

    double temp;
    for (size_t index = 0; index < pressure.size(); ++index) {
        if (mask.at(index) == 1) {
            // Compute pressure (geostrophic)
            pressure.at(index) = constants::rho0 * constants::g * ssh.at(index);

            // Convert potential temp to actual temp
            temp = depotential_temperature(pressure.at(index), thetao.at(index));

            // Compute density from UNESCO equation of state
            density.at(index) = equation_of_state(
                    temp, 
                    sal.at(index), 
                    pressure.at(index) / 1e5); // takes pressure in bar = 100 kPa = 1e5 Pa

            if ( (std::isnan(pressure.at(index))) or (std::isnan(density.at(index))) ) {
                if (std::isnan(pressure.at(index))) { fprintf( stdout, "  NaN pressure  :" ); }
                if (std::isnan(density.at( index))) { fprintf( stdout, "  NaN density  :" ); }
                fprintf(stdout, "  ssh, theta (temp), sal = %g, %g (%g), %g \n",
                        ssh.at(index), thetao.at(index), temp, sal.at(index) );
            }
        }
    }

    // Now output to interp.nc
    #if DEBUG >= 0
    fprintf(stdout, "Creating output\n");
    #endif
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("pressure");
    vars_to_write.push_back("density");

    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(myCounts.at(0)), size_t(myCounts.at(1)), 
        size_t(myCounts.at(2)), size_t(myCounts.at(3))};

    // Outputting pre-interpolated fields
    initialize_output_file(time, depth, longitude, latitude, mask, 
            vars_to_write, "pre_interp.nc", 0.);

    write_field_to_output(pressure, "pressure", starts, counts, "pre_interp.nc", &mask);
    write_field_to_output(density,  "density",  starts, counts, "pre_interp.nc", &mask);

    // Apply interpolation
    #if DEBUG >= 0
    fprintf(stdout, "Interpolating density\n");
    #endif
    interpolate_over_land_from_coast(density_interp, density, 
            time, depth, latitude, longitude, mask, myCounts);

    #if DEBUG >= 0
    fprintf(stdout, "Interpolating pressure\n");
    #endif
    interpolate_over_land_from_coast(pressure_interp, pressure, 
            time, depth, latitude, longitude, mask, myCounts);

    // Outputting interpolated fields
    initialize_output_file(time, depth, longitude, latitude, mask, vars_to_write, output_fname.c_str(), 0.);

    write_field_to_output(pressure_interp, "pressure", starts, counts, output_fname);
    write_field_to_output(density_interp,  "density",  starts, counts, output_fname);

    return 0;
}

