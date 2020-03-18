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

    // Default input file name
    const int str_len = 50;
    char input_fname[     str_len], zonal_vel_name[  str_len], merid_vel_name[str_len], 
         salinity_name[   str_len], temperature_name[str_len], sla_name[      str_len],
         flag_version[    str_len], flag_input[      str_len], flag_salinity[ str_len],
         flag_temperature[str_len], flag_sla[        str_len], flag_merid_vel[str_len],
         flag_zonal_vel[  str_len];

    snprintf(input_fname,       str_len, "input.nc");
    snprintf(zonal_vel_name,    str_len, "uo");
    snprintf(merid_vel_name,    str_len, "vo");
    snprintf(salinity_name,     str_len, "so");
    snprintf(temperature_name,  str_len, "thetao");
    snprintf(sla_name,          str_len, "zos");

    snprintf(flag_version,      str_len, "--version");
    snprintf(flag_input,        str_len, "--input_file");
    snprintf(flag_zonal_vel,    str_len, "--zonal_vel");
    snprintf(flag_merid_vel,    str_len, "--merid_vel");
    snprintf(flag_salinity,     str_len, "--salinity");
    snprintf(flag_temperature,  str_len, "--temperature");
    snprintf(flag_sla,          str_len, "--sla");

    // Parse command-line flags
    int ii = 1;
    while(ii < argc) {  
        if (wRank == 0) {
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
        }

        if (strcmp(argv[ii], flag_version) == 0) {
            // check if the flag is 'version'
            //if (wRank == 0) {
            //    print_compile_info(filter_scales);
            //}
            return 0;
        } else if (strcmp(argv[ii], flag_input) == 0) {
            // check if the flag is 'input_file' and, if it is
            //   the next input is then used as the filename
            //   of the input
            snprintf(input_fname, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_zonal_vel) == 0) {
            // check if we're given the name of the zonal velocity var
            snprintf(zonal_vel_name, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_merid_vel) == 0) {
            // check if we're given the name of the meridional velocity var
            snprintf(merid_vel_name, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_salinity) == 0) {
            // check if we're given the name of the salinity var
            snprintf(salinity_name, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_temperature) == 0) {
            // check if we're given the name of the temperature var
            snprintf(temperature_name, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_sla) == 0) {
            // check if we're given the name of the sla var
            snprintf(sla_name, str_len, argv[ii+1]);
            ++ii;
        } else {
            // Otherwise, the flag is unrecognized
            if (wRank == 0) {
                fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
            }
            return -1;
        }
        ++ii;
    }

    std::vector<double> longitude, latitude, time, depth,
                        u_lat, thetao, sal, ssh, mask;
    std::vector<int> myCounts, myStarts;

    read_var_from_file(longitude, "longitude", input_fname, NULL);
    read_var_from_file(latitude,  "latitude",  input_fname, NULL);
    read_var_from_file(time,      "time",      input_fname, NULL);
    read_var_from_file(depth,     "depth",     input_fname, NULL);

    read_var_from_file(u_lat,  merid_vel_name,      input_fname, &mask, 
            &myCounts, &myStarts);
    read_var_from_file(thetao, temperature_name,    input_fname, NULL );
    read_var_from_file(sal ,   salinity_name,       input_fname, NULL );
    read_var_from_file(ssh,    sla_name,            input_fname, NULL );

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
        }
    }

    // Apply interpolation
    #if DEBUG >= 0
    fprintf(stdout, "Interpolating density\n");
    #endif
    interpolate_over_land_from_coast(density_interp, density, 
            time, depth, latitude, longitude, mask);

    #if DEBUG >= 0
    fprintf(stdout, "Interpolating pressure\n");
    #endif
    interpolate_over_land_from_coast(pressure_interp, pressure, 
            time, depth, latitude, longitude, mask);

    // Now output to interp.nc
    #if DEBUG >= 0
    fprintf(stdout, "Creating output\n");
    #endif
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("pressure_orig");
    vars_to_write.push_back("density_orig");
    vars_to_write.push_back("pressure_interp");
    vars_to_write.push_back("density_interp");
    initialize_output_file(time, depth, longitude, latitude, mask, 
            vars_to_write, "interp.nc", 0.);

    #if DEBUG >= 0
    fprintf(stdout, "Writing appropriate variables to interp.nc\n");
    #endif
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(myCounts.at(0)), size_t(myCounts.at(1)), 
        size_t(myCounts.at(2)), size_t(myCounts.at(3))};
    write_field_to_output(pressure, "pressure_orig", starts, counts, "interp.nc", &mask);
    write_field_to_output(density,  "density_orig",  starts, counts, "interp.nc", &mask);

    write_field_to_output(pressure_interp, "pressure_interp", starts, counts, "interp.nc");
    write_field_to_output(density_interp,  "density_interp",  starts, counts, "interp.nc");

    return 0;
}

