#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "netcdf_io.hpp"
#include "preprocess.hpp"
#include "constants.hpp"
#include "functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

int main(int argc, char **argv)
{

    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> time;
    std::vector<double> depth;
    std::vector<double> scales;
    scales.resize(1);
    scales.at(0) = 0;

    std::vector<double> u_lat;

    std::vector<double> thetao;
    std::vector<double> sal;
    std::vector<double> ssh;

    std::vector<double> mask;

    read_var_from_file(longitude, "longitude", "source.nc", NULL);
    read_var_from_file(latitude,  "latitude",  "source.nc", NULL);
    read_var_from_file(time,      "time",      "source.nc", NULL);
    read_var_from_file(depth,     "depth",     "source.nc", NULL);

    int Ntime = time.size();
    int Ndepth = depth.size();
    int Nlat = latitude.size();
    int Nlon = longitude.size();

    read_var_from_file(u_lat, "vo", "source.nc", &mask);

    read_var_from_file(thetao, "thetao", "source.nc", NULL);
    read_var_from_file(sal ,   "so",     "source.nc", NULL);
    read_var_from_file(ssh,    "zos",    "source.nc", NULL);

    // Compute the density and pressure at each land point
    //    by converting first, it means that we only need
    //    to interpolate two fields instead of three, which
    //    saves us a bit.
    std::vector<double> pressure(sal.size());
    std::vector<double> density(sal.size());

    std::vector<double> pressure_interp(sal.size());
    std::vector<double> density_interp(sal.size());

    double temp;
    int index, mask_index;
    #if DEBUG >= 0
    fprintf(stdout, "Computing density and pressure at water points before interpolating.\n");
    #endif
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask.at(mask_index) == 1) {
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
            }
        }
    }

    // Apply interpolation
    #if DEBUG >= 0
    fprintf(stdout, "Interpolating density\n");
    #endif
    interpolate_over_land(density_interp, density, time, depth, latitude, longitude, mask);

    #if DEBUG >= 0
    fprintf(stdout, "Interpolating pressure\n");
    #endif
    interpolate_over_land(pressure_interp, pressure, time, depth, latitude, longitude, mask);

    // Now output to interp.nc
    #if DEBUG >= 0
    fprintf(stdout, "Creating output\n");
    #endif
    initialize_output_file(time, depth, longitude, latitude, scales, mask, "interp.nc");

    #if DEBUG >= 0
    fprintf(stdout, "Adding appropraite variables to source.nc\n");
    #endif
    const char* dim_names[] = {"time", "depth", "latitude", "longitude"};
    add_var_to_file("pressure_orig", dim_names, 4, "interp.nc");
    add_var_to_file("density_orig",  dim_names, 4, "interp.nc");

    add_var_to_file("pressure_interp", dim_names, 4, "interp.nc");
    add_var_to_file("density_interp",    dim_names, 4, "interp.nc");

    size_t starts[] = {0, 0, 0, 0};
    size_t counts[] = {time.size(),
                       depth.size(),
                       latitude.size(),
                       longitude.size()
                      };

    #if DEBUG >= 0
    fprintf(stdout, "Writing appropriate variables to source.nc\n");
    #endif
    write_field_to_output(pressure, "pressure_orig", starts, counts, "interp.nc");
    write_field_to_output(density,  "density_orig",  starts, counts, "interp.nc");

    write_field_to_output(pressure_interp, "pressure_interp", starts, counts, "interp.nc");
    write_field_to_output(density_interp,  "density_interp",  starts, counts, "interp.nc");

    return 0;
}

