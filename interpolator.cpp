#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "netcdf_io.hpp"
#include "preprocess.hpp"
#include "constants.hpp"

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

    std::vector<double> thetao_interp;
    std::vector<double> sal_interp;
    std::vector<double> ssh_interp;

    std::vector<double> mask;

    read_var_from_file(longitude, "longitude", "source.nc", NULL);
    read_var_from_file(latitude,  "latitude",  "source.nc", NULL);
    read_var_from_file(time,      "time",      "source.nc", NULL);
    read_var_from_file(depth,     "depth",     "source.nc", NULL);

    read_var_from_file(u_lat, "vo", "source.nc", &mask);

    read_var_from_file(thetao, "thetao", "source.nc", NULL);
    read_var_from_file(sal ,   "so",     "source.nc", NULL);
    read_var_from_file(ssh,    "zos",    "source.nc", NULL);

    // Convert to radians
    for (size_t Ilat = 0; Ilat < latitude.size(); Ilat++) {
        latitude.at(Ilat) = latitude.at(Ilat) * M_PI / 180.;
    }
    for (size_t Ilon = 0; Ilon < longitude.size(); Ilon++) {
        longitude.at(Ilon) = longitude.at(Ilon) * M_PI / 180.;
    }

    #if DEBUG >= 0
    fprintf(stdout, "Interpolating potential temperature (thetao)\n");
    #endif
    interpolate_over_land(thetao_interp, thetao, time, depth, latitude, longitude, mask);
    
    #if DEBUG >= 0
    fprintf(stdout, "Interpolating salinity (so)\n");
    #endif
    interpolate_over_land(sal_interp, sal, time, depth, latitude, longitude, mask);
    
    #if DEBUG >= 0
    fprintf(stdout, "Interpolating sea surface height (zos)\n");
    #endif
    interpolate_over_land(ssh_interp, ssh, time, depth, latitude, longitude, mask);
    
    #if DEBUG >= 0
    fprintf(stdout, "Creating output\n");
    #endif
    initialize_output_file(time, depth, longitude, latitude, scales, mask, "interp.nc");

    #if DEBUG >= 0
    fprintf(stdout, "Adding appropraite variables to source.nc\n");
    #endif
    const char* dim_names[] = {"time", "depth", "latitude", "longitude"};
    add_var_to_file("thetao_orig", dim_names, 4, "interp.nc");
    add_var_to_file("sal_orig",    dim_names, 4, "interp.nc");
    add_var_to_file("ssh_orig",    dim_names, 4, "interp.nc");

    add_var_to_file("thetao_interp", dim_names, 4, "interp.nc");
    add_var_to_file("sal_interp",    dim_names, 4, "interp.nc");
    add_var_to_file("ssh_interp",    dim_names, 4, "interp.nc");

    size_t starts[] = {0, 0, 0, 0};
    size_t counts[] = {time.size(),
                       depth.size(),
                       latitude.size(),
                       longitude.size()
                      };

    #if DEBUG >= 0
    fprintf(stdout, "Writing appropriate variables to source.nc\n");
    #endif
    write_field_to_output(thetao, "thetao_orig", starts, counts, "interp.nc");
    write_field_to_output(sal,    "sal_orig",    starts, counts, "interp.nc");
    write_field_to_output(ssh,    "ssh_orig",    starts, counts, "interp.nc");

    write_field_to_output(thetao_interp, "thetao_interp", starts, counts, "interp.nc");
    write_field_to_output(sal_interp,    "sal_interp",    starts, counts, "interp.nc");
    write_field_to_output(ssh_interp,    "ssh_interp",    starts, counts, "interp.nc");
    

    /*
    // Now that we've interpolated them, let's compute the density and pressure at each point
    std::vector<double> pressure(sal.size());
    std::vector<double> density(sal.size());

    double temp;
    for (size_t II = 0; II < sal.size();, II++) {
        // Compute pressure (geostrophic)
        pressure.at(II) = constants::rho0 * constants::g * ssh_interp.at(II);

        // Convert potential temp to actual temp
        temp = depotential_temperature(pressure.at(II), thetao.at(II));

        // Compute density from UNESCO equation of state
        density.at(II) = equation_of_state(temp, 
                sal_interp.at(II), 
                pressure.at(II) / 1e5); // takes pressure in bar = 100 kPa = 1e5 Pa
    }

    add_var_to_file("pressure", dim_names, 4, "interp.nc");
    add_var_to_file("density",  dim_names, 4, "interp.nc");

    write_field_to_output(pressure, "pressure", starts, counts, "interp.nc");
    write_field_to_output(density,  "density",  starts, counts, "interp.nc");
    */

    return 0;
}

