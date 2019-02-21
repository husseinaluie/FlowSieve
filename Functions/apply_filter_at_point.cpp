#include <algorithm>
#include "../functions.hpp"


void apply_filter_at_point(
        double & u_x_tmp, double & u_y_tmp, double & u_z_tmp,
        const double * u_x, const double * u_y, const double * u_z,
        const int dlon_N, const int dlat_N, 
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const double * longitude, const double * latitude,
        const double * dAreas, const double scale) {


    double kA_sum, dist, kern, area;
    int index;

    kA_sum  = 0.;
    u_x_tmp = 0.;
    u_y_tmp = 0.;
    u_z_tmp = 0.;

    for (int LAT = std::max(0, Ilat - dlat_N); LAT < Ilat + dlat_N; LAT++) {
        for (int LON = std::max(0, Ilon - dlon_N); LON < Ilon + dlon_N; LON++) {

            dist = distance(longitude[Ilon], latitude[Ilat],
                            longitude[LON],  latitude[LAT]);

            kern = kernel(dist, scale);

            index = Index(Itime, Idepth, LAT,  LON,
                          Ntime, Ndepth, Nlat, Nlon);

            area    = dAreas[index];
            kA_sum += kern * area;

            u_x_tmp += u_x[index] * kern * area;
            u_y_tmp += u_y[index] * kern * area;
            u_z_tmp += u_z[index] * kern * area;

            //fprintf(stdout, "           (u_x, u_y, u_z) = (%.4g, %.4g, %.4g)\n", u_x[index], u_y[index], u_z[index]);
        }
    }

    u_x_tmp = u_x_tmp / kA_sum;
    u_y_tmp = u_y_tmp / kA_sum;
    u_z_tmp = u_z_tmp / kA_sum;
}

