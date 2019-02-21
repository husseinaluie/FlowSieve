#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP 1

#include <stdio.h>
#include <stdlib.h>

void compute_areas(
        double * areas, 
        double * longitude, 
        double * latitude, 
        int Nlon, 
        int Nlat);

int Index( const int Itime, const int Idepth, const int Ilat, const int Ilon,
           const int Ntime, const int Ndepth, const int Nlat, const int Nlon  );

double distance(const double lon1, const double lat1, 
                const double lon2, const double lat2);

void vel_Spher_to_Cart(
            double & u_x, double & u_y, double & u_z,
            const double u_r, const double u_lon, const double u_lat,
            const double lon, const double lat );

void vel_Cart_to_Spher(
            double & u_r, double & u_lon, double & u_lat,
            const double u_x, const double u_y, const double u_z,
            const double lon, const double lat );

void filtering(const double * u_r, const double * u_lon, const double * u_lat,
               const double * scales, const int Nscales,
               const double dlon, const double dlat,
               const int Ntime, const int Ndepth, const int Nlon, const int Nlat,
               const double * dAreas, 
               const double * time, const double * depth,
               const double * longitude, const double * latitude);

void apply_filter_at_point(
        double & u_x_tmp,   double & u_y_tmp,   double & u_z_tmp,
        const double * u_x, const double * u_y, const double * u_z,
        const int dlon_N, const int dlat_N, 
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const double * longitude, const double * latitude,
        const double * dAreas, const double scale);

double kernel(const double distance, const double scale);

#endif
