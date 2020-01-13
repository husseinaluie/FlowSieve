#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

double toroidal_comp_RHS(  
        std::vector<double> & out_arr,
        std::vector<double> & work_arr_1,
        std::vector<double> & work_arr_2,
        std::vector<double> & Lap_F,
        const std::vector<double> & F,
        const std::vector<double> & curl_term,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask,
        const std::vector<double> & areas
    ) {

    /*
    def Lap_F(f, R = R_earth):
    ret =   (1/(R*cos_lat)) * ddlon( (1/(R*cos_lat)) * ddlon(f) ) \
          + (1/(R*cos_lat)) * ddlat( cos_lat / R * ddlat(f) )
    return ret
    */

    int Itime, Idepth, Ilat, Ilon, index;
    double dFdlon, dFdlat, cos_lat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&F);

    #pragma omp parallel \
    default(none) \
    shared( work_arr_1, work_arr_2, latitude, longitude, mask, F, deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, \
            dFdlon, dFdlat, lon_deriv_vals, lat_deriv_vals)
    {

        lon_deriv_vals.push_back(&dFdlon);
        lat_deriv_vals.push_back(&dFdlat);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)F.size(); index++) {

            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields,
                        longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields,
                        latitude, "lat",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                cos_lat = cos(latitude.at(Ilat));

                work_arr_1.at(index) = dFdlon / ( cos_lat * constants::R_earth );
                work_arr_2.at(index) = dFdlat *   cos_lat / constants::R_earth;
            }
        }
    }


    deriv_fields.clear();
    lon_deriv_vals.clear();
    lat_deriv_vals.clear();

    deriv_fields.push_back(&work_arr_1);
    deriv_fields.push_back(&work_arr_2);
    double lon_deriv, lat_deriv, tmp;

    #pragma omp parallel \
    default(none) \
    shared( work_arr_1, work_arr_2, F, Lap_F, curl_term, \
            latitude, longitude, mask, deriv_fields, out_arr)\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, tmp, \
            lon_deriv, lat_deriv, lon_deriv_vals, lat_deriv_vals)
    {

        lon_deriv_vals.push_back(&lon_deriv);
        lon_deriv_vals.push_back(NULL);

        lat_deriv_vals.push_back(NULL);
        lat_deriv_vals.push_back(&lat_deriv);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)F.size(); index++) {

            tmp = constants::fill_value;

            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields,
                        longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields,
                        latitude, "lat",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                cos_lat = cos(latitude.at(Ilat));
                tmp = ( lon_deriv + lat_deriv ) / ( constants::R_earth * cos_lat );
            }

            Lap_F.at(index) = tmp;
            out_arr.at(index) = Lap_F.at(index) - curl_term.at(index);

        }
    }

    double RHS_norm = 0., net_area = 0.;
    int area_index;

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, F, mask, deriv_fields, areas, out_arr ) \
    private(Itime, Idepth, Ilat, Ilon, index, area_index) \
    reduction(+ : RHS_norm, net_area)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)F.size(); index++) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            area_index = Index(0, 0, Ilat, Ilon,
                               1, 1, Nlat, Nlon);

            if (mask.at(index) == 1) { // Skip land areas
                RHS_norm += pow(out_arr.at(index), 2) * areas.at(area_index);
                net_area += areas.at(area_index);
            }
        }
    }

    return sqrt(RHS_norm / net_area);
}

