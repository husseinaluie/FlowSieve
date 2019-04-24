#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

void compute_distances(
        std::vector<double> & distances,        /**< [in] where to store the distances */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const int ref_ilat,
        const int ref_ilon,
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon                          /**< [in] Length of longitude dimension */
        ){

    double dist;
    int index;

    #if CARTESIAN
    // Grid spacing: assume uniform grid
    const double dlat_m = latitude.at( 1) - latitude.at( 0);
    const double dlon_m = longitude.at(1) - longitude.at(0);
    #endif

    const double ref_lat = latitude.at(ref_ilat);
    const double ref_lon = longitude.at(ref_ilon);

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            
            index = Index(0,     0,      Ilat, Ilon,
                          Ntime, Ndepth, Nlat, Nlon);

            #if CARTESIAN
            dist = distance(longitude.at(Ilon), latitude.at(Ilat),
                            ref_lon, ref_lat,
                            dlon_m * Nlon, dlat_m * Nlat);
            #elif not(CARTESIAN)
            dist = distance(longitude.at(Ilon), latitude.at(Ilat),
                            ref_lon, ref_lat);
            #endif

            distances.at(index) = dist;
        } // end longitude loop
    } // end latitude loop
}
