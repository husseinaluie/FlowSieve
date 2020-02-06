#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"

void compute_vorticity(
        std::vector<double> & vort_r,           /**< [in] where to store vort_r */
        std::vector<double> & vort_lon,         /**< [in] where to store vort_lon */
        std::vector<double> & vort_lat,         /**< [in] where to store vort_lat */
        const std::vector<double> & u_r,        /**< [in] full u_r for calculation */
        const std::vector<double> & u_lon,      /**< [in] full u_lon for calculation */
        const std::vector<double> & u_lat,      /**< [in] full u_lat for calculation */
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon,                         /**< [in] Length of longitude dimension */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask,       /**< [in] Mask array (2D) to distinguish land from water */
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    double vort_r_tmp, vort_lon_tmp, vort_lat_tmp;
    int Itime, Idepth, Ilat, Ilon;
    size_t index; 
    const size_t Npts = u_r.size();

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  Starting vorticity computation.\n"); }
    #endif

    #pragma omp parallel \
    default(none) \
    shared(mask, u_r, u_lon, u_lat, longitude, latitude,\
            vort_r, vort_lon, vort_lat)\
    private(Itime, Idepth, Ilat, Ilon, index, \
            vort_r_tmp, vort_lon_tmp, vort_lat_tmp)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; index++) {

            vort_r_tmp   = constants::fill_value;
            vort_lon_tmp = constants::fill_value;
            vort_lat_tmp = constants::fill_value;

            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                compute_vorticity_at_point(
                        vort_r_tmp, vort_lon_tmp, vort_lat_tmp,
                        u_r,        u_lon,        u_lat,
                        Ntime, Ndepth, Nlat, Nlon,
                        Itime, Idepth, Ilat, Ilon,
                        longitude, latitude,
                        mask);
            }

            vort_r.at(  index) = vort_r_tmp;
            vort_lon.at(index) = vort_lon_tmp;
            vort_lat.at(index) = vort_lat_tmp;

        } // end index loop
    } // end pragma
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
} // end compute_vorticity
