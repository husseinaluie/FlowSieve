#include <vector>
#include <string>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

void latitude_derivative_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const std::vector<short int> & mask,
        const bool polar_smoothing,
        const int order_of_deriv,
        const int diff_ord
        ){

    const int Nlat = source_data.Nlat;
    const int Nlon = source_data.Nlon;

    // Zonal-mean polar extrapolation is only implemented for the standard
    // 4th-order first-derivative stencil 
    const bool near_pole = polar_smoothing
                            and (order_of_deriv == 1)
                            and (diff_ord == 4)
                            and ( (Ilat == 0) or (Ilat == 1) or
                                  (Ilat == Nlat - 2) or (Ilat == Nlat - 1) );

    if ( near_pole ) {

        const int num_fields = fields.size();

        // Zero out before computing
        for (int ii = 0; ii < num_fields; ii++) {
            if (deriv_vals.at(ii) != NULL) {
                *(deriv_vals.at(ii)) = 0.;
            }
        }

        // dlat assumed uniform
        const std::vector<double> &latitude = source_data.latitude;
        const double dlat = latitude.at(1) - latitude.at(0);

        for (int ii = 0; ii < num_fields; ii++) {

            if (deriv_vals.at(ii) == NULL) { continue; }

            const std::vector<double> &F = *fields.at(ii);

            // ── Compute zonal means at the five latitude rows nearest
            //    the relevant pole. We only need the rows that actually
            //    feed into the stencil at this Ilat, but it's cheapest
            //    and clearest to just compute the full 5-row set used
            //    by the Python reference for whichever pole we're near.
            //
            //    za[0..4] correspond to:
            //      south pole side: rows 0,1,2,3,4
            //      north pole side: rows Nlat-1, Nlat-2, Nlat-3, Nlat-4, Nlat-5
            // ──────────────────────────────────────────────────────────
            auto zonal_mean_at_row = [&](int row) -> double {
                double sum = 0.;
                for (int Jlon = 0; Jlon < Nlon; Jlon++) {
                    size_t idx = source_data.local_index(Itime, Idepth, row, Jlon);
                    sum += F.at(idx);
                }
                return sum / (double) Nlon;
            };

            // Ghost-row values (zonal-mean extrapolations), matching the
            // Python reference's _ddlat_tmp[0], [1], [-2], [-1]:
            //
            //   _ddlat_tmp[0]  = 15*za[0] - 40*za[1] + 45*za[2] - 24*za[3] +  5*za[4]
            //   _ddlat_tmp[1]  =  5*za[0] - 10*za[1] + 10*za[2] -  5*za[3] +  1*za[4]
            //   _ddlat_tmp[-2] =  5*za[-1] - 10*za[-2] + 10*za[-3] -  5*za[-4] +  1*za[-5]
            //   _ddlat_tmp[-1] = 15*za[-1] - 40*za[-2] + 45*za[-3] - 24*za[-4] +  5*za[-5]
            //
            // The stencil for the derivative at row Ilat is the standard
            // 4th-order central difference:
            //   deriv = ( g[Ilat-2] - 8*g[Ilat-1] + 8*g[Ilat+1] - g[Ilat+2] ) / (12*dlat)
            // where g[] is the field padded with the above ghost rows at
            // indices -2,-1 (south side) and Nlat,Nlat+1 (north side).

            double deriv_val = 0.;

            if ( Ilat == 0 ) {
                // Need g[-2], g[-1], g[1], g[2]
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1),
                             za2 = zonal_mean_at_row(2),
                             za3 = zonal_mean_at_row(3),
                             za4 = zonal_mean_at_row(4);

                const double g_m2 = 15.*za0 - 40.*za1 + 45.*za2 - 24.*za3 +  5.*za4; // ghost row -2
                const double g_m1 =  5.*za0 - 10.*za1 + 10.*za2 -  5.*za3 +  1.*za4; // ghost row -1

                size_t idx_p1 = source_data.local_index(Itime, Idepth, 1, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, 2, Ilon);

                deriv_val = ( g_m2 - 8.*g_m1 + 8.*F.at(idx_p1) - F.at(idx_p2) ) / (12. * dlat);

            } else if ( Ilat == 1 ) {
                // Need g[-1], g[0], g[2], g[3]
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1),
                             za2 = zonal_mean_at_row(2),
                             za3 = zonal_mean_at_row(3),
                             za4 = zonal_mean_at_row(4);

                const double g_m1 =  5.*za0 - 10.*za1 + 10.*za2 -  5.*za3 +  1.*za4; // ghost row -1

                size_t idx_0  = source_data.local_index(Itime, Idepth, 0, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, 2, Ilon);
                size_t idx_p3 = source_data.local_index(Itime, Idepth, 3, Ilon);

                deriv_val = ( g_m1 - 8.*F.at(idx_0) + 8.*F.at(idx_p2) - F.at(idx_p3) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 2 ) {
                // Central stencil at row Nlat-2 needs:
                //   F[Nlat-4]  (= row-2), F[Nlat-3] (= row-1),
                //   F[Nlat-1]  (= row+1), g[Nlat]   (= row+2, one ghost row past the last real row)
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2),
                             za2 = zonal_mean_at_row(Nlat-3),
                             za3 = zonal_mean_at_row(Nlat-4),
                             za4 = zonal_mean_at_row(Nlat-5);

                const double g_p1 =  5.*za0 - 10.*za1 + 10.*za2 -  5.*za3 +  1.*za4; // ghost row Nlat (one past last)

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Nlat-4, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Nlat-3, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, Nlat-1, Ilon);

                deriv_val = ( F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*F.at(idx_p1) - g_p1 ) / (12. * dlat);

            } else if ( Ilat == Nlat - 1 ) {
                // Need g[Nlat-3], g[Nlat-2], g[Nlat], g[Nlat+1]  (ghost rows one and two past last)
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2),
                             za2 = zonal_mean_at_row(Nlat-3),
                             za3 = zonal_mean_at_row(Nlat-4),
                             za4 = zonal_mean_at_row(Nlat-5);

                const double g_p1 =  5.*za0 - 10.*za1 + 10.*za2 -  5.*za3 +  1.*za4; // ghost row Nlat
                const double g_p2 = 15.*za0 - 40.*za1 + 45.*za2 - 24.*za3 +  5.*za4; // ghost row Nlat+1

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Nlat-3, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Nlat-2, Ilon);

                deriv_val = ( F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*g_p1 - g_p2 ) / (12. * dlat);
            }

            *(deriv_vals.at(ii)) = deriv_val;
        }

        return;
    }

    // Otherwise, just a regular derivative
    spher_derivative_at_point( deriv_vals, fields, source_data.latitude, "lat",
            source_data, Itime, Idepth, Ilat, Ilon, mask, order_of_deriv, diff_ord );

}
