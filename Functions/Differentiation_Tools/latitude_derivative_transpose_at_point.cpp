#include <vector>
#include <string>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

void latitude_derivative_transpose_at_point(
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
                            and ( (Ilat <= 4) or (Ilat >= Nlat - 5) ) ;

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

            double deriv_val = 0.;

            if ( Ilat == 0 ) {
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1);

                size_t idx_p1 = source_data.local_index(Itime, Idepth, 1, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, 2, Ilon);

                deriv_val = -( (-15+8*5)*za0 + -5*za1 + 8.*F.at(idx_p1) - F.at(idx_p2) ) / (12. * dlat);

            } else if ( Ilat == 1 ) {
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1);

                size_t idx_0  = source_data.local_index(Itime, Idepth, 0, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, 2, Ilon);
                size_t idx_p3 = source_data.local_index(Itime, Idepth, 3, Ilon);

                deriv_val = -( (40-8*10)*za0 + 10*za1 - 8.*F.at(idx_0) + 8.*F.at(idx_p2) - F.at(idx_p3) ) / (12. * dlat);

            } else if ( Ilat == 2 ) {
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1);

                size_t idx_p0 = source_data.local_index(Itime, Idepth, 0, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, 1, Ilon);
                size_t idx_p3 = source_data.local_index(Itime, Idepth, 3, Ilon);
                size_t idx_p4 = source_data.local_index(Itime, Idepth, 4, Ilon);

                deriv_val = -( (-45+8*10)*za0 + -10*za1 + F.at(idx_p0) - 8.*F.at(idx_p1) + 8.*F.at(idx_p3) - F.at(idx_p4) ) / (12. * dlat);

            } else if ( Ilat == 3 ) {
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1);

                size_t idx_p0 = source_data.local_index(Itime, Idepth, 1, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, 2, Ilon);
                size_t idx_p3 = source_data.local_index(Itime, Idepth, 4, Ilon);
                size_t idx_p4 = source_data.local_index(Itime, Idepth, 5, Ilon);

                deriv_val = -( (24-8*5)*za0 + 5*za1 + F.at(idx_p0) - 8.*F.at(idx_p1) + 8.*F.at(idx_p3) - F.at(idx_p4) ) / (12. * dlat);

            } else if ( Ilat == 4 ) {
                const double za0 = zonal_mean_at_row(0),
                             za1 = zonal_mean_at_row(1);

                size_t idx_p0 = source_data.local_index(Itime, Idepth, 2, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, 3, Ilon);
                size_t idx_p3 = source_data.local_index(Itime, Idepth, 5, Ilon);
                size_t idx_p4 = source_data.local_index(Itime, Idepth, 6, Ilon);

                deriv_val = -( (-5+8)*za0 + -1*za1 + F.at(idx_p0) - 8.*F.at(idx_p1) + 8.*F.at(idx_p3) - F.at(idx_p4) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 1 ) {
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2);

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Nlat-3, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Nlat-2, Ilon);

                deriv_val = -( (15-8*5)*za0 + 5*za1 +  F.at(idx_m2) - 8.*F.at(idx_m1) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 2 ) {
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2);

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Nlat-4, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Nlat-3, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, Nlat-1, Ilon);

                deriv_val = -( (-40+8*10)*za0 + -10*za1 + F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*F.at(idx_p1) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 3 ) {
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2);

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Ilat-2, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Ilat-1, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, Ilat+1, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, Ilat+2, Ilon);

                deriv_val = -( ( 45-8*10)*za0 +  10*za1 + F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*F.at(idx_p1) - F.at(idx_p2) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 4 ) {
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2);

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Ilat-2, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Ilat-1, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, Ilat+1, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, Ilat+2, Ilon);

                deriv_val = -( (-24+8*5)*za0 + -5*za1 + F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*F.at(idx_p1) - F.at(idx_p2) ) / (12. * dlat);

            } else if ( Ilat == Nlat - 5 ) {
                const double za0 = zonal_mean_at_row(Nlat-1),
                             za1 = zonal_mean_at_row(Nlat-2);

                size_t idx_m2 = source_data.local_index(Itime, Idepth, Ilat-2, Ilon);
                size_t idx_m1 = source_data.local_index(Itime, Idepth, Ilat-1, Ilon);
                size_t idx_p1 = source_data.local_index(Itime, Idepth, Ilat+1, Ilon);
                size_t idx_p2 = source_data.local_index(Itime, Idepth, Ilat+2, Ilon);

                deriv_val = -( (5-8)*za0 + 1*za1 + F.at(idx_m2) - 8.*F.at(idx_m1) + 8.*F.at(idx_p1) - F.at(idx_p2) ) / (12. * dlat);
            }

            *(deriv_vals.at(ii)) = deriv_val;
        }

        return;
    }

    // Otherwise, just a regular derivative
    spher_derivative_transpose_at_point( deriv_vals, fields, source_data.latitude, "lat",
            source_data, Itime, Idepth, Ilat, Ilon, mask, order_of_deriv, diff_ord );

}
