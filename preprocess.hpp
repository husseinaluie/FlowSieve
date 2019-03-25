#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP 1

/*!
 *  \brief Apply ALGLIB interpolation routines to fill in masked areas
 */
void interpolate_over_land(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<double> &mask);

void interpolate_over_land_from_coast(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<double> &mask);

void get_coast(
        std::vector<double> &lon_coast,
        std::vector<double> &lat_coast,
        std::vector<double> &field_coast,
        const std::vector<double> &lon_full,
        const std::vector<double> &lat_full,
        const std::vector<double> &field_full,
        const std::vector<double> &mask,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon);

#endif
