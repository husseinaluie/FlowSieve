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

#endif
