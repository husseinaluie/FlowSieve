#include <vector>
#include <math.h>
#include <omp.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void package_field(
        std::vector<signed short> & packaged,
        double & scale_factor,
        double & add_offset,
        const std::vector<double> & original,
        const std::vector<bool> * mask,
        MPI_Comm comm
        ) {

    size_t index;

    // First, we need to compute the min and max values 
    //   to allow us to convert to signed shorts
    double fmax, fmin, fmin_loc=0, fmax_loc=0;
    #pragma omp parallel \
    default(none) private(index) \
    shared(original, mask) \
    reduction(max : fmax_loc) reduction(min : fmin_loc)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < original.size(); index++) {
            if ( mask->at(index) ) {
                fmax_loc = std::max(fmax_loc, original.at(index));
                fmin_loc = std::min(fmin_loc, original.at(index));
            }
        }
    }
    MPI_Allreduce(&fmax_loc, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&fmin_loc, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);

    // Number of Discrete Representable Values
    //   (less two for numerical reasons)
    //int ndrv = pow(2, 16) - 2;
    const int ndrv =   constants::fill_value_s < 0 
                     ? constants::fill_value_s + 2 
                     : constants::fill_value_s - 2;

    // Now that we have the min/max, we go ahead and do the conversion
    const double fmiddle = 0.5 * (fmax + fmin);
    const double frange  = fmax - fmin;
    double local_double;
    signed short local_int;
    #pragma omp parallel \
    default (none) \
    shared(mask, original, packaged) \
    private(index, local_double, local_int)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < original.size(); index++) {
            if ( mask->at(index) ) {

                // Scale original down to [-0.5,0.5] and store in local_double
                local_double = (original.at(index) - fmiddle) / frange;

                // Now convert to int in [-ndrv/2, ndrv/2]
                local_int = (signed short) round(ndrv * local_double);

            } else {
                local_int = constants::fill_value_s;
            }
            packaged.at(index) = local_int;
        }
    }

    scale_factor = frange / ndrv;
    add_offset = fmiddle;
}
