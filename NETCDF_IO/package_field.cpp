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
        const std::vector<double> * mask) {

    int index, mask_index;
    fprintf(stdout, "    Getting mask size\n");
    const int mask_len = (int) mask->size();
    fprintf(stdout, "      size is %d\n", mask_len);

    #if DEBUG >= 2
    fprintf(stdout, "    Computing the data range.\n");
    fflush(stdout);
    #endif

    // First, we need to compute the min and max values 
    //   to allow us to convert to signed shorts
    double fmax=0, fmin=0;
    #pragma omp parallel \
    private(index) \
    reduction(max : fmax) reduction(min : fmin)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int) original.size(); index++) {
            mask_index = index % mask_len;
            if (mask->at(mask_index) == 1) {
                fmax = std::max(fmax, original.at(index));
                fmin = std::min(fmin, original.at(index));
            }
        }
    }

    // Number of Discrete Representable Values
    //   (less two for numerical reasons)
    int ndrv = pow(2, 16) - 2;

    #if DEBUG >= 2
    fprintf(stdout, "    Rescaling the data (fmin, fmax) = (%g, %g)\n", fmin, fmax);
    fflush(stdout);
    #endif

    // Now that we have the min/max, we go ahead and do the conversion
    double fmean  = 0.5 * (fmax + fmin);
    double frange = fmax - fmin;
    double local_double;
    signed short local_int;
    #pragma omp parallel \
    private(index, local_double, local_int) \
    shared(frange, fmean, ndrv)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int) original.size(); index++) {
            mask_index = index % mask_len;
            if (mask->at(mask_index) == 1) {

                // Scale original down to [-0.5,0.5] and store in local_double
                local_double = (original.at(index) - fmean) / frange;

                // Now convert to int in [-ndrv/2, ndrv/2]
                local_int = (signed short) round(ndrv * local_double);

            } else {
                local_int = constants::fill_value_s;
            }
            packaged.at(index) = local_int;
        }
    }

    scale_factor = frange / ndrv;
    add_offset = fmean;
}
