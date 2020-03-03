#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

void KE_from_vels(
            std::vector<double> & KE,
            std::vector<double> * u1,
            std::vector<double> * u2,
            std::vector<double> * u3,
            const std::vector<double> & mask,
            const double rho0
        ) {

    size_t index;
    double tmp;

    #pragma omp parallel default(none) private(index, tmp) shared(KE, u1, u2, u3, mask)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u1->size(); ++index) {

            if (mask.at(index) == 1) { // Skip land areas
                tmp = 0.;
                if (u1 != NULL) { tmp += 0.5 * rho0 * pow(u1->at(index), 2); }
                if (u2 != NULL) { tmp += 0.5 * rho0 * pow(u2->at(index), 2); }
                if (u3 != NULL) { tmp += 0.5 * rho0 * pow(u3->at(index), 2); }
                KE.at(index) = tmp;
            } else {
                KE.at(index) = constants::fill_value;
            }
        }
    }
}
