#include <algorithm>
#include <stdlib.h>
#include <omp.h>

/*! 
 * \brief Compute openmp chunk size for OMP for loops
 *
 * The point of using a function for this is that it standardizes
 *    the chunk size across the different functions.
 *
 * This assumes that the loop being split is a Lat-Lon loop.
 *
 * @param[in]   Nlat,Nlon   size of the latitude and longitude grids
 *
 * @returns chunksize that can be passed to OpenMP pragma commands
 *
 */
int get_omp_chunksize(const int Nlat, const int Nlon) {

    int nthreads = omp_get_num_threads();

    // Set the chunk size to be 1% of what each thread would get 
    //    in a uniform distribution
    int chunksize = std::max(1, (int) (Nlat * Nlon / 100. / nthreads));

    return chunksize;

}
