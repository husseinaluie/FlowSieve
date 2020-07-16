#include <algorithm>
#include <stdlib.h>
#include <omp.h>

int get_omp_chunksize(const int Nlat, const int Nlon) {

    int nthreads = omp_get_num_threads();

    // Set the chunk size to be 1% of what each thread would get 
    //    in a uniform distribution
    int chunksize = std::max(1, (int) (Nlat * Nlon / 100. / nthreads));

    return chunksize;

}
