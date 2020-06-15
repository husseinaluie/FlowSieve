#include <math.h>
#include <algorithm>
#include <mpi.h>
#include <cassert>
#include "../constants.hpp"
#include "../functions.hpp"

// This file provides the implementation details for the Timing_Records class

// Class constructor
//   create a key-value pair for each desired record
//   initialize them to zero
//   order here does not matter
Timing_Records::Timing_Records() {

    time_records["kernel_precomputation"] = 0.;

    time_records["filter_main"] = 0.;
    time_records["filter_for_Pi"] = 0.;
    time_records["filter_for_Lambda"] = 0.;
    time_records["land"] = 0.;

    time_records["filter_tor_pot"] = 0.;

    time_records["compute_vorticity"] = 0.;
    time_records["compute_Pi"] = 0.;
    time_records["compute_Lambda"] = 0.;
    time_records["compute_transport"] = 0.;

    time_records["writing"] = 0.;

    time_records["postprocess"] = 0.;

}

// Reset all records to zero.
//   This will autotically hit all records 
void Timing_Records::reset() {
    for(auto& entry : time_records) {
        entry.second = 0.;
    }
}

// Update a record with a new time delta
//    Given a delta and a record_name,
//    update the record indicated by record_name
//    by adding delta to the previous value.
//
//    If record_name does not map to a valid record,
//    then print an error statement and halt.
void Timing_Records::add_to_record( double delta, std::string record_name ) {

    int wRank = -1;

    if (time_records.count(record_name) == 1) {
        time_records[record_name] += delta;
    } else {
        MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
        if (wRank == 0) {
            fprintf(stderr, "%s is not a valid key of time_records."
                    " Aborting\n", record_name.c_str());
        }
        assert(false);
    }
}

// Print the results.
//    Compute mean and standard deviations (across processors)
//    and print the results. 
void Timing_Records::print() {

    int wRank = -1, wSize = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    double time_val, mean_val, std_val, tmp;
    if (wRank == 0) {
        fprintf(stdout, "\n\n## Internal Timings : mean ( standard deviation )\n\n");
    }

    for(const auto& entry : time_records) {
        time_val = entry.second;

        // Get mean timing value
        MPI_Allreduce( &time_val, &mean_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        mean_val = mean_val / wSize;

        // Get standard deviation timing value
        tmp = pow(time_val - mean_val, 2);

        MPI_Allreduce( &tmp, &std_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std_val = pow(std_val, 0.5) / wSize;

        // Print information
        if (wRank == 0) {
            fprintf(stdout, "  %-25s : %8.6e ( %8.6e )\n", 
                    entry.first.c_str(), mean_val, std_val);
        }
    }
}
