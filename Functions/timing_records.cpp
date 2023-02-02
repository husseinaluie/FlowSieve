#include <math.h>
#include <algorithm>
#include <mpi.h>
#include <cassert>
#include "../constants.hpp"
#include "../functions.hpp"

// This file provides the implementation details for the Timing_Records class

// Class constructor
Timing_Records::Timing_Records() {
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
//    then create a new record for that name
void Timing_Records::add_to_record( const double delta, const std::string record_name ) {
    if (time_records.count(record_name) == 1) {
        // If the record already exists, add on to it
        time_records[record_name] += delta;
    } else {
        // Otherwise, just create it
        time_records.insert( std::pair< std::string, double >( record_name, delta ) );
    }
}

// Print the results.
//    Compute mean and standard deviations (across processors)
//    and print the results. 
void Timing_Records::print() const{

    int wRank = -1, wSize = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    double time_val, mean_val, std_val, tmp, total_time = 0., total_variance = 0.;
    if (wRank == 0) {
        fprintf(stdout, "\n\n## Internal Timings : mean ( standard deviation )\n\n");
    }

    for(const auto& entry : time_records) {
        time_val = entry.second;

        // Get mean timing value
        //MPI_Allreduce( &time_val, &mean_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Reduce( &time_val, &mean_val, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        mean_val = mean_val / wSize;

        // Get standard deviation timing value
        tmp = pow(time_val - mean_val, 2);

        //MPI_Allreduce( &tmp, &std_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Reduce( &tmp, &std_val, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        std_val = pow(std_val, 0.5) / wSize;

        // Also accumulate total values
        total_time += mean_val;
        total_variance += pow( std_val, 2.);

        // Print information
        if (wRank == 0) {
            fprintf(stdout, "  %-35s : %8.4e ( %8.4e )\n", entry.first.c_str(), mean_val, std_val);
        }
    }
    if (wRank == 0) {
        fprintf(stdout, " --------------------------------------- \n" );
        fprintf(stdout, "  Total : %8.4e ( %8.4e )\n", total_time, sqrt( total_variance ) );
    }
}
