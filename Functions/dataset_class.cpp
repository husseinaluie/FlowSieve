#include <math.h>
#include <algorithm>
#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include <cassert>


// Class constructor
dataset::dataset() {
};

void dataset::load_time( const std::string dim_name, const std::string filename ) {
    if ( ( dim_name == "DNE" ) or ( dim_name == "DOES_NOT_EXIST" ) ) {
        time.resize(1);
        time[0] = 0.;
        #if DEBUG >= 1
        int wRank=-1;
        MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
        if (wRank == 0) { fprintf(stdout, "Time dimension DNE, so setting as singleton.\n"); }
        #endif
    } else {
        read_var_from_file(time, dim_name, filename);
    }
    full_Ntime = time.size();
};

void dataset::load_depth( const std::string dim_name, const std::string filename ) {
    if ( ( dim_name == "DNE" ) or ( dim_name == "DOES_NOT_EXIST" ) ) {
        depth.resize(1);
        depth[0] = 0.;
        #if DEBUG >= 1
        int wRank=-1;
        MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
        if (wRank == 0) { 
            fprintf(stdout, "Depth dimension DNE, so setting as singleton.\n\n");
            fflush( stdout );
        }
        #endif
    } else {
        read_var_from_file(depth, dim_name, filename);
    }
    full_Ndepth = depth.size();

    // Now also check if our depth grid is increasing or decreasing
    //  i.e. is the bottom the first or last point
    if (full_Ndepth > 1) {
        depth_is_increasing = ( depth[0] < depth[1] ) ? true : false;
    }
};

void dataset::load_latitude( const std::string dim_name, const std::string filename ) {
    read_var_from_file(latitude, dim_name, filename);
    Nlat = latitude.size();
};

void dataset::load_longitude( const std::string dim_name, const std::string filename ) {
    read_var_from_file(longitude, dim_name, filename);
    Nlon = longitude.size();
};

void dataset::compute_cell_areas() {

    assert( (Nlat > 0) and (Nlon > 0) );    

    areas.resize( Nlat * Nlon );
    compute_areas( areas, longitude, latitude );
}

void dataset::load_variable( 
        const std::string var_name, 
        const std::string var_name_in_file, 
        const std::string filename,
        const bool read_mask,
        const bool load_counts,
        const bool do_splits
        ) {

    // Add a new entry to the variables dictionary with an empty array
    variables.insert( std::pair< std::string, std::vector<double> >( var_name, std::vector<double>() ) );

    // Now read in from the file and store in the variables dictionary
    read_var_from_file( variables.at( var_name ), var_name_in_file, filename, 
                        read_mask ? &mask : NULL, 
                        load_counts ? &myCounts : NULL, 
                        load_counts ? &myStarts : NULL, 
                        Nprocs_in_time, Nprocs_in_depth,
                        do_splits );
};

void dataset::check_processor_divisions(    const int Nprocs_in_time_input, 
                                            const int Nprocs_in_depth_input, 
                                            const MPI_Comm ) {

    assert( (full_Ntime > 0) and (full_Ndepth > 0) and (Nlon > 0) and (Nlat > 0) ); // Must read in dimensions before checking processor divisions.

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Apply some cleaning to the processor allotments if necessary. 
    Nprocs_in_time  = ( full_Ntime  == 1 ) ? 1 : 
                      ( full_Ndepth == 1 ) ? wSize : 
                                             Nprocs_in_time_input;
    Nprocs_in_depth = ( full_Ndepth == 1 ) ? 1 : 
                      ( full_Ntime  == 1 ) ? wSize : 
                                             Nprocs_in_depth_input;

    #if DEBUG >= 0
    if (Nprocs_in_time != Nprocs_in_time_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in time to %'d from %'d\n", Nprocs_in_time, Nprocs_in_time_input); }
    }
    if (Nprocs_in_depth != Nprocs_in_depth_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in depth to %'d from %'d\n", Nprocs_in_depth, Nprocs_in_depth_input); }
    }
    if (wRank == 0) { fprintf(stdout, " Nproc(time, depth) = (%'d, %'d)\n\n", Nprocs_in_time, Nprocs_in_depth); }
    #endif

    assert( Nprocs_in_time * Nprocs_in_depth == wSize );

    // Now that processor divisions have been tested, also create the sub-communicator items
    int color, key;

    // communicator for ranks with the same times
    color = wRank / Nprocs_in_depth;
    key   = wRank % Nprocs_in_depth;
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_sametimes); 
    #if DEBUG >= 2
    fprintf( stdout, "Processor %d has time-color %d and depth-color %d.\n", wRank, color, key );
    #endif

    // communicator for ranks with the same depths
    color = wRank % Nprocs_in_depth;
    key   = wRank / Nprocs_in_depth;
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_samedepths); 
}


void dataset::compute_region_areas() {

    #if DEBUG >= 2
    int wRank=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );

    if (wRank == 0) { fprintf(stdout, "  Computing geographic region areas\n"); fflush(stdout); }
    #endif

    assert( mask.size() > 0 ); // must read in mask before computing region areas
    assert( areas.size() > 0 ); // must compute cell areas before computing region areas
    const size_t num_regions = region_names.size();
    region_areas.resize( num_regions * Ntime * Ndepth );
    if (constants::FILTER_OVER_LAND) { region_areas_water_only.resize( region_areas.size() ); }

    double local_area, local_area_water_only;
    size_t Ilat, Ilon, reg_index, index, area_index;

    for (size_t Iregion = 0; Iregion < num_regions; ++Iregion) {
        for (size_t Itime = 0; Itime < (size_t) Ntime; ++Itime) {
            for (size_t Idepth = 0; Idepth < (size_t) Ndepth; ++Idepth) {

                local_area = 0.;
                local_area_water_only = 0.;

                #pragma omp parallel default(none)\
                private( Ilat, Ilon, index, area_index )\
                shared( mask, areas, Iregion, Itime, Idepth ) \
                reduction(+ : local_area, local_area_water_only)
                { 
                    #pragma omp for collapse(2) schedule(guided)
                    for (Ilat = 0; Ilat < (size_t) Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < (size_t) Nlon; ++Ilon) {

                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( mask.at(index) ) { // Skip land areas

                                area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                                if ( regions.at( region_names.at( Iregion ) ).at(area_index) ) {
                                    local_area += areas.at(area_index); 
                                }
                            }

                            if (constants::FILTER_OVER_LAND) {
                                if ( reference_mask.at(index) ) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                                    if ( regions.at( region_names.at( Iregion ) ).at(area_index) ) {
                                        local_area_water_only += areas.at(area_index); 
                                    }
                                }
                            }
                        }
                    }
                }
                reg_index = Index(0, Itime, Idepth, Iregion, 1, Ntime, Ndepth, num_regions);

                region_areas.at(reg_index) = local_area;
                if (constants::FILTER_OVER_LAND) { region_areas_water_only.at(reg_index) = local_area_water_only; }
                
            }
        }
    }
}




void dataset::gather_variable_across_depth( const std::vector<double> & var,
                                            std::vector<double> & gathered_var
                                            ) const {

    const MPI_Comm &comm = MPI_subcomm_sametimes;

    int wSize, wRank;
    MPI_Comm_size( comm, &wSize );
    MPI_Comm_rank( comm, &wRank );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Preparing to merge variables across depth. Pre-merge size is %'zu\n", var.size() ); }
    #endif

    // resize the new variable
    size_t gathered_size = Ntime * full_Ndepth * Nlat * Nlon;
    gathered_var.resize( gathered_size, 0. );
    
    // Get the Ndepth assignments for each rank
    std::vector<int> Ndepths( wSize ),
                     rec_counts( wSize, 1 ),
                     offsets( wSize );
    for ( int i = 0; i < wSize; ++i ){ offsets[i] = i; }
    MPI_Allgatherv( &Ndepth, 1, MPI_INT, &Ndepths[0], &rec_counts[0], &offsets[0], MPI_INT, comm );

    // Now set up for and do the actual gather
    int sendcount = Ntime * Ndepth * Nlat * Nlon;
    for ( int i = 0; i < wSize; ++i ){
        rec_counts.at(i) = Ntime * Ndepths.at(i) * Nlat * Nlon;
        offsets.at(i) = (i==0) ? 0 : (offsets.at(i-1) + rec_counts.at(i-1));
    }

    MPI_Allgatherv( &var[0],          sendcount,                   MPI_DOUBLE, 
                    &gathered_var[0], &rec_counts[0], &offsets[0], MPI_DOUBLE, 
                    comm );
}

void dataset::gather_mask_across_depth( const std::vector<bool> & var,
                                        std::vector<bool> & gathered_var
                                      ) const {

    const MPI_Comm &comm = MPI_subcomm_sametimes;

    int wSize, wRank;
    MPI_Comm_size( comm, &wSize );
    MPI_Comm_rank( comm, &wRank );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Preparing to merge mask across depth. Pre-merge size is %'zu\n", var.size() ); }
    #endif

    // resize the new variable
    size_t gathered_size = Ntime * full_Ndepth * Nlat * Nlon;
    gathered_var.resize( gathered_size, false );
    
    // Get the Ndepth assignments for each rank
    std::vector<int> Ndepths( wSize ),
                     rec_counts( wSize, 1 ),
                     offsets( wSize );
    for ( int i = 0; i < wSize; ++i ){ offsets[i] = i; }
    MPI_Allgatherv( &Ndepth, 1, MPI_INT, &Ndepths[0], &rec_counts[0], &offsets[0], MPI_INT, comm );

    // Now set up for and do the actual gather
    int sendcount = Ntime * Ndepth * Nlat * Nlon;
    for ( int i = 0; i < wSize; ++i ){
        rec_counts.at(i) = Ntime * Ndepths.at(i) * Nlat * Nlon;
        offsets.at(i) = (i==0) ? 0 : (offsets.at(i-1) + rec_counts.at(i-1));
    }

    // For **really** annoying reasons, std::vector<bool> is a special thing, and
    //      we can't do things like &var[0] like we would with a std::vector<double>.
    //      This means we need to make an int array first, and the communicate that.
    std::vector<int> src(var.size()), dest(gathered_var.size());
    for ( size_t index = 0; index < var.size(); index++ ) { src.at(index) = (int)var.at(index); }
    MPI_Allgatherv( &src[0],  sendcount,                   MPI_INT, 
                    &dest[0], &rec_counts[0], &offsets[0], MPI_INT, 
                    comm );
    for ( size_t index = 0; index < gathered_var.size(); index++ ) { gathered_var.at(index) = (bool)dest.at(index); }
}


size_t dataset::local_index(  const int Itime, const int Idepth, const int Ilat, const int Ilon 
        ) const {
    return Index( Itime, Idepth, Ilat, Ilon, 
                  Ntime, Ndepth, Nlat, Nlon );
}


size_t dataset::global_index(   const int Itime, const int Idepth, const int Ilat, const int Ilon,
                                const std::string merge_kind 
                            ) const{
    const bool merged_time  = ((merge_kind == "Time")  or (merge_kind == "TimeDepth")) ? true : false;
    const bool merged_depth = ((merge_kind == "Depth") or (merge_kind == "TimeDepth")) ? true : false;

    int Itime_global = Itime  + ( merged_time  ? myStarts[0] : 0 ),
        Ntime_global = merged_time ? full_Ntime : Ntime,
        Idepth_global = Idepth + ( merged_depth ? myStarts[1] : 0 ),
        Ndepth_global = merged_depth ? full_Ndepth : Ndepth;
    return Index( Itime_global, Idepth_global, Ilat, Ilon, 
                  Ntime_global, Ndepth_global, Nlat, Nlon );
}

size_t dataset::index_local_to_global(  const size_t index, const std::string merge_kind 
        ) const {
    int Itime, Idepth, Ilat, Ilon;
    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
    return global_index( Itime, Idepth, Ilat, Ilon, merge_kind );
}

size_t dataset::index_global_to_local(  const size_t index, const std::string merge_kind 
        ) const {
    int Itime, Idepth, Ilat, Ilon;

    const bool merged_time  = ((merge_kind == "Time")  or (merge_kind == "TimeDepth")) ? true : false;
    const bool merged_depth = ((merge_kind == "Depth") or (merge_kind == "TimeDepth")) ? true : false;

    int Ntime_global  = merged_time  ? full_Ntime  : Ntime,
        Ndepth_global = merged_depth ? full_Ndepth : Ndepth;

    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime_global, Ndepth_global, Nlat, Nlon );

    int Itime_global  = Itime  - ( merged_time  ? myStarts[0] : 0 ),
        Idepth_global = Idepth - ( merged_depth ? myStarts[1] : 0 );

    return local_index( Itime_global, Idepth_global, Ilat, Ilon );
}

void dataset::index1to4_local(  const size_t index, 
                                int & Itime, int & Idepth, int & Ilat, int & Ilon 
                              ) const {
    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
}

void dataset::index1to4_global( const size_t index, 
                                int & Itime, int & Idepth, int & Ilat, int & Ilon,
                                const std::string merge_kind
        ) const {
    const bool merged_time  = ((merge_kind == "Time")  or (merge_kind == "TimeDepth")) ? true : false;
    const bool merged_depth = ((merge_kind == "Depth") or (merge_kind == "TimeDepth")) ? true : false;

    int Ntime_global  = merged_time  ? full_Ntime  : Ntime,
        Ndepth_global = merged_depth ? full_Ndepth : Ndepth;

    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime_global, Ndepth_global, Nlat, Nlon );
}

void dataset::write_adjacency(
        const std::string filename,
        const MPI_Comm comm
        ) {

    std::vector<std::string> vars_to_write;

    vars_to_write.push_back("adjacency_indices");
    vars_to_write.push_back("adjacency_proj_x");
    vars_to_write.push_back("adjacency_proj_y");
    vars_to_write.push_back("adjacency_distances");
    vars_to_write.push_back("adjacency_ddlon_weights");
    vars_to_write.push_back("adjacency_ddlat_weights");
    vars_to_write.push_back("adjacency_d2dlon2_weights");
    vars_to_write.push_back("adjacency_d2dlat2_weights");

    initialize_adjacency_file( *this, vars_to_write, filename.c_str() );

    std::vector<double> var_to_write;

    const size_t Npts = longitude.size(),
                 Nneighbours = num_neighbours;

    size_t starts[2] = { 0,   0  };
    size_t counts[2] = { Npts, Nneighbours+1 };
    var_to_write.resize( Npts*(Nneighbours+1) );

    // Indices
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = (double) adjacency_indices[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_indices", starts, counts, filename.c_str() );

    // Projected x
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_projected_x[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_proj_x", starts, counts, filename.c_str() );

    // Projected y
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_projected_y[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_proj_y", starts, counts, filename.c_str() );

    // Distances
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_distances[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_distances", starts, counts, filename.c_str() );


    // 1st lon deriv weights
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_ddlon_weights[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_ddlon_weights", 
            starts, counts, filename.c_str() );

    // 1st lat deriv weights
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_ddlat_weights[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_ddlat_weights", 
            starts, counts, filename.c_str() );

    // 2nd lon deriv weights
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_d2dlon2_weights[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_d2dlon2_weights", 
            starts, counts, filename.c_str() );

    // 2nd lat deriv weights
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = adjacency_d2dlat2_weights[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_d2dlat2_weights", 
            starts, counts, filename.c_str() );

}

void dataset::load_adjacency(
        const std::string filename,
        const MPI_Comm comm
        ) {

    std::vector<double> tmp_var;

    const size_t Npts = longitude.size(),
                 Nneighbours = num_neighbours;

    // Indices
    read_var_from_file( tmp_var, 
                        "adjacency_indices", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_indices.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_indices.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_indices[II][JJ] = (size_t) round(tmp_var.at(II*(Nneighbours+1) + JJ));
        } 
    }

    // projected x
    read_var_from_file( tmp_var,
                        "adjacency_proj_x", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_projected_x.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_projected_x.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_projected_x[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // projected y
    read_var_from_file( tmp_var,
                        "adjacency_proj_y", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_projected_y.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_projected_y.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_projected_y[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // distances
    read_var_from_file( tmp_var,
                        "adjacency_distances", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_distances.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_distances.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_distances[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // 1st lon deriv weights
    read_var_from_file( tmp_var,
                        "adjacency_ddlon_weights", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_ddlon_weights.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_ddlon_weights.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_ddlon_weights[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // 1st lat deriv weights
    read_var_from_file( tmp_var,
                        "adjacency_ddlat_weights", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_ddlat_weights.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_ddlat_weights.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_ddlat_weights[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // 2nd lon deriv weights
    read_var_from_file( tmp_var,
                        "adjacency_d2dlon2_weights", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_d2dlon2_weights.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_d2dlon2_weights.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_d2dlon2_weights[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }

    // 2nd lat deriv weights
    read_var_from_file( tmp_var,
                        "adjacency_d2dlat2_weights", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_d2dlat2_weights.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_d2dlat2_weights.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_d2dlat2_weights[II][JJ] = tmp_var.at(II*(Nneighbours+1) + JJ);
        } 
    }


}
