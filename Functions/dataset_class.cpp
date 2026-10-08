#include <math.h>
#include <algorithm>
#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include <omp.h>
#include <Eigen/QR>



// Class constructor
dataset::dataset() {
};

// Copy constructor
dataset::dataset( dataset &set_to_copy ) {
};

// Reset
void dataset::clear() {

    time.clear();
    depth.clear();
    latitude.clear();
    longitude.clear();

    areas.clear();

    variables.clear();

    regions.clear();
    region_names.clear();
    region_areas.clear();
    region_areas_water_only.clear();

    coarse_map_lat.clear();
    coarse_map_lon.clear();
    coarse_map_areas.clear();

    mask.clear();
    reference_mask.clear();
    mask_DEPTH.clear();

    myCounts.clear();
    myStarts.clear();


    adjacency_indices.clear();
    adjacency_indices_adjoint.clear();
    adjacency_projected_x.clear();
    adjacency_projected_y.clear();
    adjacency_distances.clear();
    adjacency_ddlon_weights.clear();
    adjacency_ddlat_weights.clear();
    adjacency_d2dlon2_weights.clear();
    adjacency_d2dlat2_weights.clear();

};

// Copy from pointer
void dataset::copy_from_ptr( const dataset *source ) {
    full_Ntime = source->full_Ntime;
    full_Ndepth = source->full_Ndepth;

    Ntime = source->Ntime;
    Ndepth = source->Ndepth;
    Nlat = source->Nlat;
    Nlon = source->Nlon;

    time = source->time;
    depth = source->depth;

    myCounts = source->myCounts;
    myStarts = source->myStarts;

    longitude = source->longitude;
    latitude = source->latitude;

    mask = source->mask;
    areas = source->areas;

    adjacency_indices = source->adjacency_indices;

    adjacency_ddlon_weights = source->adjacency_ddlon_weights;
    adjacency_ddlat_weights = source->adjacency_ddlat_weights;

    adjacency_d2dlon2_weights = source->adjacency_d2dlon2_weights;
    adjacency_d2dlat2_weights = source->adjacency_d2dlat2_weights;
}

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

void dataset::initialize_index_functions() {
    stride_lat   = Nlon;
    stride_depth = Nlat   * stride_lat;
    stride_time  = Ndepth * stride_depth;
    #if DEBUG>0
    fprintf(stdout, "Index strides are: (time) %zu, (depth) %zu, (latitude) %zu.\n",
            stride_time, stride_depth, stride_lat);
    #endif
}

void dataset::compute_cell_areas() {

    if (not( (Nlat > 0) and (Nlon > 0) )) {
        throw std::runtime_error("Unable to compute cell sizes. Size zero dimensions exist.");
    }


    areas.resize( Nlat * Nlon );
    compute_areas( areas, longitude, latitude );
}


void dataset::compute_sin_cos_lat() {

    sin_lat.resize( latitude.size() );
    cos_lat.resize( latitude.size() );

    sin_lon.resize( longitude.size() );
    cos_lon.resize( longitude.size() );

    const size_t Nlat = latitude.size(),
                 Nlon = longitude.size();
    size_t Ipt;
    #pragma omp parallel default(none) \
    shared( cos_lat, sin_lat, latitude ) \
    private( Ipt ) firstprivate( Nlat )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Nlat; Ipt++ ) {
            cos_lat[Ipt] = cos(latitude[Ipt]);
            sin_lat[Ipt] = sin(latitude[Ipt]);
        }
    }
    #pragma omp parallel default(none) \
    shared( cos_lon, sin_lon, longitude ) \
    private( Ipt ) firstprivate( Nlon )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Nlon; Ipt++ ) {
            cos_lon[Ipt] = cos(longitude[Ipt]);
            sin_lon[Ipt] = sin(longitude[Ipt]);
        }
    }

    sin_cos_lat_pre_comp = true;

}

// This is a wrapper for Functions/distance 
// that takes in indices
double dataset::distance_from_indices(
        const size_t Ipt,
        const size_t Jpt
        ) const {

    if ( Ipt == Jpt ) { return 0.; }

    const bool pre_comp = 
        constants::PRECOMPUTE_SIN_COS_LATITUDE and sin_cos_lat_pre_comp;
    return distance(
                longitude[Ipt],
                latitude[ Ipt],
                longitude[Jpt],
                latitude[ Jpt],
                Cartesian_x_range,
                Cartesian_y_range,
                pre_comp ? cos_lat[Ipt] : 0.,
                pre_comp ? cos_lat[Jpt] : 0.,
                pre_comp ? sin_lat[Ipt] : 0.,
                pre_comp ? sin_lat[Jpt] : 0.,
                pre_comp ? cos_lon[Ipt] : 0.,
                pre_comp ? cos_lon[Jpt] : 0.,
                pre_comp ? sin_lon[Ipt] : 0.,
                pre_comp ? sin_lon[Jpt] : 0.,
                pre_comp
            );

}

// This is a wrapper for Functions/distance 
// that takes in indices
double dataset::distance_between_datasets(
        const size_t Ipt,
        const size_t Jpt,
        const dataset & Jset
        ) const {

    const bool pre_comp = 
        constants::PRECOMPUTE_SIN_COS_LATITUDE 
        and sin_cos_lat_pre_comp
        and Jset.sin_cos_lat_pre_comp;

    return distance(
                longitude[Ipt],
                latitude[ Ipt],
                Jset.longitude[Jpt],
                Jset.latitude[ Jpt],
                Cartesian_x_range,
                Cartesian_y_range,
                pre_comp ?      cos_lat[Ipt] : 0.,
                pre_comp ? Jset.cos_lat[Jpt] : 0.,
                pre_comp ?      sin_lat[Ipt] : 0.,
                pre_comp ? Jset.sin_lat[Jpt] : 0.,
                pre_comp ?      cos_lon[Ipt] : 0.,
                pre_comp ? Jset.cos_lon[Jpt] : 0.,
                pre_comp ?      sin_lon[Ipt] : 0.,
                pre_comp ? Jset.sin_lon[Jpt] : 0.,
                pre_comp
            );

}

// This is a wrapper for Functions/distance 
// that takes in Ilat,Ilon indices
double dataset::distance_from_IlatIlon(
        const size_t Ilat,
        const size_t Ilon,
        const size_t Jlat,
        const size_t Jlon
        ) const {

    if (( Ilat == Jlat ) and (Ilon == Jlon))  { return 0.; }

    const bool pre_comp = 
        constants::PRECOMPUTE_SIN_COS_LATITUDE and sin_cos_lat_pre_comp;
    return distance(
                longitude[Ilon],
                latitude[ Ilat],
                longitude[Jlon],
                latitude[ Jlat],
                Cartesian_x_range,
                Cartesian_y_range,
                pre_comp ? cos_lat[Ilat] : 0.,
                pre_comp ? cos_lat[Jlat] : 0.,
                pre_comp ? sin_lat[Ilat] : 0.,
                pre_comp ? sin_lat[Jlat] : 0.,
                pre_comp ? cos_lon[Ilon] : 0.,
                pre_comp ? cos_lon[Jlon] : 0.,
                pre_comp ? sin_lon[Ilon] : 0.,
                pre_comp ? sin_lon[Jlon] : 0.,
                pre_comp
            );
}

// This is a wrapper for Functions/distance 
// that takes in Ilat,Ilon indices
double dataset::distance_between_datasets_from_IlatIlon(
        const size_t Ilat,
        const size_t Ilon,
        const size_t Jlat,
        const size_t Jlon,
        const dataset & Jset
        ) const {

    const bool pre_comp = 
        constants::PRECOMPUTE_SIN_COS_LATITUDE 
        and sin_cos_lat_pre_comp
        and Jset.sin_cos_lat_pre_comp;
    return distance(
                longitude[Ilon],
                latitude[ Ilat],
                Jset.longitude[Jlon],
                Jset.latitude[ Jlat],
                Cartesian_x_range,
                Cartesian_y_range,
                pre_comp ?      cos_lat[Ilat] : 0.,
                pre_comp ? Jset.cos_lat[Jlat] : 0.,
                pre_comp ?      sin_lat[Ilat] : 0.,
                pre_comp ? Jset.sin_lat[Jlat] : 0.,
                pre_comp ?      cos_lon[Ilon] : 0.,
                pre_comp ? Jset.cos_lon[Jlon] : 0.,
                pre_comp ?      sin_lon[Ilon] : 0.,
                pre_comp ? Jset.sin_lon[Jlon] : 0.,
                pre_comp
            );
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
    
    const int force_split_dim = -1;
    const double land_fill_value = 0;

    // Now read in from the file and store in the variables dictionary
    read_var_from_file( variables.at( var_name ), var_name_in_file, filename, 
                        read_mask ? &mask : NULL, 
                        load_counts ? &myCounts : NULL, 
                        load_counts ? &myStarts : NULL, 
                        Nprocs_in_time, Nprocs_in_depth,
                        do_splits, force_split_dim, land_fill_value, 
                        (Nprocs_in_quadrature == 1) ? MPI_COMM_WORLD : MPI_subcomm_samequadrature );
};

void dataset::check_processor_divisions(    const int Nprocs_in_time_input, 
                                            const int Nprocs_in_depth_input, 
                                            const int Nprocs_in_quad_input, 
                                            const MPI_Comm ) {

    if (constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        if (not( (full_Ntime > 0) and (full_Ndepth > 0) and (Nlon > 0) and (Nlat > 0) )) {
            throw std::runtime_error("Must read in dimensions before checking processor divisions.");
        }
    }

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Apply some cleaning to the processor allotments if necessary. 
    Nprocs_in_time  = ( full_Ntime  == 1 ) ? 1 : 
                      ( full_Ndepth == 1 ) ? wSize / Nprocs_in_quad_input : 
                                             Nprocs_in_time_input;
    Nprocs_in_depth = ( full_Ndepth == 1 ) ? 1 : 
                      ( full_Ntime  == 1 ) ? wSize / Nprocs_in_quad_input : 
                                             Nprocs_in_depth_input;
    Nprocs_in_quadrature = ( full_Ntime * full_Ndepth == 1 ) ? wSize : Nprocs_in_quad_input;

    #if DEBUG >= 0
    if (Nprocs_in_time != Nprocs_in_time_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in time to %'d from %'d\n", 
                Nprocs_in_time, Nprocs_in_time_input); }
    }
    if (Nprocs_in_depth != Nprocs_in_depth_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in depth to %'d from %'d\n", 
                Nprocs_in_depth, Nprocs_in_depth_input); }
    }
    if (Nprocs_in_quadrature != Nprocs_in_quad_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in quadrature to %'d from %'d\n", 
                Nprocs_in_quadrature, Nprocs_in_quad_input); }
    }
    if (wRank == 0) { fprintf(stdout, " Nproc(time, depth) = (%'d, %'d)\n\n", Nprocs_in_time, Nprocs_in_depth); }
    #endif

    if ( Nprocs_in_time * Nprocs_in_depth * Nprocs_in_quadrature != wSize ) {
        throw std::runtime_error("Processor division mismatch. Nprocs_in_time * Nprocs_in_depth != total MPI ranks");
    }

    // Now that processor divisions have been tested, also create the sub-communicator items
    int color, key;

    // communicator for ranks with the same times
    color = (wRank % (Nprocs_in_depth*Nprocs_in_time)) / Nprocs_in_depth;
    key   = (wRank % (Nprocs_in_depth*Nprocs_in_time)) % Nprocs_in_depth;
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_sametimes); 
    #if DEBUG >= 2
    fprintf( stdout, "Processor %d has time color %d and rank %d.\n", wRank, color, key );
    #endif

    // communicator for ranks with the same depths
    color = (wRank % (Nprocs_in_depth*Nprocs_in_time)) % Nprocs_in_depth;
    key   = (wRank % (Nprocs_in_depth*Nprocs_in_time)) / Nprocs_in_depth;
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_samedepths); 
    #if DEBUG >= 2
    fprintf( stdout, "Processor %d has depth color %d and rank %d.\n", wRank, color, key );
    #endif


    // If we have processors for quadrature, we also need to split the communicator
    // there was well. In particular, we want processors to have the same colour
    // if they correspond to the same times and depths
    color = wRank % ( Nprocs_in_time * Nprocs_in_depth );
    key   = wRank / ( Nprocs_in_time * Nprocs_in_depth );
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_sametimedepths); 
    #if DEBUG >= 2
    fprintf( stdout, "Processor %d has time-depth color %d and rank %d.\n", wRank, color, key );
    #endif

    
    color = wRank / ( Nprocs_in_time * Nprocs_in_depth );
    key   = wRank % ( Nprocs_in_time * Nprocs_in_depth );
    MPI_Comm_split( MPI_Comm_Global, color, key, &MPI_subcomm_samequadrature); 
    #if DEBUG >= 2
    if (Nprocs_in_quadrature > 1) {
        fprintf( stdout, "Processor %d has quadrature color %d and rank %d.\n", wRank, color, key );
    }
    #endif

}


void dataset::compute_region_areas() {

    #if DEBUG >= 2
    int wRank=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );

    if (wRank == 0) { fprintf(stdout, "  Computing geographic region areas\n"); fflush(stdout); }
    #endif

    if ( mask.size() == 0 ) { throw std::runtime_error("Must read in mask before computing region areas."); }
    if ( areas.size() == 0 ) { throw std::runtime_error("Must compute cell areas before computing region areas."); }
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


void dataset::prepare_for_coarsened_grids(
        const std::string filename,
        const MPI_Comm comm
        ){


    read_var_from_file( coarse_map_lat, "latitude", filename);
    read_var_from_file( coarse_map_lon, "longitude", filename);

    convert_coordinates( coarse_map_lon, coarse_map_lat );

    // also do the areas of the coarsened map grid
    coarse_map_areas.resize( coarse_map_lat.size() * coarse_map_lon.size() );
    compute_areas( coarse_map_areas, coarse_map_lat, coarse_map_lon );
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

void dataset::gather_mask_across_depth( const std::vector<short int> & var,
                                        std::vector<short int> & gathered_var
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
    /*
    return Index( Itime, Idepth, Ilat, Ilon, 
                  Ntime, Ndepth, Nlat, Nlon );
    */
    size_t index = Itime * stride_time + Idepth * stride_depth + Ilat * stride_lat + Ilon;
    return index;
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
    //Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
    size_t tmp_index = index;

    Itime      = tmp_index / stride_time;
    tmp_index -= Itime * stride_time;

    Idepth     = tmp_index / stride_depth;
    tmp_index -= Idepth * stride_depth;

    Ilat = tmp_index / stride_lat;
    Ilon = tmp_index - Ilat * stride_lat;
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
    //vars_to_write.push_back("adjacency_ddlon_weights");
    //vars_to_write.push_back("adjacency_ddlat_weights");
    //vars_to_write.push_back("adjacency_d2dlon2_weights");
    //vars_to_write.push_back("adjacency_d2dlat2_weights");

    initialize_adjacency_file( *this, vars_to_write, filename.c_str() );

    std::vector<double> var_to_write;

    const size_t Npts = longitude.size(),
                 Nneighbours = num_neighbours;

    // Cell areas
    size_t starts_1D[1] = { 0    };
    size_t counts_1D[1] = { Npts };
    write_field_to_output( areas, "dA", starts_1D, counts_1D, filename.c_str() );


    size_t starts[2] = { 0,   0  };
    size_t counts[2] = { Npts, Nneighbours+1 };
    var_to_write.resize( Npts*(Nneighbours+1) );

    // Indices
    for ( size_t II = 0; II < Npts; II++ ) { for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            var_to_write.at( (Nneighbours+1)*II + JJ ) = (double) adjacency_indices[II][JJ];
    } }
    write_field_to_output( var_to_write, "adjacency_indices", starts, counts, filename.c_str() );

    /*
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
    */

    /*
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
    */

}

void dataset::load_adjacency_indices(
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
}

void dataset::load_adjacency(
        const std::string filename,
        const MPI_Comm comm
        ) {

    std::vector<double> tmp_var;

    /*
    // Get dimensions
    int ncid = nc_open( filename.c_str(), NC_NOWRITE );
    int neighbour_dimid;
    nc_inq_dimid( ncid, "neighbour", &neighbour_dimid );
    nc_inq_dimlen( ncid, dimid, &num_neighbours );
    */

    const size_t Npts = longitude.size(),
                 Nneighbours = num_neighbours;

    // Indices
    read_var_from_file( tmp_var, 
                        "adjacency_indices", 
                        filename, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, false );
    adjacency_indices.resize(Npts);
    adjacency_ddlon_weights.resize(Npts);
    adjacency_ddlat_weights.resize(Npts);
    for ( size_t II = 0; II < Npts; II++ ) { 
        adjacency_ddlon_weights.at(II).resize(Nneighbours+1);
        adjacency_ddlat_weights.at(II).resize(Nneighbours+1);
        adjacency_indices.at(II).resize(Nneighbours+1);
        for ( size_t JJ = 0; JJ < Nneighbours+1; JJ++ ) {
            adjacency_indices[II][JJ] = (size_t) round(tmp_var.at(II*(Nneighbours+1) + JJ));
        } 
    }
}

/*
    size_t Ipt;
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate(Nneighbours, Npts) 
    {
        const size_t size_neighbourhood = Nneighbours+1;
        const size_t Nderiv = size_neighbourhood-1;
        //const int Nconds = 10; // cubic
        const int Nconds = 3; // linear
        const bool keep_latlon = false;

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            double centre_lat = latitude[Ipt];

            long double mean_dist = 0,
                        min_lon =  M_PI,
                        max_lon = -M_PI,
                        min_lat =  M_PI,
                        max_lat = -M_PI,
                        min_x   =  10,
                        max_x   = -10,
                        min_y   =  10,
                        max_y   = -10;
            for ( size_t In = 0; In < Nderiv; In++ ) {

                long double local_lon = longitude[ adjacency_indices[Ipt][In] ] - longitude[Ipt];
                if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                long double local_lat = latitude[ adjacency_indices[Ipt][In] ];

                min_lon = std::fmin( local_lon, min_lon );
                max_lon = std::fmax( local_lon, max_lon );

                min_lat = std::fmin( local_lat, min_lat );
                max_lat = std::fmax( local_lat, max_lat );

                long double x,y;
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    x = local_lon;
                    y = local_lat - centre_lat;
                } else {
                    // Convert to Cartesian
                    // Rotate to have x = 0 at centre
                    // Rotate to have z = 1 at centre
                    // Now a local tangent plane centred at the north pole
                    x =   cos( local_lat ) * sin( local_lon );
                    y =   cos( centre_lat ) * sin( local_lat ) 
                        - sin( centre_lat ) * cos( local_lat ) * cos( local_lon );
                }

                mean_dist += sqrt( x*x + y*y ) / Nderiv;

                min_x = std::fmin( x, min_x );
                min_y = std::fmin( y, min_y );

                max_x = std::fmax( x, max_x );
                max_y = std::fmax( y, max_y );
            }
            const long double x_scale = max_x - min_x,
                              y_scale = max_y - min_y;

            // Regression
            std::map< std::string, long double > Sij;
            std::vector< std::string > S_comps = 
                (Nconds == 10) ? std::vector<std::string> {"", "x", "y", "xy", "xx", "yy", "xxx", "xxy", "xyy", "yyy"}
              : (Nconds ==  6) ? std::vector<std::string> {"", "x", "y", "xy", "xx", "yy"}
              :                  std::vector<std::string> {"", "x", "y"};
            for ( size_t i = 0; i < S_comps.size(); i++ ) {
                for ( size_t j = 0; j <= i; j++ ) {
                    std::string name = S_comps[i] + S_comps[j];
                    Sij.insert( std::pair< std::string, long double >(name, 0.) );
                }
            }
            for ( size_t In = 0; In < Nderiv; In++ ) {

                long double local_lon = longitude[ adjacency_indices[Ipt][In] ] - longitude[Ipt];
                if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                long double local_lat = latitude[ adjacency_indices[Ipt][In] ];

                long double x,y;
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    x = local_lon;
                    y = local_lat - centre_lat;
                } else {
                    // Convert to Cartesian
                    // Rotate to have x = 0 at centre
                    // Rotate to have z = 1 at centre
                    // Now a local tangent plane centred at the north pole
                    x =   cos( local_lat ) * sin( local_lon );
                    y =   cos( centre_lat ) * sin( local_lat ) 
                        - sin( centre_lat ) * cos( local_lat ) * cos( local_lon );
                }

                x /= mean_dist;
                y /= mean_dist;
                
                long double dA = areas[adjacency_indices[Ipt][In]];

                for (auto const& sij : Sij ) {
                    std::string name = sij.first;
                    long double value = 1.;
                    for ( size_t i = 0; i < name.size(); i++ ) {
                        if ( name[i] == 'x' ) { value *= x; }
                        else {                  value *= y; }
                    }
                    Sij[name] += value * dA;
                }
            }

            Eigen::Matrix<long double, Nconds, Nconds> A;
            for ( size_t i = 0; i < S_comps.size(); i++ ) {
                for ( size_t j = 0; j <= i; j++ ) {
                    std::string name = S_comps[i] + S_comps[j];
                    A(i,j) = Sij[name];
                    A(j,i) = Sij[name];
                }
            }


            for ( size_t In = 0; In < Nderiv; In++ ) {

                long double local_lon = longitude[ adjacency_indices[Ipt][In] ] - longitude[Ipt];
                if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                long double local_lat = latitude[ adjacency_indices[Ipt][In] ];

                long double x,y;
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    x = local_lon;
                    y = local_lat - centre_lat;
                } else {
                    // Convert to Cartesian
                    // Rotate to have x = 0 at centre
                    // Rotate to have z = 1 at centre
                    // Now a local tangent plane centred at the north pole
                    x =   cos( local_lat ) * sin( local_lon );
                    y =   cos( centre_lat ) * sin( local_lat ) 
                        - sin( centre_lat ) * cos( local_lat ) * cos( local_lon );
                }

                x /= mean_dist;
                y /= mean_dist;
                
                long double dA = areas[adjacency_indices[Ipt][In]];

                Eigen::Vector<long double, Nconds> RHS, sol;
                RHS[0] = 1.  * dA;
                RHS[1] = x   * dA;
                RHS[2] = y   * dA;
                if ( Nconds > 3 ) {
                    RHS[3] = x*y * dA;
                    RHS[4] = x*x * dA;
                    RHS[5] = y*y * dA;
                }
                if ( Nconds > 6 ) {
                    RHS[6] = x*x*x * dA;
                    RHS[7] = x*x*y * dA;
                    RHS[8] = x*y*y * dA;
                    RHS[9] = y*y*y * dA;
                }
                sol = A.completeOrthogonalDecomposition().solve(RHS);

                const bool DDLON_IS_ZERO = false;//std::fabs(max_lon - min_lon) < 1e-4;
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    adjacency_ddlon_weights[Ipt][In] = DDLON_IS_ZERO ? 0. : sol[1]                   / mean_dist;
                } else {
                    adjacency_ddlon_weights[Ipt][In] = DDLON_IS_ZERO ? 0. : sol[1] * cos(centre_lat) / mean_dist;
                }

                const bool DDLAT_IS_ZERO = false;//std::fabs(max_lat - min_lat) < 1e-4;
                adjacency_ddlat_weights[Ipt][In] = DDLAT_IS_ZERO ? 0. : sol[2] / mean_dist;
            }
*/

            /* Polynomial fit
            // use a local Cartesian projection to get derivative weights near the poles
            // we explicitely exclude the centre point from the weighting for stability
            Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A;
            A.resize(Nconds,Nderiv);
            for ( size_t In = 0; In < Nderiv; In++ ) {

                long double local_lon = longitude[ adjacency_indices[Ipt][In] ] - longitude[Ipt];
                if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                long double local_lat = latitude[ adjacency_indices[Ipt][In] ];

                long double x,y;
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    x = local_lon;
                    y = local_lat - centre_lat;
                } else {
                    // Convert to Cartesian
                    // Rotate to have x = 0 at centre
                    // Rotate to have z = 1 at centre
                    // Now a local tangent plane centred at the north pole
                    x =   cos( local_lat ) * sin( local_lon );
                    y =   cos( centre_lat ) * sin( local_lat ) 
                        - sin( centre_lat ) * cos( local_lat ) * cos( local_lon );
                }

                x /= mean_dist;
                y /= mean_dist;

                // We construct the local polynomial approximation
                A(0,In) = 1.;
                A(1,In) = x;
                A(2,In) = y;
                if ( Nconds == 3 ) { continue; }

                A(3,In) = x * y; if ( Nconds == 4 ) { continue; }
                A(4,In) = 0.5 * pow( x, 2. );
                A(5,In) = 0.5 * pow( y, 2. );
                if ( Nconds == 6 ) { continue; }

                A(6,In) = (1./6.) * pow( x, 3. );
                A(7,In) = (1./2.) * pow( x, 2. ) * y;
                A(8,In) = (1./2.) * pow( y, 2. ) * x;
                A(9,In) = (1./6.) * pow( y, 3. );
                if ( Nconds == 10 ) { continue; }

                A(10,In) = (1./24.) * pow( x, 4. ) * pow( y, 0. );
                A(11,In) = (1./6. ) * pow( x, 3. ) * pow( y, 1. );
                A(12,In) = (1./4. ) * pow( x, 2. ) * pow( y, 2. );
                A(13,In) = (1./6. ) * pow( x, 1. ) * pow( y, 3. );
                A(14,In) = (1./24.) * pow( x, 0. ) * pow( y, 4. );
            }

            Eigen::Vector<long double, Eigen::Dynamic> RHS, ddx, ddy;
            RHS.resize(Nconds);
            for ( int II = 0; II < Nconds; II++ ) { RHS[II] = 0.; }

            const bool DDLON_IS_ZERO = false; //std::fabs(max_lon - min_lon) < 1e-3 * M_PI / 180.;
            if ( DDLON_IS_ZERO ) {
                for ( size_t In = 0; In < Nderiv; In++ ) {
                    adjacency_ddlon_weights[Ipt][In] = 0.;
                }
            } else {
                RHS[1] = 1.;
                ddx = A.completeOrthogonalDecomposition().solve(RHS);
                if ( keep_latlon and (std::fabs(centre_lat) < 88. * M_PI / 180.) ) {
                    for ( size_t In = 0; In < Nderiv; In++ ) {
                        adjacency_ddlon_weights[Ipt][In] = ddx(In) / mean_dist;
                    }
                } else {
                    for ( size_t In = 0; In < Nderiv; In++ ) {
                        adjacency_ddlon_weights[Ipt][In] = cos(centre_lat) * ddx(In) / mean_dist;
                    }
                }
            }

            RHS[1] = 0.;
            RHS[2] = 1.;
            ddy = A.completeOrthogonalDecomposition().solve(RHS);

            const bool DDLAT_IS_ZERO = false; //std::fabs(max_lat - min_lat) < 1e-3 * M_PI / 180.;
            for ( size_t In = 0; In < Nderiv; In++ ) {
                adjacency_ddlat_weights[Ipt][In] = DDLAT_IS_ZERO ? 0. : ddy(In) / mean_dist;
            }
            */


//        }
//    }
        
    //
    //// Setting up a least-squares adjoint system for div and grad
    //
    /*

    grad_weights_lon.resize(Npts);
    grad_weights_lat.resize(Npts);

    perp_grad_weights_lon.resize(Npts);
    perp_grad_weights_lat.resize(Npts);

    for ( size_t II = 0; II < Npts; II++ ) {
        grad_weights_lon.at(II).resize(Nneighbours);
        grad_weights_lat.at(II).resize(Nneighbours);
        perp_grad_weights_lon.at(II).resize(Nneighbours);
        perp_grad_weights_lat.at(II).resize(Nneighbours);
    }
    //size_t Ipt;
    size_t tikhov_regularized = 0;
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate(Nneighbours, Npts) \
    reduction( +:tikhov_regularized )
    {
        const size_t size_neighbourhood = Nneighbours;

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            const double lat_i = latitude[ Ipt];
            const double lon_i = longitude[Ipt];

            // outward normal vector
            Vec3 x_i = { cos(lat_i) * cos(lon_i), 
                         cos(lat_i) * sin(lon_i), 
                         sin(lat_i) 
            };

            // get tangent plane basis vectrs
            Vec3 e_lon = { -sin(lon_i), 
                            cos(lon_i), 
                            0.0 
            };
            Vec3 e_lat = cross( x_i, e_lon );
            //Vec3 e_lat = { - cos(lon_i) * sin(lat_i),
            //               - sin(lon_i) * sin(lat_i),
            //                 cos(lat_i)
            //};

            std::vector<double> A(4, 0.);

            for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                const size_t Jpt = adjacency_indices[Ipt][In];

                const double lat_j = latitude[ Jpt];
                const double lon_j = longitude[Jpt];

                Vec3 x_j = { cos(lat_j) * cos(lon_j), 
                             cos(lat_j) * sin(lon_j), 
                             sin(lat_j) 
                };

                // get displacement in local tangent plane
                Vec3 r_ij_3d = x_j - dot( x_j, x_i ) * x_i;

                // convert displacement to a 2-vector in local tangent plane
                Vec2 r_ij_2d = { dot(r_ij_3d, e_lon ), dot(r_ij_3d, e_lat) };

                // Accumulate A = sum( r_ij r_ij^T )
                A[0] += r_ij_2d[0] * r_ij_2d[0];
                A[1] += r_ij_2d[1] * r_ij_2d[0];
                A[2] += r_ij_2d[0] * r_ij_2d[1];
                A[3] += r_ij_2d[1] * r_ij_2d[1];
            }

            //
            //// Invert A (Tikhov regularize if necessary)
            //
            double det = A[0]*A[3] - A[1]*A[2];
            const double trace = A[0] + A[3];

            const double eigen_1 = 0.5 * ( trace + sqrt( trace*trace - 4 * det ) );
            const double eigen_2 = 0.5 * ( trace - sqrt( trace*trace - 4 * det ) );
            const double kappa = (eigen_2 == 0) ? 1e100 : eigen_1 / eigen_2;

            if ( kappa > 1e2 ) {
                tikhov_regularized += 1;
                const double alpha = 1e-6;
                const double eps_abs = 1e-12;

                const double eps_tikhov = std::fmax( alpha * trace, eps_abs );

                // change from A to A + eps*I
                A[0] += eps_tikhov;
                A[3] += eps_tikhov;
                det = A[0]*A[3] - A[1]*A[2];
            }

            const double invdet = 1.0 / det;
            std::vector<double> invA(4, 0.);
            invA[0] =  A[3] * invdet;
            invA[1] = -A[1] * invdet;
            invA[2] = -A[2] * invdet;
            invA[3] =  A[0] * invdet;

            for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                const size_t Jpt = adjacency_indices[Ipt][In];

                const double lat_j = latitude[ Jpt];
                const double lon_j = longitude[Jpt];

                Vec3 x_j = { cos(lat_j) * cos(lon_j), 
                             cos(lat_j) * sin(lon_j), 
                             sin(lat_j) 
                };

                // get displacement in local tangent plane
                Vec3 r_ij_3d = x_j - dot( x_j, x_i ) * x_i;

                // convert displacement to a 2-vector in local tangent plane
                Vec2 r_ij_2d = { dot(r_ij_3d, e_lon ), dot(r_ij_3d, e_lat) };

                // Compute and store w_ij as gradient weights
                grad_weights_lon[Ipt][In] = ( invA[0] * r_ij_2d[0] + invA[1] * r_ij_2d[1] ) / constants::R_earth;
                grad_weights_lat[Ipt][In] = ( invA[2] * r_ij_2d[0] + invA[3] * r_ij_2d[1] ) / constants::R_earth;

                perp_grad_weights_lon[Ipt][In] = - grad_weights_lat[Ipt][In];
                perp_grad_weights_lat[Ipt][In] =   grad_weights_lon[Ipt][In];

                // get the perpendicular weights, as a 3-vector
                //const Vec3 w_ij_perp_3d =   grad_weights_lon[Ipt][In] * cross( x_i, e_lon )
                //                          + grad_weights_lat[Ipt][In] * cross( x_i, e_lat );

                // convert to a 2-vector in local tangent plane and store
                //perp_grad_weights_lon[Ipt][In] = dot(w_ij_perp_3d, e_lat );
                //perp_grad_weights_lat[Ipt][In] = dot(w_ij_perp_3d, e_lon );
            }
        }
    }
    fprintf(stdout, "  !! %'zu points were tikhov regularized.\n", tikhov_regularized);

    // THIS ONLY POSSIBLE IF GUARANTEED MUTUAL NEIGHBOURS
    // Symmetrize the gradient weights
    #pragma omp parallel default(none) private( Ipt ) shared(stdout) \
    firstprivate(Nneighbours, Npts)
    {
        const size_t size_neighbourhood = Nneighbours;

        #pragma omp for collapse(1) schedule(static,4)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                const size_t Jpt = adjacency_indices[Ipt][In];

                if ( Jpt >= Ipt ) { continue; }

                size_t Jn;
                for ( Jn = 0; Jn < size_neighbourhood; Jn++ ) {
                    if ( Ipt == adjacency_indices[Jpt][Jn] ) { break; }
                }
                if ( Jn >= size_neighbourhood ) {
                    fprintf( stdout, "!! Pt %'zu N %'zu (p %'zu) are not mutual neighbours.",
                            Ipt, In, Jpt );
                    throw std::runtime_error("Unable to symmetric gradient weights. Adjacency is not mutual.");
                }

                double w_ij_lon = 0.5 * ( grad_weights_lon[Ipt][In] + grad_weights_lon[Jpt][Jn] );
                double w_ij_lat = 0.5 * ( grad_weights_lat[Ipt][In] + grad_weights_lat[Jpt][Jn] );

                grad_weights_lon[Ipt][In] = w_ij_lon;
                grad_weights_lon[Jpt][Jn] = w_ij_lon;

                grad_weights_lat[Ipt][In] = w_ij_lat;
                grad_weights_lat[Jpt][Jn] = w_ij_lat;

                perp_grad_weights_lon[Ipt][In] = - w_ij_lat;
                perp_grad_weights_lon[Jpt][Jn] = - w_ij_lat;

                perp_grad_weights_lat[Ipt][In] =   w_ij_lon;
                perp_grad_weights_lat[Jpt][Jn] =   w_ij_lon;
            }
        }
    }
    */


void dataset::build_differentiation_weights(
        const MPI_Comm comm
        ) {

    const size_t Npts = longitude.size(),
                 Nneighbours = num_neighbours;

    //
    //// A DEC-type scheme for differentiation weights
    //
    grad_weights_lon.resize(Npts);
    grad_weights_lat.resize(Npts);

    perp_grad_weights_lon.resize(Npts);
    perp_grad_weights_lat.resize(Npts);

    for ( size_t II = 0; II < Npts; II++ ) {
        grad_weights_lon.at(II).resize(Nneighbours);
        grad_weights_lat.at(II).resize(Nneighbours);
        perp_grad_weights_lon.at(II).resize(Nneighbours);
        perp_grad_weights_lat.at(II).resize(Nneighbours);
    }
    size_t Ipt;
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate(Nneighbours, Npts)
    {
        const size_t size_neighbourhood = Nneighbours;

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            const double lat_i = latitude[ Ipt];
            const double lon_i = longitude[Ipt];

            // outward normal vector
            Vec3 x_i = { cos(lat_i) * cos(lon_i), 
                         cos(lat_i) * sin(lon_i), 
                         sin(lat_i) 
            };

            // eastward basis vector
            Vec3 e_lon = { -sin(lon_i), 
                            cos(lon_i), 
                            0.0 
            };

            // northward basis vector
            //Vec3 e_lat = cross( x_i, e_lon );
            Vec3 e_lat = { - cos(lon_i) * sin(lat_i),
                           - sin(lon_i) * sin(lat_i),
                             cos(lat_i)
            };

            const double A_i = areas[Ipt];

            for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                const size_t Jpt = adjacency_indices[Ipt][In];

                // Distance between points
                const double ell_ij = distance_from_indices( Ipt, Jpt );

                const double lat_j = latitude[ Jpt];
                const double lon_j = longitude[Jpt];

                Vec3 x_j = { cos(lat_j) * cos(lon_j), 
                             cos(lat_j) * sin(lon_j), 
                             sin(lat_j) 
                };

                // get displacement in local tangent plane
                Vec3 v_ij = x_j - dot( x_j, x_i ) * x_i;

                Vec3 e_ij_hat = ( 1. / sqrt( dot( v_ij, v_ij ) ) ) * v_ij;
                Vec3 t_ij_hat = cross( x_i, e_ij_hat );

                // Compute and store w_ij as gradient weights
                grad_weights_lon[Ipt][In] = (ell_ij / A_i) * dot( e_ij_hat, e_lon ); // / constants::R_earth;
                grad_weights_lat[Ipt][In] = (ell_ij / A_i) * dot( e_ij_hat, e_lat ); // / constants::R_earth;

                perp_grad_weights_lon[Ipt][In] = - grad_weights_lat[Ipt][In];
                perp_grad_weights_lat[Ipt][In] =   grad_weights_lon[Ipt][In];
            }
        }
    }

}

