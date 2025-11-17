#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <utility> // std::swap
#include <vector>
#include <deque>
#include <omp.h>
#include <math.h>
#include <numeric> // std::accumulate
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/QR>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_long_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void SouthPolar_Bivariate_Regression(
        double & alpha_const,
        double & beta_x,
        double & beta_y,
        double & gamma_xx,
        double & gamma_xy,
        double & gamma_yy,
        double & delta_xxx,
        double & delta_xxy,
        double & delta_xyy,
        double & delta_yyy,
        double & x0,
        double & y0,
        const std::vector<double> & lons,
        const std::vector<double> & lats,
        const std::vector<double> & vals,
        const bool SouthPolarProjection = true
        ){

    const size_t Nvals = vals.size();
    const size_t Nconds = 6; //3;

    std::map< std::string, long double > Sij;
    std::vector< std::string > S_comps = 
        (Nconds == 10) ? std::vector<std::string> {"", "x", "y", "xx", "xy", "yy", "xxx", "xxy", "xyy", "yyy"}
      : (Nconds ==  6) ? std::vector<std::string> {"", "x", "y", "xx", "xy", "yy"}
      :                  std::vector<std::string> {"", "x", "y"};
    for ( size_t i = 0; i < S_comps.size(); i++ ) {
        for ( size_t j = 0; j <= i; j++ ) {
            std::string name = S_comps[i] + S_comps[j];
            Sij.insert( std::pair< std::string, long double >(name, 0.) );
        }
    }


    // Compute the Sx, Sxx, Sxy, etc. values
    //  and build the RHS
    Eigen::VectorXd RHS, sol;
    RHS.resize(Nconds);
    for ( size_t In = 0; In < Nconds; In++ ) { RHS[In] = 0.; }


    for ( size_t In = 0; In < Nvals; In++ ) {

        long double local_lon = lons[In];
        long double local_lat = lats[In];

        long double x,y;
        x = sin( local_lon ) * cos( local_lat );
        y = cos( local_lon ) * cos( local_lat );

        long double dA = 1.;

        for (auto const& sij : Sij ) {
            std::string name = sij.first;
            long double value = 1.;
            for ( size_t i = 0; i < name.size(); i++ ) {
                if ( name[i] == 'x' ) { value *= x; }
                else {                  value *= y; }
            }
            Sij[name] += value * dA;
        }

        RHS[0] += dA * vals[In];
        RHS[1] += dA * vals[In] * x;
        RHS[2] += dA * vals[In] * y;
        if ( Nconds > 3 ) {
            RHS[3] += dA * vals[In] * x * x;
            RHS[4] += dA * vals[In] * x * y;
            RHS[5] += dA * vals[In] * y * y;
        }
        if ( Nconds > 6 ) {
            RHS[6] += dA * vals[In] * x * x * x;
            RHS[7] += dA * vals[In] * x * x * y;
            RHS[8] += dA * vals[In] * x * y * y;
            RHS[9] += dA * vals[In] * y * y * y;
        }
    }

    // Assemble info a matrix
    Eigen::MatrixXd A;
    A.resize(Nconds, Nconds);
    for ( size_t i = 0; i < S_comps.size(); i++ ) {
        for ( size_t j = 0; j <= i; j++ ) {
            std::string name = S_comps[i] + S_comps[j];
            A(i,j) = Sij[name];
            A(j,i) = Sij[name];
        }
    }

    // Get the coefficients
    sol = A.completeOrthogonalDecomposition().solve(RHS);

    alpha_const = sol[0];

    beta_x = sol[1];
    beta_y = sol[2];

    gamma_xx = (Nconds > 3) ? sol[3] : 0.;
    gamma_xy = (Nconds > 3) ? sol[4] : 0.;
    gamma_yy = (Nconds > 3) ? sol[5] : 0.;

    delta_xxx = (Nconds > 6) ? sol[6] : 0.;
    delta_xxy = (Nconds > 6) ? sol[7] : 0.;
    delta_xyy = (Nconds > 6) ? sol[8] : 0.;
    delta_yyy = (Nconds > 6) ? sol[9] : 0.;

    x0 = 0.; //centre_lon;
    y0 = 0.; //centre_lat;

}

void SouthPolar_ApplyHaloRegression(
        const std::vector<size_t> & halo_indices,
        const std::vector<size_t> & core_indices,
        const dataset & grid,
        std::vector<double> & data,
        const bool SouthPolarProjection = true
        ){

    const size_t NHalo = halo_indices.size();

    std::vector<double> lons(NHalo, 0.), lats(NHalo, 0.), vals(NHalo, 0.);
    for ( size_t IHalo = 0; IHalo < NHalo; IHalo++ ) {
        size_t HaloPt = halo_indices[IHalo];

        lons[IHalo] = grid.longitude[HaloPt];
        lats[IHalo] = grid.latitude[HaloPt];
        vals[IHalo] = data[HaloPt];
    }

    double alpha_const, beta_x, beta_y, gamma_xx, gamma_xy, gamma_yy, x0, y0,
           delta_xxx, delta_xxy, delta_xyy, delta_yyy;

    SouthPolar_Bivariate_Regression(
            alpha_const, beta_x, beta_y,
            gamma_xx, gamma_xy, gamma_yy,
            delta_xxx, delta_xxy, delta_xyy, delta_yyy,
            x0, y0,
            lons, lats, vals
            );

    //fprintf( stdout, "Regression params: %.2e (%.2e %.2e) (%.2e %.2e %.2e) (%.2e %.2e %.2e %.2e)\n", 
    //        alpha_const, beta_x, beta_y, gamma_xx, gamma_xy, gamma_yy, 
    //        delta_xxx, delta_xxy, delta_xyy, delta_yyy );


    const size_t NCore = core_indices.size();
    for ( size_t ICore = 0; ICore < NCore; ICore++ ) {
        size_t CorePt = core_indices[ICore];

        double x, y;
        double local_lon = grid.longitude[CorePt];
        double local_lat = grid.latitude[ CorePt];

        x = sin( local_lon ) * cos( local_lat );
        y = cos( local_lon ) * cos( local_lat );

        double f = alpha_const 
                    + beta_x * x
                    + beta_y * y 
                    + gamma_xx * x * x
                    + gamma_xy * x * y
                    + gamma_yy * y * y
                    + delta_xxx * x * x * x
                    + delta_xxy * x * x * y
                    + delta_xyy * x * y * y
                    + delta_yyy * y * y * y;
        data[CorePt] = f;
    }

}

void Scalar_Land_Filler(
        const std::string output_fname,
        dataset & source_data,
        const std::string scalar_name,
        const double rel_tol,
        const int num_refinements,
        const double CFL,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    const size_t Npts = source_data.mask.size();

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    dataset solution_grid;
    dataset *coarsened_grid = new dataset();

    // Beginning the v-cycle [or, more acurrately, just a / cycle?]
    // Start at the coarsest resolution, solve, refine, repeat
    int Nlat_coarse = 0, Nlon_coarse = 0;
    size_t Npts_coarse = 0, Ipt;
    for ( int refine_level = num_refinements; refine_level >= 0; refine_level-- ) {

        #if DEBUG >= 0
        fprintf(stdout, "\n\nBeginning refinement level %d.\n", refine_level);
        #endif
        std::vector<short int> saved_mask;

        // If refinement is turned on, do it
        if ( refine_level > 0 ) {
            Nlat_coarse = floor(sqrt( 0.5 * (Npts / pow(4., 2*refine_level)) ));
            Nlon_coarse = 2 * Nlat_coarse;
            Npts_coarse = ((size_t) Nlat_coarse) * ((size_t) Nlon_coarse);

            #if DEBUG >= 2
            fprintf(stdout, "Refinement level has Nlat x Nlon (%zu) of %d x %d.\n", 
                    Npts_coarse, Nlat_coarse, Nlon_coarse);
            #endif


            // Downsample the velocities onto the coarse grid
            std::vector<std::string> vars_to_map = { scalar_name };

            initialize_coarsened_grid( *coarsened_grid, source_data, Nlat_coarse, Nlon_coarse );

            //Npts_coarse = coarsened_grid->longitude.size(); // The coarsened grid might have a different size
            for ( size_t Ivar = 0; Ivar < vars_to_map.size(); Ivar++ ) {
                coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                            vars_to_map[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
            }

            #if DEBUG >= 2
            fprintf(stdout, "Mapping velocities onto the coarse grid.\n");
            #endif
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( source_data, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

            coarsened_grid->adjacency_indices_adjoint.resize(0);
        } else {
            // Otherwise we're on the original grid, so just use the source
            coarsened_grid = &source_data;
            Npts_coarse = Npts;

            Nlat_coarse = source_data.Nlat;
            Nlon_coarse = source_data.Nlon;
        }

        // Add land-filled scalar variable
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "land_filled", std::vector<double>(Npts_coarse, 0.) ) );

        // Now, up-sample the previous solutions onto the new grid
        if ( refine_level < num_refinements ) {

            std::vector<std::string> vars_to_map = { "land_filled" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( solution_grid, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Upsampling"); }
        }

        // Re-initialize the solution grid on the new computational grid
        // To avoid having to re-build, just copy from the coarsened grid
        solution_grid.clear();
        solution_grid.copy_from_ptr( coarsened_grid );
        std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 

        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "land_filled", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.Nlat = coarsened_grid->Nlat;
        solution_grid.Nlon = coarsened_grid->Nlon;

        #pragma omp parallel default(none) private(Ipt) \
        shared( solution_grid, coarsened_grid ) \
        firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(guided)
            for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                if ( coarsened_grid->mask[Ipt] ) {
                    solution_grid.variables["land_filled"][Ipt] = coarsened_grid->variables["scalar"][Ipt];
                } else {
                    solution_grid.variables["land_filled"][Ipt] = coarsened_grid->variables["land_filled"][Ipt];
                }
            }
        }

        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.iteration_max = Npts_coarse * 10;
        Helmholtz_data.iterations_per_cycle = std::fmin( 10000, Npts_coarse / 100 );

        // Seed with zonal average on coarsest grid
        if ( refine_level == num_refinements ) {
            int Ilat;
            std::vector<int> zonal_counts(Nlat_coarse, 0);
            std::vector<double> zonal_means(Nlat_coarse, 0.);
            #pragma omp parallel default(none) private( Ilat ) \
            shared( solution_grid, coarsened_grid, zonal_counts, zonal_means ) \
            firstprivate( Nlat_coarse, Nlon_coarse )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                    long double zonal_mean = 0.;
                    size_t zonal_count = 0;
                    for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                        size_t index = Ilat*Nlon_coarse + Ilon;
                        if (coarsened_grid->mask[index]) {
                            zonal_mean += solution_grid.variables["land_filled"][index];
                            zonal_count++;
                        }
                    }
                    zonal_counts[Ilat] = zonal_count;
                    zonal_mean = (zonal_count == 0) ? 0. : zonal_mean / zonal_count;
                    zonal_means[Ilat] = zonal_mean;
                    for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                        size_t index = Ilat*Nlon_coarse + Ilon;
                        if (not(coarsened_grid->mask[index])) {
                            solution_grid.variables["land_filled"][index] = zonal_mean;
                        }
                    }
                }
            }

            // Account for south polar region, which has no zonal average to use
            size_t good_index = std::upper_bound (zonal_counts.begin(), zonal_counts.end(), 0) - zonal_counts.begin();
            for ( Ilat = 0; Ilat < good_index; Ilat++ ) {
                for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                    size_t index = Ilat*Nlon_coarse + Ilon;
                    if (not(coarsened_grid->mask[index])) {
                        solution_grid.variables["land_filled"][index] = zonal_means[good_index];
                    }
                }
            }

            coarsened_grid->variables.at("land_filled") = solution_grid.variables.at("land_filled");
        }

        std::vector<double> min_delta(Npts_coarse, 0.);
        {
            int Ilat, Ilon;
            #pragma omp parallel default(none) private( Ilat, Ilon ) shared( min_delta, coarsened_grid ) \
            firstprivate( Nlat_coarse, Nlon_coarse )
            {
                #pragma omp for collapse(2) schedule(static)
                for ( Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                    for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                        size_t index = Ilat*Nlon_coarse + Ilon;
                        double cos_lat = cos( coarsened_grid->latitude[Ilat] );
                        min_delta[index] = constants::R_earth * ( 2 * M_PI / Nlon_coarse ) * cos_lat;
                    }
                }
            }
        }

        bool keep_solving = true;
        unsigned long int iter_cycle = 0;

        std::vector<double> Laplace_of_Scalar(Npts_coarse, 0.);
        scalar_laplacian( Laplace_of_Scalar, solution_grid.variables["land_filled"],
                *coarsened_grid, solution_grid.mask );

        // Write seed
        {
            std::ostringstream name_stream;
            name_stream << "seed_R" << (refine_level) << "_" << (scalar_name) <<  ".nc";
            const std::string output_filename = name_stream.str();
            std::vector<std::string> vars_to_write = { scalar_name, "Laplacian" };
            initialize_output_file( *coarsened_grid, vars_to_write, output_filename, -1 );

            const int ndims = 4;
            size_t starts[ndims] = { 0, 0, 0, 0 };
            size_t counts[ndims] = { 1, 1, size_t(coarsened_grid->Nlat), size_t(coarsened_grid->Nlon) };
            write_field_to_output( solution_grid.variables.at("land_filled"), 
                    scalar_name, starts, counts, output_filename, NULL );
            write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
        }

        #if DEBUG >= 0
        fprintf( stdout, "Beginning solver.\n");
        #endif
        while ( keep_solving ) {


            scalar_laplacian( Laplace_of_Scalar, solution_grid.variables["land_filled"],
                              *coarsened_grid, solution_grid.mask );

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            long double laplace_norm_numer = 0., laplace_norm_denom = 0.;
            double max_filled_scalar = -1e100, min_filled_scalar = 1e100;

            #pragma omp parallel default(none) private(Ipt) \
            shared( coarsened_grid, solution_grid, Laplace_of_Scalar, min_delta ) \
            firstprivate( Npts_coarse, Nlon_coarse, CFL ) \
            reduction( +:laplace_norm_numer,laplace_norm_denom ) \
            reduction( max:max_filled_scalar ) \
            reduction( min:min_filled_scalar )
            {

                laplace_norm_numer = 0.;
                laplace_norm_denom = 0.;

                #pragma omp for collapse(1) schedule(guided)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                    if ( coarsened_grid->mask[Ipt] ) {
                        continue;
                    } else if (Ipt < 2 * Nlon_coarse) {
                        // If we're in the two south-most rings, just set equal to the ring above that
                        size_t ref_point = (Ipt < Nlon_coarse) ? Ipt + 2 * Nlon_coarse : Ipt + Nlon_coarse;
                        #ifdef NDEBUG
                        solution_grid.variables["land_filled"][Ipt] = solution_grid.variables["land_filled"][ref_point];
                        #else
                        solution_grid.variables.at("land_filled").at(Ipt) = solution_grid.variables.at("land_filled").at(ref_point);
                        #endif
                    } else {
                        #ifdef NDEBUG
                        double dScalar_dt = Laplace_of_Scalar[Ipt];
                        #else
                        double dScalar_dt = Laplace_of_Scalar.at(Ipt);
                        #endif

                        double dt = CFL * pow( min_delta[Ipt], 2. );

                        // Euler for Diffusion
                        #ifdef NDEBUG
                        solution_grid.variables["land_filled"][Ipt] += dt * dScalar_dt;
                        #else
                        solution_grid.variables.at("land_filled").at(Ipt) += dt * dScalar_dt;
                        #endif

                        // and accumulate for the 2-norm of the laplacian
                        long double dA = coarsened_grid->areas[Ipt];
                        laplace_norm_numer += dA * pow( (dt/CFL) * dScalar_dt, 2. );
                        laplace_norm_denom += dA;

                        max_filled_scalar = std::max( max_filled_scalar, solution_grid.variables["land_filled"][Ipt] );
                        min_filled_scalar = std::min( min_filled_scalar, solution_grid.variables["land_filled"][Ipt] );
                    }

                }
            }
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            iter_cycle++;
            if ( iter_cycle > Helmholtz_data.iteration_max ) { keep_solving = false; }

            if ( iter_cycle % Helmholtz_data.iterations_per_cycle == 0 ) {
                long double Laplace_2_norm = (laplace_norm_denom == 0) ? -1 : sqrt( laplace_norm_numer / laplace_norm_denom );

                #if DEBUG >= 0
                fprintf( stdout,"  Solver cycle %'zu complete. 2-norm of d(scalar)/dt is %.2Le [sqrt(%.2Le / %.2Le)]. Filled range is %.4g to %.4g\n", 
                        iter_cycle, Laplace_2_norm, laplace_norm_numer, laplace_norm_denom, min_filled_scalar, max_filled_scalar
                       );
                #endif

                const double ref_val = 0.5 * ( std::fabs(min_filled_scalar) + std::fabs(max_filled_scalar) );
                if ( Laplace_2_norm < rel_tol * ref_val ) { keep_solving = false; }
            }
        }

        coarsened_grid->variables.at("land_filled") = solution_grid.variables.at("land_filled");

        std::ostringstream name_stream;
        name_stream << "scalar_landfill_R" << (refine_level) << "_" << (scalar_name) <<  ".nc";
        const std::string output_filename = name_stream.str();
        std::vector<std::string> vars_to_write = { scalar_name, "Laplacian" };
        initialize_output_file( *coarsened_grid, vars_to_write, output_filename, -1 );

        const int ndims = 4;
        size_t starts[ndims] = { 0, 0, 0, 0 };
        size_t counts[ndims] = { 1, 1, size_t(coarsened_grid->Nlat), size_t(coarsened_grid->Nlon) };
        write_field_to_output( coarsened_grid->variables.at("land_filled"), 
                               scalar_name, starts, counts, output_filename, NULL );
        write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
        add_attr_to_file( "CFL", CFL, output_filename );

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
