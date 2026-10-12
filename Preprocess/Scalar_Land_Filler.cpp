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
#include <random>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
//#include <Eigen/QR>
//#include <Eigen/Sparse>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_long_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void Laplace_of_Scalar(
        std::vector<double> & input,
        std::vector<double> & output,
        const dataset * coarsened_grid,
        const std::vector<short int> & mask,
        const std::vector<double> & fixed_water_values,
        std::vector<double> & work1,
        std::vector<double> & work2
            ) {

    size_t Ipt, Npts = fixed_water_values.size();
    Npts = coarsened_grid->variables.at("land_filled").size();
    work1 = fixed_water_values;

    size_t Ifill = 0, Nfill = input.size();
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        if ( not(mask[Ipt]) ) { 
            work1[Ipt] = input[Ifill];
            Ifill++;
        } else {
            work1[Ipt] = coarsened_grid->variables.at("land_filled")[Ipt];
        }
        if ( Ifill >= Nfill ) { break; }
    }

    /*
    // fix land values [lazy approach]
    #pragma omp parallel default(none) private( Ipt ) \
    shared( coarsened_grid, mask, input ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            if ( mask[Ipt] ) { input[Ipt] = coarsened_grid->variables.at("land_filled")[Ipt]; }
        }
    }
    */

    // get the scalar Laplacian
    //scalar_laplacian( output, input, *coarsened_grid, coarsened_grid->mask );
    scalar_laplacian( work2, work1, *coarsened_grid, coarsened_grid->mask );

    // copy to output
    Ifill = 0;
    double max_output = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        if ( not(mask[Ipt]) ) { 
            output[Ifill] = work2[Ipt];
            max_output = std::fmax( max_output, std::fabs(work2[Ipt]) );
            Ifill++;
        }
        if ( Ifill >= Nfill ) { break; }
    }
    fprintf( stdout, "Max output = %g\n", max_output );

    // zero it out over water
    /*
    #pragma omp parallel default(none) private( Ipt ) \
    shared( mask, output ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            if ( mask[Ipt] ) { output[Ipt] = 0.; }
        }
    }
    */
}

class ScalarLaplaceOperator; // declaring existence
namespace Eigen {
    namespace internal {

        template<>
            struct traits<ScalarLaplaceOperator>
            {
                using Scalar = double;
                using StorageKind = Dense;
                using StorageIndex = Eigen::Index;
                using XprKind = MatrixXpr;

                enum {
                    RowsAtCompileTime = Dynamic,
                    ColsAtCompileTime = Dynamic,
                    MaxRowsAtCompileTime = Dynamic,
                    MaxColsAtCompileTime = Dynamic,
                    Flags = 0
                };
            };

    } // namespace internal
} // namespace Eigen



class ScalarLaplaceOperator : public Eigen::EigenBase<ScalarLaplaceOperator>
{
    public:

        using Scalar = double;
        using RealScalar = double;
        using StorageIndex = Eigen::Index;

        enum {
            ColsAtCompileTime = Eigen::Dynamic,
            RowsAtCompileTime = Eigen::Dynamic,
            MaxColsAtCompileTime = Eigen::Dynamic,
            MaxRowsAtCompileTime = Eigen::Dynamic
        };

        ScalarLaplaceOperator(
                const dataset* grid_,
                const std::vector<short int>& mask_,
                const std::vector<double> fixed_water_values_,
                const size_t num_land_
                )
            :
                num_land(num_land_),
                grid(grid_),
                mask(mask_),
                fixed_water_values(fixed_water_values_),
                vin(num_land_),
                vout(num_land_),
                work1(grid_->areas.size()),
                work2(grid_->areas.size())
    {}

        int rows() const { return grid->areas.size(); }
        int cols() const { return grid->areas.size(); }

        //template<typename Vec>
        //    void perform_op(const Vec& x, Vec& y) const
        //void _mul_vector(const Eigen::VectorXd& x,
        //                       Eigen::VectorXd& y) const
        // Matrix-vector product
        Eigen::VectorXd operator*(const Eigen::VectorXd& x) const 
            {

                /* copy Eigen → std::vector */
                for( size_t i=0; i<num_land; i++) {
                    vin[i] = x[i];
                }

                /* apply Laplacian */
                Laplace_of_Scalar(
                        vin,
                        vout,
                        grid,
                        mask,
                        fixed_water_values,
                        work1,
                        work2
                       );

                //const size_t N = grid->areas.size();
                Eigen::VectorXd y(num_land);
                for( size_t i=0; i<num_land; i++) {
                    y[i] = vout[i];
                }
                return y;
            }

    private:

        const size_t num_land;

        const dataset* grid;

        const std::vector<short int>& mask;
        const std::vector<double>& fixed_water_values;

        /* persistent buffers */
        mutable std::vector<double> vin;
        mutable std::vector<double> vout;

        mutable std::vector<double> work1;
        mutable std::vector<double> work2;
};

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
        //std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 

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
        coarsened_grid->variables.at("land_filled") = solution_grid.variables.at("land_filled");

        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        //Helmholtz_data.iteration_max = Npts_coarse * 10;
        //Helmholtz_data.iterations_per_cycle = std::fmin( 10000, Npts_coarse / 100 );
        Helmholtz_data.iteration_max = 10;
        Helmholtz_data.iterations_per_cycle = 10;

        // Seed with zonal average on coarsest grid
        if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
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
        }

        // Impose boundary conditions with one layer of replicate padding
        size_t num_coast = 0;
        double max_boundary_fill = 0, min_boundary_fill = 1e10;
        for (Ipt = 0; Ipt < Npts; Ipt++) {

            if ( coarsened_grid->mask[Ipt]) { 
                // if water, skip
                continue;
            } else {
                double neigh_sum = 0, neigh_count = 0;
                size_t Jpt;
                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    /*
                    // Get coordinates
                    int Itime, Idepth, Ilat, Ilon;
                    Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                    // check each neighbour
                    Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                    Jpt = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                    if ( coarsened_grid->mask[Jpt] ) { 
                        neigh_count++;
                        neigh_sum += coarsened_grid->variables.at("land_filled")[Jpt];
                    }

                    Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                    Jpt = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                    if ( coarsened_grid->mask[Jpt] ) { 
                        neigh_count++;
                        neigh_sum += coarsened_grid->variables.at("land_filled")[Jpt];
                    }

                    if (Ilat > 0) {
                        Jpt = Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon);
                        if ( coarsened_grid->mask[Jpt] ) { 
                            neigh_count++;
                            neigh_sum += coarsened_grid->variables.at("land_filled")[Jpt];
                        }
                    }
                    if (Ilat < Nlat-1) {
                        Jpt = Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon);
                        if ( coarsened_grid->mask[Jpt] ) { 
                            neigh_count++;
                            neigh_sum += coarsened_grid->variables.at("land_filled")[Jpt];
                        }
                    }
                    */
                } else {
                    for (size_t neighbour_ind = 0; neighbour_ind < coarsened_grid->num_neighbours; neighbour_ind++ ) {
                        Jpt = coarsened_grid->adjacency_indices.at(Ipt)[neighbour_ind];
                        if ( coarsened_grid->mask[Jpt] ) { 
                            neigh_count++;
                            neigh_sum += coarsened_grid->variables.at("land_filled")[Jpt];
                        }
                    }
                }

                // coastal points will be updated and flagged as water to fix them
                /*
                if ( neigh_count > 0 ) {
                    coarsened_grid->variables.at("land_filled")[Ipt] = neigh_sum / neigh_count;
                    coarsened_grid->mask[Ipt] = true;

                    num_coast++;
                    max_boundary_fill = std::fmax( max_boundary_fill, std::fabs(neigh_sum / neigh_count) );
                    min_boundary_fill = std::fmin( min_boundary_fill, std::fabs(neigh_sum / neigh_count) );
                }
                */
            }
        }
        fprintf( stdout, "Filled %'zu coastal values. Max fill was %g. Min fill was %g.\n", 
                num_coast, max_boundary_fill, min_boundary_fill );
        solution_grid.variables.at("land_filled") = coarsened_grid->variables.at("land_filled");

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

            if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
                const int ndims = 4;
                size_t starts[ndims] = { 0, 0, 0, 0 };
                size_t counts[ndims] = { 1, 1, size_t(coarsened_grid->Nlat), size_t(coarsened_grid->Nlon) };
                write_field_to_output( solution_grid.variables.at("land_filled"), 
                        scalar_name, starts, counts, output_filename, NULL );
                write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
            } else {
                const int ndims = 3;
                size_t starts[ndims] = { 0, 0, 0 };
                size_t counts[ndims] = { 1, 1, size_t(Npts_coarse) };
                write_field_to_output( solution_grid.variables.at("land_filled"), 
                        scalar_name, starts, counts, output_filename, NULL );
                write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
            }
        }

        size_t num_land = 0;
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
            if ( not( coarsened_grid->mask[Ipt] ) ) { num_land++; }
        }
        fprintf(stdout, "Will solver over the remaining %'zu land points.\n", num_land);

        // Declare the solver
        ScalarLaplaceOperator A( 
                coarsened_grid, 
                coarsened_grid->mask,
                solution_grid.variables.at("land_filled"),
                num_land 
                );
        Eigen::GMRES<ScalarLaplaceOperator, Eigen::IdentityPreconditioner> GMRES_solver;
        GMRES_solver.compute(A);
        GMRES_solver.setTolerance(1e-100);
        GMRES_solver.setMaxIterations(Helmholtz_data.iterations_per_cycle);

        // And make some mappings between std::vector objects and Eigen::VectorXd objects
        /*
        std::vector<double> output_vec(Npts_coarse, 0.), rhs( Npts_coarse, 0. );
        Eigen::Map<Eigen::VectorXd> solver_rhs( rhs.data(), Npts_coarse ),
                                    solver_seed( coarsened_grid->variables.at("land_filled").data(), Npts_coarse ),
                                    solver_output( output_vec.data(),   Npts_coarse );
        */
        std::vector<double> output_vec(num_land, 0.), rhs( num_land, 0. ), seed( num_land, 0 );
        Eigen::Map<Eigen::VectorXd> solver_rhs( rhs.data(), num_land ),
                                    solver_seed( seed.data(), num_land ),
                                    solver_output( output_vec.data(), num_land );


        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(4.0, 40.0);
        for ( size_t Ifill = 0; Ifill < num_land; Ifill++ ) {
            seed[Ifill] = dis(gen);
        }

        #if DEBUG >= 0
        fprintf( stdout, "Beginning solver.\n");
        #endif
        while ( keep_solving ) {

            // apply GMRERS solver
            solver_output = GMRES_solver.solveWithGuess( solver_rhs, solver_seed );
            //solver_output = GMRES_solver.solve( solver_rhs );
            fprintf( stdout, "    GMRES : solver halted after %ld iterations to error of %.5e\n",
                  GMRES_solver.iterations(), GMRES_solver.error() );

            size_t Ifill = 0, Nfill = num_land;
            double max_abs_output = 0.;
            for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                if ( not(coarsened_grid->mask[Ipt]) ) { 
                    //coarsened_grid->variables.at("land_filled")[Ipt] = solver_output[Ifill];
                    //solution_grid.variables.at("land_filled")[Ipt] = solver_output[Ifill];
                    double rand_val = dis(gen);
                    coarsened_grid->variables.at("land_filled")[Ipt] = rand_val;
                    solution_grid.variables.at("land_filled")[Ipt] = rand_val;
                    max_abs_output = std::fmax( max_abs_output, std::fabs( solver_output[Ifill] ) );
                    Ifill++;
                }
                if ( Ifill >= Nfill ) { break; }
            }
            fprintf( stdout, "max(abs(fill)) = %g\n", max_abs_output);
            //std::copy( solver_output.data(), solver_output.data() + solver_output.size(), 
            //           solution_grid.variables.at("land_filled").begin() );


            /*
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
            */

            //iter_cycle++;
            iter_cycle += Helmholtz_data.iterations_per_cycle;
            if ( iter_cycle >= Helmholtz_data.iteration_max ) { keep_solving = false; }

            #if DEBUG >= 0
            fprintf( stdout,"  Solver cycle %'zu complete.\n", iter_cycle );
            #endif

            /*
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
            */
        }

        coarsened_grid->variables.at("land_filled") = solution_grid.variables.at("land_filled");
        scalar_laplacian( Laplace_of_Scalar, solution_grid.variables["land_filled"],
                *coarsened_grid, solution_grid.mask );

        std::ostringstream name_stream;
        name_stream << "scalar_landfill_R" << (refine_level) << "_" << (scalar_name) <<  ".nc";
        const std::string output_filename = name_stream.str();
        std::vector<std::string> vars_to_write = { scalar_name, "Laplacian" };
        initialize_output_file( *coarsened_grid, vars_to_write, output_filename, -1 );

        if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
            const int ndims = 4;
            size_t starts[ndims] = { 0, 0, 0, 0 };
            size_t counts[ndims] = { 1, 1, size_t(coarsened_grid->Nlat), size_t(coarsened_grid->Nlon) };
            write_field_to_output( coarsened_grid->variables.at("land_filled"), 
                                   scalar_name, starts, counts, output_filename, NULL );
            write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
        } else {
            const int ndims = 4;
            size_t starts[ndims] = { 0, 0, 0 };
            size_t counts[ndims] = { 1, 1, Npts_coarse };
            write_field_to_output( coarsened_grid->variables.at("land_filled"), 
                                   scalar_name, starts, counts, output_filename, NULL );
            write_field_to_output( Laplace_of_Scalar, "Laplacian", starts, counts, output_filename, NULL );
        }

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
