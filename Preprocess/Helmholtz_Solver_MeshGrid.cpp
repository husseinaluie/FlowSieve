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

#include <Eigen/Core>
//#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_long_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

void filter_vort_and_div( 
        dataset * coarsened_grid,
        const double ell,
        const bool do_Psi_Phi_instead
        ){

    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   //Ntime  = myCounts.at(0),
                //Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt;

    // some short-hand vars
    std::vector<double> &vort = do_Psi_Phi_instead ? 
                                     coarsened_grid->variables.at("Psi") :
                                     coarsened_grid->variables.at("vort");
    std::vector<double> &div = do_Psi_Phi_instead ? 
                                     coarsened_grid->variables.at("Phi") :
                                     coarsened_grid->variables.at("div");

    const size_t Npts = vort.size();
    const size_t Nvars = 2;

    std::vector<double> coarse_vort( Npts, 0. ), coarse_div( Npts, 0. );

    double dl_kernel_val, dll_kernel_val;
    std::vector<double > filter_values_doubles, filter_dl_values_doubles, filter_dll_values_doubles, 
        local_kernel(Nlat * Nlon, 0.), local_dl_kernel( 0, 0.), local_dll_kernel( 0, 0.);
    std::vector<double*> filter_values_ptrs, filter_dl_values_ptrs, filter_dll_values_ptrs;
    std::vector<const std::vector<double>*> filter_fields;
    filter_fields.push_back( &vort );
    filter_fields.push_back( &div );

    int Itime, Idepth, Ilat, Ilon, LAT_lb, LAT_ub;

    #pragma omp parallel default(none) \
    shared( coarsened_grid, filter_fields, coarse_vort, coarse_div ) \
    private( Ipt, filter_values_doubles, filter_dl_values_doubles, filter_dll_values_doubles, \
                  filter_values_ptrs,    filter_dl_values_ptrs,    filter_dll_values_ptrs, \
             dl_kernel_val, dll_kernel_val, Itime, Idepth, Ilat, Ilon, \
             LAT_lb, LAT_ub, local_kernel, local_dl_kernel, local_dll_kernel ) \
    firstprivate( Npts, ell, Nvars, Nlat, Nlon )
    {
        local_kernel.resize( Nlat*Nlon );
        #pragma omp for collapse(1) schedule(static,1)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            coarsened_grid->index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

            // re initialize vars
            filter_dl_values_doubles.clear();
            filter_dll_values_doubles.clear();
            filter_dl_values_ptrs.clear();
            filter_dll_values_ptrs.clear();

            filter_values_doubles.clear();
            filter_values_ptrs.clear();
            filter_values_doubles.resize( Nvars );
            filter_values_ptrs.resize( Nvars );
            for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) { 
                filter_values_ptrs.at(Ivar) = &(filter_values_doubles.at(Ivar)); 
            }

            // get the kernel
            get_lat_bounds(LAT_lb, LAT_ub, coarsened_grid->latitude, Ilat, ell); 

            // If our longitude grid is uniform, and spans the full periodic domain,
            // then we can just compute it once and translate it at each lon index
            if ( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                std::fill(local_kernel.begin(), local_kernel.end(), 0);
                compute_local_kernel( 
                        local_kernel, local_dl_kernel, local_dll_kernel, 
                        ell, *coarsened_grid, Ilat, 0, LAT_lb, LAT_ub );
            } else {
                throw std::runtime_error("Filter not set up for these settings.");   
            }

            apply_filter_at_point(  
                    filter_values_ptrs, filter_dl_values_ptrs, filter_dll_values_ptrs,
                    dl_kernel_val, dll_kernel_val,
                    filter_fields, *coarsened_grid, Itime, Idepth, Ilat, Ilon, 
                    LAT_lb, LAT_ub, ell, std::vector<bool>(Nvars,false), 
                    local_kernel, local_dl_kernel, local_dll_kernel );

            coarse_vort[Ipt] = filter_values_doubles[0];
            coarse_div[ Ipt] = filter_values_doubles[1];
        }
    }

    // replace fields with coarse versions
    vort = coarse_vort;
    div  = coarse_div;

}

void AM_demodulation( 
        dataset * coarsened_grid,
        const double ell
        ){

    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   //Ntime  = myCounts.at(0),
                //Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt;

    // some short-hand vars
    std::vector<double> &Psi = coarsened_grid->variables.at("Psi");
    std::vector<double> &Phi = coarsened_grid->variables.at("Phi");

    const size_t Npts = Psi.size();
    const size_t Nvars = 2;

    // We're going to multiply by (-1)^(i+j) for logical indices i,j to demodulate
    std::vector<double> demod_Psi( Npts, 0. ), demod_Phi( Npts, 0. );
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        double factor = ( Ipt % 2 == 0 ) ? 1 : -1;
        demod_Psi[Ipt] = factor * Psi[Ipt];
        demod_Phi[Ipt] = factor * Phi[Ipt];
    }

    std::vector<double> coarse_demod_Psi( Npts, 0. ), coarse_demod_Phi( Npts, 0. );

    double dl_kernel_val, dll_kernel_val;
    std::vector<double > filter_values_doubles, filter_dl_values_doubles, filter_dll_values_doubles, 
        local_kernel(Nlat * Nlon, 0.), local_dl_kernel( 0, 0.), local_dll_kernel( 0, 0.);
    std::vector<double*> filter_values_ptrs, filter_dl_values_ptrs, filter_dll_values_ptrs;
    std::vector<const std::vector<double>*> filter_fields;
    filter_fields.push_back( &demod_Psi );
    filter_fields.push_back( &demod_Phi );

    int Itime, Idepth, Ilat, Ilon, LAT_lb, LAT_ub;

    #pragma omp parallel default(none) \
    shared( coarsened_grid, filter_fields, coarse_demod_Psi, coarse_demod_Phi ) \
    private( Ipt, filter_values_doubles, filter_dl_values_doubles, filter_dll_values_doubles, \
                  filter_values_ptrs,    filter_dl_values_ptrs,    filter_dll_values_ptrs, \
             dl_kernel_val, dll_kernel_val, Itime, Idepth, Ilat, Ilon, \
             LAT_lb, LAT_ub, local_kernel, local_dl_kernel, local_dll_kernel ) \
    firstprivate( Npts, ell, Nvars, Nlat, Nlon )
    {
        local_kernel.resize( Nlat*Nlon );
        #pragma omp for collapse(1) schedule(static,1)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            coarsened_grid->index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

            // re initialize vars
            filter_dl_values_doubles.clear();
            filter_dll_values_doubles.clear();
            filter_dl_values_ptrs.clear();
            filter_dll_values_ptrs.clear();

            filter_values_doubles.clear();
            filter_values_ptrs.clear();
            filter_values_doubles.resize( Nvars );
            filter_values_ptrs.resize( Nvars );
            for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) { 
                filter_values_ptrs.at(Ivar) = &(filter_values_doubles.at(Ivar)); 
            }

            // get the kernel
            get_lat_bounds(LAT_lb, LAT_ub, coarsened_grid->latitude, Ilat, ell); 

            // If our longitude grid is uniform, and spans the full periodic domain,
            // then we can just compute it once and translate it at each lon index
            if ( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                std::fill(local_kernel.begin(), local_kernel.end(), 0);
                compute_local_kernel( 
                        local_kernel, local_dl_kernel, local_dll_kernel, 
                        ell, *coarsened_grid, Ilat, 0, LAT_lb, LAT_ub );
            } else {
                throw std::runtime_error("Filter not set up for these settings.");   
            }

            apply_filter_at_point(  
                    filter_values_ptrs, filter_dl_values_ptrs, filter_dll_values_ptrs,
                    dl_kernel_val, dll_kernel_val,
                    filter_fields, *coarsened_grid, Itime, Idepth, Ilat, Ilon, 
                    LAT_lb, LAT_ub, ell, std::vector<bool>(Nvars,false), 
                    local_kernel, local_dl_kernel, local_dll_kernel );

            coarse_demod_Psi[Ipt] = filter_values_doubles[0];
            coarse_demod_Phi[Ipt] = filter_values_doubles[1];
        }
    }

    // replace fields with coarse versions
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        double factor = ( Ipt % 2 == 0 ) ? 1 : -1;
        Psi[Ipt] = Psi[Ipt] - factor * coarse_demod_Psi[Ipt];
        Phi[Ipt] = Phi[Ipt] - factor * coarse_demod_Phi[Ipt];
    }

}

void Laplace_of_field(
        std::vector<double> & input,
        std::vector<double> & output,
        const dataset * coarsened_grid,
        const std::vector<double> & cos_ip,
        const std::vector<double> & cos_im,
        std::vector<double> & work_1,
        std::vector<double> & work_2
            ) {

    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    // for now, assume uniform grid
    const double dlon_inv = 1. / ( longitude[1] - longitude[0]);
    const double dlat_inv = 1. / ( latitude[1]  - latitude[0] );


    const std::vector<int>  &myCounts = coarsened_grid->myCounts;
    const int   //Ntime  = myCounts.at(0),
                //Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt;

    // some short-hand vars
    std::vector<double> &proj_uo = work_1,
                        &proj_vo = work_2;

    const size_t Npts = input.size();

    // get mean polar values of Psi and Phi
    double S_pole_Phi = 0, N_pole_Phi = 0.;
    const double per_Nlon = 1. / Nlon;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

        S_pole_Phi += input[Ipt_S] * per_Nlon;
        N_pole_Phi += input[Ipt_N] * per_Nlon;
    }

    /*
    // Set polar rings to mean value
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

        input[Ipt_S] = S_pole_Phi;
        input[Ipt_N] = N_pole_Phi;
    }
    */

    // zero out the fields
    #pragma omp parallel default(none) private( Ipt ) \
    shared( proj_uo, proj_vo, output ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            #ifdef NDEBUG
            proj_uo[Ipt] = 0; 
            proj_vo[Ipt] = 0;  
            output[ Ipt] = 0;
            #else
            proj_uo.at(Ipt) = 0; 
            proj_vo.at(Ipt) = 0;  
            output.at( Ipt) = 0;
            #endif
        }
    }

    const double R_E_inv = 1. / constants::R_earth;


    //
    // First loop set gets the velocity components
    //


    #pragma omp parallel default(none) private( Ipt ) \
    shared( coarsened_grid, input, proj_uo, proj_vo ) \
    firstprivate( R_E_inv, Npts, Nlat, Nlon, dlat_inv, dlon_inv, N_pole_Phi, S_pole_Phi )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );

            // E-W indices
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // N-S indices
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );


            #ifdef NDEBUG
            const double Phi_N = ( Ilat < Nlat - 1 ) ? input[N_ind] : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? input[S_ind] : S_pole_Phi,
                         Phi_E = input[E_ind],
                         Phi_W = input[W_ind];
            #else
            const double Phi_N = ( Ilat < Nlat - 1 ) ? input.at(N_ind) : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? input.at(S_ind) : S_pole_Phi,
                         Phi_E = input.at(E_ind),
                         Phi_W = input.at(W_ind);
            #endif

            // NOTE we are computing cos(lat)*uo and cos(lat)*vo for stability
            proj_uo[Ipt] = ( Phi_E - Phi_W ) * R_E_inv *           dlon_inv * 0.5;
            proj_vo[Ipt] = ( Phi_N - Phi_S ) * R_E_inv * cos_lat * dlat_inv * 0.5;

            //proj_uo[Ipt] = ( Phi_E - Phi_W ) * R_E_cos_lat_inv * dlon_inv * 0.5;
            //proj_vo[Ipt] = ( Phi_N - Phi_S ) * R_E_inv         * dlat_inv * 0.5;
        }
    }


    //
    // Second loop sets gets vorticity and divergence
    //

    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Nlat, Nlon, Npts, dlat_inv, dlon_inv )\
    shared( coarsened_grid, output, proj_uo, proj_vo )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            /*
            double cos_plus = 
                ( Ilat == Nlat - 1 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat+1] :
                cos( coarsened_grid->latitude[Ilat+1] );
            double cos_minus = 
                ( Ilat == 0 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat-1] :
                cos( coarsened_grid->latitude[Ilat-1] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            */

            // E-W indices
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // N-S indices
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_E = proj_uo[E_ind];
            double uo_W = proj_uo[W_ind];

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo[N_ind] : 0.;
            double vo_S = ( Ilat > 0        ) ? proj_vo[S_ind] : 0.;
            #else
            double uo_E = proj_uo.at(E_ind);
            double uo_W = proj_uo.at(W_ind);

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo.at(N_ind) : 0.;
            double vo_S = ( Ilat > 0        ) ? proj_vo.at(S_ind) : 0.;
            #endif

            // Divergence
            double tmp_div =    ( uo_E - uo_W ) * dlon_inv * 0.5
                              + ( vo_N - vo_S ) * dlat_inv * 0.5 * cos_lat;
                            //+ ( cos_plus * vo_N - cos_minus * vo_S ) * dlat_inv * 0.5;

            // scaling
            tmp_div *= R_E_inv; // we are getting (cos(lat))^2 * Lap, for stability
            //tmp_div *= R_E_cos_lat_inv;
            //output[Ipt] = tmp_div;

            // area weight the Laplace to area-bias the solver
            double dA = 1.; //coarsened_grid->areas[Ipt];

            #ifdef NDEBUG
            output[Ipt] = tmp_div * dA;
            #else
            output.at(Ipt) = tmp_div * dA;
            #endif
        }
    }

}

class LaplaceOperator; // declaring existence
namespace Eigen {
    namespace internal {

        template<>
            struct traits<LaplaceOperator>
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



class LaplaceOperator : public Eigen::EigenBase<LaplaceOperator>
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

        LaplaceOperator(
                const dataset* grid_,
                const std::vector<double>& cos_ip_,
                const std::vector<double>& cos_im_
                )
            :
                grid(grid_),
                cos_ip(cos_ip_),
                cos_im(cos_im_),
                work1(grid_->areas.size()),
                work2(grid_->areas.size()),
                vin(grid_->areas.size()),
                vout(grid_->areas.size())
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
                const size_t N = grid->areas.size();

                /* copy Eigen → std::vector */
                for( size_t i=0; i<N; i++) {
                    vin[i] = x[i];
                }

                /* apply Laplacian */
                Laplace_of_field(
                        vin,
                        vout,
                        grid,
                        cos_ip,
                        cos_im,
                        work1,
                        work2
                       );

                Eigen::VectorXd y(N);
                for( size_t i=0; i<N; i++) {
                    y[i] = vout[i];
                }
                return y;
            }

    private:

        const dataset* grid;

        const std::vector<double>& cos_ip;
        const std::vector<double>& cos_im;

        /* persistent buffers */

        mutable std::vector<double> work1;
        mutable std::vector<double> work2;

        mutable std::vector<double> vin;
        mutable std::vector<double> vout;
};

void get_true_vort_div_Eigen(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid
        ) {

    //
    // get vorticity and divergence
    //

    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   //Ntime  = myCounts.at(0),
                //Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    std::vector<double> &uo = coarsened_grid->variables.at("u_lon"),
                        &vo = coarsened_grid->variables.at("u_lat"),
                        &vort = coarsened_grid->variables.at("vort"),
                        &div  = coarsened_grid->variables.at("div" );

    const size_t Npts = div.size();

    // for now, assume uniform grid
    const double dlon_inv = 1. / ( longitude[1] - longitude[0]);
    const double dlat_inv = 1. / ( latitude[1]  - latitude[0] );
    const double R_E_inv = 1. / constants::R_earth;

    size_t Ipt;

    // get mean polar values of vo and uo
    /*
    double S_pole_vo = 0, N_pole_vo = 0., S_pole_uo = 0, N_pole_uo = 0;
    const double per_Nlon = 1. / Nlon;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0,      Ilon ); 

        S_pole_uo += uo[Ipt_S] * per_Nlon;
        N_pole_uo += uo[Ipt_N] * per_Nlon;

        S_pole_vo += vo[Ipt_S] * per_Nlon;
        N_pole_vo += vo[Ipt_N] * per_Nlon;
    }
    */

    //
    // Second loop sets gets vorticity and divergence
    //
    // this time using south and west neighbours

    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Nlat, Nlon, Npts, dlat_inv, dlon_inv )\
    shared( coarsened_grid, PsiPhi_grid, vort, div, uo, vo )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            /*
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double cos_plus = 
                ( Ilat == Nlat - 1 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat+1] :
                cos( coarsened_grid->latitude[Ilat+1] );
            double cos_minus = 
                ( Ilat == 0 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat-1] :
                cos( coarsened_grid->latitude[Ilat-1] );
            */

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_N = ( Ilat < Nlat - 1 ) ? uo[N_ind] : 0.; //N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? uo[S_ind] : 0.; //S_pole_uo;
            double uo_E = uo[E_ind];
            double uo_W = uo[W_ind];

            double vo_N = ( Ilat < Nlat - 1 ) ? vo[N_ind] : 0.; //N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? vo[S_ind] : 0.; //S_pole_vo;
            double vo_E = vo[E_ind];
            double vo_W = vo[W_ind];
            #else
            double uo_N = ( Ilat < Nlat - 1 ) ? uo.at(N_ind) : 0.; //N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? uo.at(S_ind) : 0.; //S_pole_uo;
            double uo_E = uo.at(E_ind);
            double uo_W = uo.at(W_ind);

            double vo_N = ( Ilat < Nlat - 1 ) ? vo.at(N_ind) : 0.; //N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? vo.at(S_ind) : 0.; //S_pole_vo;
            double vo_E = vo.at(E_ind);
            double vo_W = vo.at(W_ind);
            #endif

            // NOTE: velocitie have been scaled by cos(lat)

            // Divergence
            double tmp_div =    ( uo_E - uo_W ) * dlon_inv * 0.5
                              + ( vo_N - vo_S ) * dlat_inv * 0.5 * cos_lat;
            //double tmp_div =    (            uo_E -             uo_W ) * dlon_inv * 0.5;
            //                  + ( cos_plus * vo_N - cos_minus * vo_S ) * dlat_inv * 0.5;

            // Vorticity
            double tmp_vort =   ( vo_E - vo_W ) * dlon_inv * 0.5
                              - ( uo_N - uo_S ) * dlat_inv * 0.5 * cos_lat;
            //double tmp_vort =   (            vo_E -             vo_W ) * dlon_inv * 0.5;
            //                  - ( cos_plus * uo_N - cos_minus * uo_S ) * dlat_inv * 0.5;

            // apply geometric scaling
            //tmp_div  *= R_E_cos_lat_inv;
            //tmp_vort *= R_E_cos_lat_inv;
            tmp_div  *= R_E_inv;// * cos_lat; // getting (cos(lat))^2 * div
            tmp_vort *= R_E_inv;// * cos_lat; // getting (cos(lat))^2 * vort
            

            #ifdef NDEBUG
            vort[Ipt] = tmp_vort;
            div[ Ipt] = tmp_div;
            #else
            vort.at(Ipt) = tmp_vort;
            div.at( Ipt) = tmp_div;
            #endif
        }
    }
}

void get_vels_vort_and_div_Eigen(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid,
        const bool cos_prefactor
            ) {

    // Similar to get_vels_vort_and_div, except this manually implements the derivatives
    // By lumping all of the derivative calculations together we can be more efficient, especially
    // with the threading.
    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    // for now, assume uniform grid
    const double dlon_inv = 1. / ( longitude[1] - longitude[0]);
    const double dlat_inv = 1. / ( latitude[1]  - latitude[0] );


    const std::vector<int>  &myCounts = coarsened_grid->myCounts;
    const int   //Ntime  = myCounts.at(0),
                //Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt;

    // some short-hand vars
    std::vector<double> &proj_uo_tor = coarsened_grid->variables.at("proj_uo_tor"),
                        &proj_uo_pot = coarsened_grid->variables.at("proj_uo_pot"),
                        &proj_vo_tor = coarsened_grid->variables.at("proj_vo_tor"),
                        &proj_vo_pot = coarsened_grid->variables.at("proj_vo_pot"),
                        &proj_vort   = coarsened_grid->variables.at("proj_vort"),
                        &proj_div    = coarsened_grid->variables.at("proj_div" ),
                        &Phi_vort    = coarsened_grid->variables.at("proj_vort_of_pot"),
                        &Psi_div     = coarsened_grid->variables.at("proj_div_of_tor" ),
                        &Psi         = PsiPhi_grid.variables["Psi"],
                        &Phi         = PsiPhi_grid.variables["Phi"];

    const size_t Npts = Psi.size();

    // get mean polar values of Psi and Phi
    double S_pole_Phi = 0, N_pole_Phi = 0., S_pole_Psi = 0, N_pole_Psi = 0;
    const double per_Nlon = 1. / Nlon;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0,      Ilon ); 

        S_pole_Phi += Phi[Ipt_S] * per_Nlon;
        N_pole_Phi += Phi[Ipt_N] * per_Nlon;

        S_pole_Psi += Psi[Ipt_S] * per_Nlon;
        N_pole_Psi += Psi[Ipt_N] * per_Nlon;
    }

    // zero out the fields
    #pragma omp parallel default(none) private( Ipt ) \
    shared( proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot, proj_vort, proj_div ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            #ifdef NDEBUG
            proj_uo_tor[Ipt] = 0; 
            proj_vo_tor[Ipt] = 0;  
            proj_uo_pot[Ipt] = 0;
            proj_vo_pot[Ipt] = 0;

            proj_vort[Ipt] = 0;
            proj_div[ Ipt] = 0;
            #else
            proj_uo_tor.at(Ipt) = 0; 
            proj_vo_tor.at(Ipt) = 0;  
            proj_uo_pot.at(Ipt) = 0;
            proj_vo_pot.at(Ipt) = 0;

            proj_vort.at(Ipt) = 0;
            proj_div.at( Ipt) = 0;
            #endif
        }
    }

    const double R_E_inv = 1. / constants::R_earth;


    //
    // First loop set gets the velocity components
    //



    #pragma omp parallel default(none) private( Ipt ) \
    shared( coarsened_grid, PsiPhi_grid, Psi, Phi, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot ) \
    firstprivate( R_E_inv, Npts, Nlat, Nlon, dlat_inv, dlon_inv,\
           N_pole_Psi, N_pole_Phi, S_pole_Psi, S_pole_Phi, cos_prefactor )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            /*
            double cos_plus = 
                ( Ilat == Nlat - 1 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat+1] :
                cos( coarsened_grid->latitude[Ilat+1] );
            double cos_minus = 
                ( Ilat == 0 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat-1] :
                cos( coarsened_grid->latitude[Ilat-1] );
            */

            // E-W indices
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // N-S indices
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );


            #ifdef NDEBUG
            const double Psi_N = ( Ilat < Nlat - 1 ) ? Psi[N_ind] : N_pole_Psi,
                         Psi_S = ( Ilat > 0        ) ? Psi[S_ind] : S_pole_Psi,
                         Psi_E = Psi[E_ind],
                         Psi_W = Psi[W_ind];
            const double Phi_N = ( Ilat < Nlat - 1 ) ? Phi[N_ind] : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? Phi[S_ind] : S_pole_Phi,
                         Phi_E = Phi[E_ind],
                         Phi_W = Phi[W_ind];
            #else
            const double Psi_N = ( Ilat < Nlat - 1 ) ? Psi.at(N_ind) : N_pole_Psi,
                         Psi_S = ( Ilat > 0        ) ? Psi.at(S_ind) : S_pole_Psi,
                         Psi_E = Psi.at(E_ind),
                         Psi_W = Psi.at(W_ind);
            const double Phi_N = ( Ilat < Nlat - 1 ) ? Phi.at(N_ind) : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? Phi.at(S_ind) : S_pole_Phi,
                         Phi_E = Phi.at(E_ind),
                         Phi_W = Phi.at(W_ind);
            #endif


            if ( cos_prefactor ) {
                // NOTE: computing cos(lat) * vel for stability
                proj_uo_tor[Ipt] = - ( Psi_N - Psi_S ) * R_E_inv * cos_lat * dlat_inv * 0.5;
                proj_vo_tor[Ipt] =   ( Psi_E - Psi_W ) * R_E_inv *           dlon_inv * 0.5;

                proj_uo_pot[Ipt] =   ( Phi_E - Phi_W ) * R_E_inv           * dlon_inv * 0.5;
                proj_vo_pot[Ipt] =   ( Phi_N - Phi_S ) * R_E_inv * cos_lat * dlat_inv * 0.5;
            } else {
                proj_uo_tor[Ipt] = - ( Psi_N - Psi_S ) * R_E_inv         * dlat_inv * 0.5;
                proj_vo_tor[Ipt] =   ( Psi_E - Psi_W ) * R_E_cos_lat_inv * dlon_inv * 0.5;

                proj_uo_pot[Ipt] =   ( Phi_E - Phi_W ) * R_E_cos_lat_inv * dlon_inv * 0.5;
                proj_vo_pot[Ipt] =   ( Phi_N - Phi_S ) * R_E_inv         * dlat_inv * 0.5;
            }
        }
    }

    // get mean polar values of vo and uo
    //double S_pole_vo = 0, N_pole_vo = 0., S_pole_uo = 0, N_pole_uo = 0;
    /*
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

        S_pole_uo += (proj_uo_tor[Ipt_S] + proj_uo_pot[Ipt_S]) * per_Nlon;
        N_pole_uo += (proj_uo_tor[Ipt_N] + proj_uo_pot[Ipt_N]) * per_Nlon;

        S_pole_vo += (proj_vo_tor[Ipt_S] + proj_vo_pot[Ipt_S]) * per_Nlon;
        N_pole_vo += (proj_vo_tor[Ipt_N] + proj_vo_pot[Ipt_N]) * per_Nlon;
    }
    */


    //
    // Second loop sets gets vorticity and divergence
    //
    // this time using south and west neighbours

    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Nlat, Nlon, Npts, dlat_inv, dlon_inv ) \
    shared( coarsened_grid, PsiPhi_grid, proj_vort, proj_div, Phi_vort, Psi_div, proj_uo_tor, \
           proj_vo_tor, proj_uo_pot, proj_vo_pot, cos_prefactor )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double cos_plus = 
                ( Ilat == Nlat - 1 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat+1] :
                cos( coarsened_grid->latitude[Ilat+1] );
            double cos_minus = 
                ( Ilat == 0 ) ? 0 :
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat-1] :
                cos( coarsened_grid->latitude[Ilat-1] );

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot[N_ind] + proj_uo_tor[N_ind] : 0.; //N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? proj_uo_pot[S_ind] + proj_uo_tor[S_ind] : 0.; //S_pole_uo;
            double uo_E = proj_uo_pot[E_ind] + proj_uo_tor[E_ind];
            double uo_W = proj_uo_pot[W_ind] + proj_uo_tor[W_ind];

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_pot[N_ind] + proj_vo_tor[N_ind] : 0.; //N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? proj_vo_pot[S_ind] + proj_vo_tor[S_ind] : 0.; //S_pole_vo;
            double vo_E = proj_vo_pot[E_ind] + proj_vo_tor[E_ind];
            double vo_W = proj_vo_pot[W_ind] + proj_vo_tor[W_ind];
            #else
            double uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot.at(N_ind) + proj_uo_tor.at(N_ind) : 0.; //N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? proj_uo_pot.at(S_ind) + proj_uo_tor.at(S_ind) : 0.; //S_pole_uo;
            double uo_E = proj_uo_pot.at(E_ind) + proj_uo_tor.at(E_ind);
            double uo_W = proj_uo_pot.at(W_ind) + proj_uo_tor.at(W_ind);

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_pot.at(N_ind) + proj_vo_tor.at(N_ind) : 0.; //N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? proj_vo_pot.at(S_ind) + proj_vo_tor.at(S_ind) : 0.; //S_pole_vo;
            double vo_E = proj_vo_pot.at(E_ind) + proj_vo_tor.at(E_ind);
            double vo_W = proj_vo_pot.at(W_ind) + proj_vo_tor.at(W_ind);
            #endif

            double tmp_div = 0, tmp_vort = 0;
            if ( cos_prefactor ) {
                // Divergence
                tmp_div =   ( uo_E - uo_W ) * dlon_inv * 0.5
                          + ( vo_N - vo_S ) * dlat_inv * 0.5 * cos_lat;

                // Vorticity
                tmp_vort =   ( vo_E - vo_W ) * dlon_inv * 0.5
                           - ( uo_N - uo_S ) * dlat_inv * 0.5 * cos_lat;

                // geometric scalings
                tmp_div  *= R_E_inv; // getting (cos(lat))^2 * div
                tmp_vort *= R_E_inv; // getting (cos(lat))^2 * vort
            } else {
                // Divergence
                tmp_div =   (            uo_E -             uo_W ) * dlon_inv * 0.5
                          + ( cos_plus * vo_N - cos_minus * vo_S ) * dlat_inv * 0.5;

                // Vorticity
                tmp_vort =   (            vo_E -             vo_W ) * dlon_inv * 0.5
                           - ( cos_plus * uo_N - cos_minus * uo_S ) * dlat_inv * 0.5;

                // geometric scalings
                tmp_div  *= R_E_cos_lat_inv;
                tmp_vort *= R_E_cos_lat_inv;
            }


            #ifdef NDEBUG
            proj_vort[Ipt] = tmp_vort;
            proj_div[ Ipt] = tmp_div;
            #else
            proj_vort.at(Ipt) = tmp_vort;
            proj_div.at( Ipt) = tmp_div;
            #endif


            //
            //// And repeat for get div(tor_vel) and vort(pot_vel) for sanity checking
            //
            #ifdef NDEBUG
            uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot[N_ind] : 0.; //N_pole_uo;
            uo_S = ( Ilat > 0        ) ? proj_uo_pot[S_ind] : 0.; //S_pole_uo;
            uo_E = proj_uo_tor[E_ind];
            uo_W = proj_uo_tor[W_ind];

            vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_tor[N_ind] : 0.; //N_pole_vo;
            vo_S = ( Ilat > 0        ) ? proj_vo_tor[S_ind] : 0.; //S_pole_vo;
            vo_E = proj_vo_pot[E_ind];
            vo_W = proj_vo_pot[W_ind];
            #else
            uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot.at(N_ind) : 0.; //N_pole_uo;
            uo_S = ( Ilat > 0        ) ? proj_uo_pot.at(S_ind) : 0.; //S_pole_uo;
            uo_E = proj_uo_tor.at(E_ind);
            uo_W = proj_uo_tor.at(W_ind);

            vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_tor.at(N_ind) : 0.; //N_pole_vo;
            vo_S = ( Ilat > 0        ) ? proj_vo_tor.at(S_ind) : 0.; //S_pole_vo;
            vo_E = proj_vo_pot.at(E_ind);
            vo_W = proj_vo_pot.at(W_ind);
            #endif

            if ( cos_prefactor ) {
                // Divergence
                tmp_div =    ( uo_E - uo_W ) * dlon_inv * 0.5
                    + ( vo_N - vo_S ) * dlat_inv * 0.5 * cos_lat;

                // Vorticity
                tmp_vort =   ( vo_E - vo_W ) * dlon_inv * 0.5
                    - ( uo_N - uo_S ) * dlat_inv * 0.5 * cos_lat;

                // geometric scalings
                tmp_div  *= R_E_inv; // getting (cos(lat))^2 * div
                tmp_vort *= R_E_inv; // getting (cos(lat))^2 * vort
            } else {
                // Divergence
                tmp_div =    (            uo_E -             uo_W ) * dlon_inv * 0.5
                           + ( cos_plus * vo_N - cos_minus * vo_S ) * dlat_inv * 0.5;

                // Vorticity
                tmp_vort =   (            vo_E -             vo_W ) * dlon_inv * 0.5
                           - ( cos_plus * uo_N - cos_minus * uo_S ) * dlat_inv * 0.5;

                // geometric scalings
                tmp_div  *= R_E_cos_lat_inv;
                tmp_vort *= R_E_cos_lat_inv;
            }

            #ifdef NDEBUG
            Phi_vort[Ipt] = tmp_vort;
            Psi_div[ Ipt] = tmp_div;
            #else
            Phi_vort.at(Ipt) = tmp_vort;
            Psi_div.at( Ipt) = tmp_div;
            #endif
        }
    }


}

void Helmholtz_Solver_MeshGrid(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const int num_refinements,
        const bool seed_provided,
        const std::vector<double> & Psi_seed,
        const std::vector<double> & Phi_seed,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    const size_t Npts = source_data.mask.size();

    /*
    const int Ntime = source_data.Ntime;
    const int Ndepth = source_data.Ndepth;
    */

    // Fill in the land areas with zero velocity
    size_t index, Ipt;
    #pragma omp parallel default(none) \
    shared( source_data ) \
    private( index ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; index++) {
            if (not(source_data.mask[index])) {
                source_data.variables["u_lon"][index] = 0.;
                source_data.variables["u_lat"][index] = 0.;
            }
        }
    }

    double mean_uo = get_mean( source_data.variables.at("u_lon"), &source_data );
    double mean_vo = get_mean( source_data.variables.at("u_lat"), &source_data );
    fprintf( stdout, "Mean values of zonal and meridional velocity are %.2e, %.2e\n", mean_uo, mean_vo );

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    dataset solution_grid;
    dataset *coarsened_grid = new dataset();

    // Beginning the v-cycle [or, more acurrately, just a / cycle?]
    // Start at the coarsest resolution, solve, refine, repeat
    int Nlat_coarse = 0, Nlon_coarse = 0;
    size_t Npts_coarse = 0;
    for ( int refine_level = num_refinements; refine_level >= 0; refine_level-- ) {

        #if DEBUG >= 0
        fprintf(stdout, "\n\nBeginning refinement level %d.\n", refine_level);
        #endif
        std::vector<short int> saved_mask;

        // If refinement is turned on, do it
        if ( refine_level > 0 ) {
            Nlat_coarse = floor(sqrt( 0.5 * (Npts / pow(4., 2*refine_level )) ));
            Nlon_coarse = 2 * (Nlat_coarse + 1);
            Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;

            #if DEBUG >= 2
            fprintf(stdout, "Refinement level has Nlat x Nlon of %d x %d.\n", Nlat_coarse, Nlon_coarse);
            #endif


            // Downsample the velocities onto the coarse grid
            std::vector<std::string> vars_to_map = { "u_lon", "u_lat", "vort", "div" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            initialize_coarsened_grid( *coarsened_grid, source_data, Nlat_coarse, Nlon_coarse );
            coarsened_grid->initialize_index_functions();
            coarsened_grid->compute_sin_cos_lat();
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "InitializeCoarseGrid"); }

            Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;
            for ( size_t Ivar = 0; Ivar < vars_to_map.size(); Ivar++ ) {
                coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                            vars_to_map[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
            }

            #if DEBUG >= 2
            fprintf(stdout, "Mapping velocities onto the coarse grid.\n");
            #endif
            vars_to_map = { "u_lon", "u_lat" };
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( source_data, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

            coarsened_grid->adjacency_indices_adjoint.resize(0);


        } else {
            // Otherwise we're on the original grid, so just use the source
            coarsened_grid = &source_data;
            Npts_coarse = Npts;
            Nlat_coarse = coarsened_grid->Nlat;
            Nlon_coarse = coarsened_grid->Nlon;

            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "vort", std::vector<double>(Npts_coarse, 0.) ) );
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "div", std::vector<double>(Npts_coarse, 0.) ) );

            coarsened_grid->initialize_index_functions();
            coarsened_grid->compute_sin_cos_lat();
            coarsened_grid->compute_cell_areas();

        }

        fprintf( stdout, "Npts_coarse = %'zu, len(areas) = %'zu\n", Npts_coarse, coarsened_grid->areas.size() );

        // We're done with the land mask, so remove it. This causes vort and
        // div to be computed on land cells to ensure agreement
        saved_mask = coarsened_grid->mask; // store a copy
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );
        
        // Solver coefficients
        std::vector<double> cos_ip(Nlat_coarse, 0.),
                            cos_im(Nlat_coarse, 0.);
        //const double R_E_sqr = pow( constants::R_earth, 2. );

        // Get vorticity and divergence
        #if DEBUG >= 2
        fprintf(stdout, "Computing vorticity and divergence on coarse grid.\n");
        #endif
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        // apply cosine weighting to vels
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            coarsened_grid->variables.at("u_lon")[Ipt] = 
                cos_lat * coarsened_grid->variables.at("u_lon")[Ipt];
            coarsened_grid->variables.at("u_lat")[Ipt] = 
                cos_lat * coarsened_grid->variables.at("u_lat")[Ipt];
        }

        get_true_vort_div_Eigen( coarsened_grid, *coarsened_grid );

        double mean_vort = get_mean( coarsened_grid->variables.at("vort"), coarsened_grid );
        double mean_div  = get_mean( coarsened_grid->variables.at("div" ), coarsened_grid );
        /*
        double mean_vort, mean_div;
        demean( coarsened_grid->variables.at("vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("div" ), coarsened_grid, mean_div  );
        */
        fprintf( stdout, "Mean (coarse-grid) values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

        // We need to de-noise the vorticity and divergence
        //  we don't believe the grid-scale features on this grid anyways
        // We'll set the filter scale (kernel diameter) at 2*delta, where delta is the
        //   latitudinal grid spacing
        /*
        if ( refine_level > 0 ) {
            const double filter_ell = 2 * constants::R_earth * ( M_PI / Nlat_coarse );
            fprintf( stdout, "Filtering vorticity and divergence at %.0fkm before projecting.\n", 
                    filter_ell / 1e3 );
            filter_vort_and_div( coarsened_grid, filter_ell, false );
        }
        */

        double mean_uo_coarse = get_mean( coarsened_grid->variables.at("u_lon"), coarsened_grid );
        double mean_vo_coarse = get_mean( coarsened_grid->variables.at("u_lat"), coarsened_grid );
        fprintf( stdout, "Mean (coarse-grid) values of zonal and meridional velocity are %.2e, %.2e\n", mean_uo_coarse, mean_vo_coarse );

        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

        // Add Psi and Phi variables to our coarsened grid
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );

        // Now, up-sample the previous solutions onto the new grid
        if ( refine_level < num_refinements ) {

            std::vector<std::string> vars_to_map = { "Psi", "Phi" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( solution_grid, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Upsampling"); }

            // Filter our grid-scale nonsense from the previous solution
            //  note: *4 is because we're filtering based on the coarser grid (4x coarsening)
            /*
            const double filter_ell = 2 * 4 * constants::R_earth * ( M_PI / Nlat_coarse );
            coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
            coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");
            filter_vort_and_div( coarsened_grid, filter_ell, true );
            solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
            solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
            */
        }

        // Re-initialize the solution grid on the new computational grid
        // To avoid having to re-build, just copy from the coarsened grid
        solution_grid.clear();
        solution_grid.copy_from_ptr( coarsened_grid );
        solution_grid.initialize_index_functions();
        solution_grid.compute_sin_cos_lat();
        std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 

        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );


        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.tolerance = rel_tol;

        //Helmholtz_data.iterations_per_cycle = Nlon_coarse *   10;
        //Helmholtz_data.iterations_per_cycle = Nlon_coarse *   1;
        //Helmholtz_data.iteration_max        = Nlon_coarse * 1000;
        //Helmholtz_data.iterations_per_cycle = 500;
        //Helmholtz_data.iteration_max        = 50*500;
        Helmholtz_data.iterations_per_cycle = 20000;
        Helmholtz_data.iteration_max        = 
            ( refine_level == num_refinements ) ? Helmholtz_data.iterations_per_cycle
            : 5*Helmholtz_data.iterations_per_cycle;

        // Use spherical harmonics to set the seed
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );
        if ( (num_refinements == 0) and (seed_provided) ) {
            coarsened_grid->variables["Psi"] = Psi_seed;
            coarsened_grid->variables["Phi"] = Phi_seed;
        } else if ( refine_level == num_refinements ) {
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            SphericalHarmonicSolver(  
                    coarsened_grid->variables.at("Psi"),
                    coarsened_grid->variables.at("Phi"),
                   *coarsened_grid,
                    coarsened_grid->variables.at("vort"),
                    coarsened_grid->variables.at("div") 
                    );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "SphericalHarmonicSolver"); }
        }


        //
        //// Output the seed / starting guess
        //
        std::vector<std::string> vars_to_add = { 
            "projected_u_lon", "projected_u_lat", "projected_vort", "projected_div",
            "proj_uo_tor", "proj_uo_pot", "proj_vo_tor", "proj_vo_pot" ,
            "proj_vort", "proj_div"
        };
        for ( size_t Ivar = 0; Ivar < vars_to_add.size(); Ivar++ ) {
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        vars_to_add[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
        }

        // As a sanity check, get the divergence of the toroidal field and vorticity of the potential field
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "proj_vort_of_pot", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "proj_div_of_tor", std::vector<double>(Npts_coarse, 0.) ) );

        // Compute velocity components
        get_vels_vort_and_div_Eigen( coarsened_grid, *coarsened_grid, true );
        //demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
        //demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );

        for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
            coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                  coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
            coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                  coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
        }

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div")  = coarsened_grid->variables.at("proj_div");

        // Write the seed to output
        std::ostringstream seed_stream;
        seed_stream << "seed_R" << (refine_level) << ".nc";
        const std::string seed_filename = seed_stream.str();
        solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
        solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
        write_Helmholtz_output( seed_filename, Helmholtz_data, *coarsened_grid, solution_grid );

        // Make some short-hand references
        std::vector<double> &soln_Psi = solution_grid.variables.at("Psi"),
                            &soln_Phi = solution_grid.variables.at("Phi"); 

        std::vector<double> //&cgrid_proj_vort = coarsened_grid->variables.at("proj_vort"),
                            //&cgrid_proj_div  = coarsened_grid->variables.at("proj_div"),
                            &cgrid_vort = coarsened_grid->variables.at("vort"),
                            &cgrid_div  = coarsened_grid->variables.at("div");


        // compute errors and check convergence
        long double vel_2_err  = 0, vel_inf_err  = 0,
                    vort_2_err = 0, vort_inf_err = 0,
                    div_2_err  = 0, div_inf_err  = 0;
        long double vel_2_norm  = 0, vel_inf_norm  = 0,
                    vort_2_norm = 0, vort_inf_norm = 0,
                    div_2_norm  = 0, div_inf_norm  = 0;
        long double vort_2_ener = 0, div_2_ener = 0;
        long double Psi_2_norm = 0, Phi_2_norm = 0;

        get_norms_GS(  vel_2_err,  vort_2_err,  div_2_err,
                vel_2_norm, vort_2_norm, div_2_norm,
                vel_inf_err,  vort_inf_err,  div_inf_err,
                vel_inf_norm, vort_inf_norm, div_inf_norm,
                vort_2_ener, div_2_ener,
                Psi_2_norm, Phi_2_norm,
                coarsened_grid
                );
        fprintf( stdout,"Solver starting with relative 2-norm errors of vels, vort, and div at  %.2Le, %.2Le, and %.2Le\n", 
                vel_2_err / vel_2_norm,
                vort_2_err / vort_2_norm,
                div_2_err / div_2_norm
               );

        double min_div_err = 1e10, min_vort_err = 1e10;
        int vort_stagnation_count = 0, div_stagnation_count = 0;
        std::vector<double> dPsidt(Npts_coarse, 0), dPhidt(Npts_coarse, 0);

        bool keep_solving = true;
        unsigned long int iter_cycle = 0;

        // Declare the solver
        LaplaceOperator A( coarsened_grid, cos_ip, cos_im );

        Eigen::BiCGSTAB<LaplaceOperator, Eigen::IdentityPreconditioner> BiCGSTAB_solver;
        BiCGSTAB_solver.compute(A);
        BiCGSTAB_solver.setTolerance(1e-10);
        BiCGSTAB_solver.setMaxIterations(Helmholtz_data.iterations_per_cycle);

        // a small (in # of iters) GMRES solve to clean up grid-scale issues
        Eigen::GMRES<LaplaceOperator, Eigen::IdentityPreconditioner> GMRES_solver;
        GMRES_solver.compute(A);
        GMRES_solver.setTolerance(1e-10);
        GMRES_solver.setMaxIterations(Helmholtz_data.iterations_per_cycle);

        // And make some mappings between std::vector objects and Eigen::VectorXd objects
        std::vector<double> Psi_output_vec(Npts_coarse, 0.), Phi_output_vec(Npts_coarse, 0.);
        Eigen::Map<Eigen::VectorXd> rhs_vort( cgrid_vort.data(), Npts_coarse ),
                                    rhs_div(  cgrid_div.data(),  Npts_coarse ),
                                    Psi_seed( soln_Psi.data(),   Npts_coarse ),
                                    Phi_seed( soln_Phi.data(),   Npts_coarse ),
                                    Psi_output( Psi_output_vec.data(),   Npts_coarse ),
                                    Phi_output( Phi_output_vec.data(),   Npts_coarse );

        while ( keep_solving ) {

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            // primary convergence with BiCGSTAB
            /*
            Psi_output = BiCGSTAB_solver.solveWithGuess( rhs_vort, Psi_seed );
            //fprintf( stdout, "    BiCGSTAB : Psi solver halted after %ld iterations to error of %.5e\n",
            //      BiCGSTAB_solver.iterations(), BiCGSTAB_solver.error() );
            Phi_output = BiCGSTAB_solver.solveWithGuess( rhs_div,  Phi_seed );
            //fprintf( stdout, "    BiCGSTAB : Phi solver halted after %ld iterations to error of %.5e\n",
            //      BiCGSTAB_solver.iterations(), BiCGSTAB_solver.error() );

            std::copy( Psi_output.data(), Psi_output.data()+Psi_output.size(), soln_Psi.begin() );
            std::copy( Phi_output.data(), Phi_output.data()+Phi_output.size(), soln_Phi.begin() );
            */

            // now do a GMRES run to clean up grid scales
            Psi_output = GMRES_solver.solveWithGuess( rhs_vort, Psi_seed );
            //fprintf( stdout, "    GMRES : Psi solver halted after %ld iterations to error of %.5e\n",
            //      GMRES_solver.iterations(), GMRES_solver.error() );
            Phi_output = GMRES_solver.solveWithGuess( rhs_div,  Phi_seed );
            //fprintf( stdout, "    GMRES : Phi solver halted after %ld iterations to error of %.5e\n",
            //      GMRES_solver.iterations(), GMRES_solver.error() );

            std::copy( Psi_output.data(), Psi_output.data()+Psi_output.size(), soln_Psi.begin() );
            std::copy( Phi_output.data(), Phi_output.data()+Phi_output.size(), soln_Phi.begin() );

            // Denoise the signal by with amplitude-modulated demodulation of Nyquist frequency noise
            /*
            const double filter_ell = 3 * constants::R_earth * ( M_PI / Nlat_coarse );
            coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
            coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");
            AM_demodulation( coarsened_grid, filter_ell );
            solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
            solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
            */

            // Filter Psi and Phi to remove grid-scale issues
            /*
            const double filter_ell = 3 * constants::R_earth * ( M_PI / Nlat_coarse );
            coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
            coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");
            filter_vort_and_div( coarsened_grid, filter_ell, true );
            solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
            solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
            */

            /*
            // Remove grid-scale oscillations along lines of constant latitude
            for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                long double alternating_sum_Psi = 0.,
                            alternating_sum_Phi = 0.;

                // Get the alternating mode magnitude
                for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                    size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                    int sign = (Ilon % 2 == 0) ? 1 : -1;
                    alternating_sum_Psi += sign * soln_Psi[Ipt] / Nlon_coarse;
                    alternating_sum_Phi += sign * soln_Phi[Ipt] / Nlon_coarse;
                }

                // Remove the mode
                for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                    size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                    int sign = (Ilon % 2 == 0) ? 1 : -1;
                    soln_Psi[Ipt] -= sign * alternating_sum_Psi;
                    soln_Phi[Ipt] -= sign * alternating_sum_Phi;
                }
            }

            // Remove grid-scale oscillations along lines of constant longitude
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                long double alternating_sum_Psi = 0.,
                            alternating_sum_Phi = 0.;

                // Get the alternating mode magnitude
                for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                    size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                    int sign = (Ilat % 2 == 0) ? 1 : -1;
                    alternating_sum_Psi += sign * soln_Psi[Ipt] / Nlat_coarse;
                    alternating_sum_Phi += sign * soln_Phi[Ipt] / Nlat_coarse;
                }

                // Remove the mode
                for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                    size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                    int sign = (Ilat % 2 == 0) ? 1 : -1;
                    soln_Psi[Ipt] -= sign * alternating_sum_Psi;
                    soln_Phi[Ipt] -= sign * alternating_sum_Phi;
                }
            }
            */

            /*
            // get mean polar values
            double S_pole_Psi = 0, S_pole_Phi = 0, N_pole_Psi = 0, N_pole_Phi = 0.;
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat_coarse-1, Ilon ); 
                size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

                S_pole_Psi += soln_Psi[Ipt_S] / Nlon_coarse;
                S_pole_Phi += soln_Phi[Ipt_S] / Nlon_coarse;

                N_pole_Psi += soln_Psi[Ipt_N] / Nlon_coarse;
                N_pole_Phi += soln_Phi[Ipt_N] / Nlon_coarse;
            }

            // and set the polar rings to those values
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat_coarse-1, Ilon ); 
                size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

                soln_Psi[Ipt_S] = S_pole_Psi;
                soln_Psi[Ipt_N] = N_pole_Psi;

                soln_Phi[Ipt_S] = S_pole_Phi;
                soln_Phi[Ipt_N] = N_pole_Phi;
            }

            fprintf( stdout, "North polar values are: %.4e and %.4e\n", N_pole_Psi, N_pole_Phi );
            fprintf( stdout, "South polar values are: %.4e and %.4e\n", S_pole_Psi, S_pole_Phi );
            */

            /*
            // Apply some damping to meridional grid-scale oscillations
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                    size_t index_im1 = coarsened_grid->local_index( 0, 0, (Ilat == 0) ? 0 : Ilat - 1, Ilon );
                    size_t index_i   = coarsened_grid->local_index( 0, 0, Ilat, Ilon );
                    size_t index_ip1 = coarsened_grid->local_index( 0, 0, (Ilat == Nlat_coarse-1) ? Nlat_coarse-1 : Ilat + 1, Ilon );

                    double Psi_im1 = (Ilat == 0) ? S_pole_Psi : soln_Psi[index_im1];
                    double Psi_i   = soln_Psi[index_i];
                    double Psi_ip1 = (Ilat == Nlat_coarse-1) ? N_pole_Psi : soln_Psi[index_ip1];

                    Psi_output_vec[index_i] = 0.25 * ( Psi_im1 + 2 * Psi_i + Psi_ip1 );

                    double Phi_im1 = (Ilat == 0) ? S_pole_Phi : soln_Phi[index_im1];
                    double Phi_i   = soln_Phi[index_i];
                    double Phi_ip1 = (Ilat == Nlat_coarse-1) ? N_pole_Phi : soln_Phi[index_ip1];

                    Phi_output_vec[index_i] = 0.25 * ( Phi_im1 + 2 * Phi_i + Phi_ip1 );
                }
            }
            std::copy( Psi_output.data(), Psi_output.data()+Psi_output.size(), soln_Psi.begin() );
            std::copy( Phi_output.data(), Phi_output.data()+Phi_output.size(), soln_Phi.begin() );
            */

            /*
            long double accum_Psi = 0, accum_Phi = 0, denom = 0.;
            #pragma omp parallel default(none) private(Ipt) \
            shared( soln_Phi, soln_Psi, coarsened_grid ) \
            reduction( +:accum_Psi,accum_Phi,denom ) \
            firstprivate( Npts_coarse )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                    accum_Psi += soln_Psi[Ipt];
                    accum_Phi += soln_Phi[Ipt];
                    denom += coarsened_grid->areas[Ipt];
                }
            }
            */

            // De-drift Psi, Phi by setting Psi[0] = Phi[0] = 0
            const double Psi_shift = soln_Psi[0],
                         Phi_shift = soln_Phi[0];
            //const double Psi_shift = accum_Psi / denom,
            //             Phi_shift = accum_Phi / denom;
            #pragma omp parallel default(none) private(Ipt) \
            shared( soln_Phi, soln_Psi ) \
            firstprivate( Npts_coarse, Psi_shift, Phi_shift )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                    soln_Phi[Ipt] -= Phi_shift;
                    soln_Psi[Ipt] -= Psi_shift;
                }
            }

            //iter_cycle++;
            iter_cycle += Helmholtz_data.iterations_per_cycle;
            if ( iter_cycle >= Helmholtz_data.iteration_max ) { keep_solving = false; }
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            if ( iter_cycle % Helmholtz_data.iterations_per_cycle == 0 ) {

                // Compute velocity components and vort,div
                get_vels_vort_and_div_Eigen( coarsened_grid, solution_grid, true );

                //double mean_vort, mean_div;
                //demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
                //demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );

                //double mean_Psi, mean_Phi;
                //demean( solution_grid.variables.at("Psi"),  coarsened_grid, mean_Psi );
                //demean( solution_grid.variables.at("Phi" ), coarsened_grid, mean_Phi );

                // Copy the solution over to the other storage set
                coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
                coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");

                // compute errors and check convergence
                long double vel_2_err  = 0, vel_inf_err  = 0,
                            vort_2_err = 0, vort_inf_err = 0,
                            div_2_err  = 0, div_inf_err  = 0;
                long double vel_2_norm  = 0, vel_inf_norm  = 0,
                            vort_2_norm = 0, vort_inf_norm = 0,
                            div_2_norm  = 0, div_inf_norm  = 0;
                long double vort_2_ener = 0, div_2_ener = 0;
                long double Psi_2_norm = 0, Phi_2_norm = 0;

                get_norms_GS(  vel_2_err,  vort_2_err,  div_2_err,
                            vel_2_norm, vort_2_norm, div_2_norm,
                            vel_inf_err,  vort_inf_err,  div_inf_err,
                            vel_inf_norm, vort_inf_norm, div_inf_norm,
                            vort_2_ener, div_2_ener,
                            Psi_2_norm, Phi_2_norm,
                            coarsened_grid
                          );

                // If we've hit tolerance, halt
                if (      (vort_2_err / vort_2_norm < rel_tol) 
                      and (div_2_err  / div_2_norm  < rel_tol) ) {
                    keep_solving = false;
                }

                // If diverged or stagnated, halt
                if ( vort_2_err / vort_2_norm > 500 * min_vort_err ) {
                    //keep_solving = false;
                } else if ( std::fabs( (vort_2_err / vort_2_norm) - min_vort_err ) < std::fmax(1e-4 * min_vort_err, 1e-14) ) {
                    vort_stagnation_count++;
                } else {
                    vort_stagnation_count = 0;
                }
                min_vort_err = std::fmin( min_vort_err, vort_2_err / vort_2_norm );

                // If diverged or stagnated, halt
                if ( div_2_err / div_2_norm > 500 * min_div_err ) {
                    //keep_solving = false;
                } else if ( std::fabs( (div_2_err / div_2_norm) - min_div_err ) < std::fmax(1e-4 * min_div_err, 1e-14) ) {
                    div_stagnation_count++;
                } else {
                    div_stagnation_count = 0;
                }
                min_div_err = std::fmin( min_div_err, div_2_err / div_2_norm );

                // If the solver has stagnated, move on to the next resolution tier
                if ( ( div_stagnation_count > 5 ) and ( vort_stagnation_count > 5 ) ) { keep_solving = false; }

                // Save the convergence records
                Helmholtz_data.vel_2_errors.push_back(  (double) vel_2_err  );
                Helmholtz_data.vort_2_errors.push_back( (double) vort_2_err );
                Helmholtz_data.div_2_errors.push_back(  (double) div_2_err  );

                Helmholtz_data.vel_inf_errors.push_back(  (double) vel_inf_err  );
                Helmholtz_data.vort_inf_errors.push_back( (double) vort_inf_err );
                Helmholtz_data.div_inf_errors.push_back(  (double) div_inf_err  );

                // Save the convergence reference values
                Helmholtz_data.vel_2_norms.push_back(  (double) vel_2_norm  );
                Helmholtz_data.vort_2_norms.push_back( (double) vort_2_norm );
                Helmholtz_data.div_2_norms.push_back(  (double) div_2_norm  );

                Helmholtz_data.vel_inf_norms.push_back(  (double) vel_inf_norm  );
                Helmholtz_data.vort_inf_norms.push_back( (double) vort_inf_norm );
                Helmholtz_data.div_inf_norms.push_back(  (double) div_inf_norm  );


                // If sufficient convergence, stop
                if ( Helmholtz_data.IsConverged() ) { 
                    keep_solving = false;
                    fprintf( stdout,"  Halting at solver cycle %'zu. Solver has converged to desired tolerance.\n", iter_cycle );
                }

                #if DEBUG >= 0
                fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels, vort, and div are %.2e, %.2e%s, and %.2e%s\n", 
                        iter_cycle,
                        Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                        Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                        vort_stagnation_count > 0 ? "(*)" : "",
                        Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() ,
                        div_stagnation_count > 0 ? "(*)" : ""
                       );
                #endif
            }
        }

        fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
               );

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
        coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");

        // Compute velocity components
        get_vels_vort_and_div_Eigen( coarsened_grid, *coarsened_grid, true );

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                      coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                      coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }

        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();
        write_Helmholtz_output( output_filename, Helmholtz_data, *coarsened_grid, solution_grid );

        mean_vort = get_mean( coarsened_grid->variables.at("proj_vort"), coarsened_grid );
        mean_div  = get_mean( coarsened_grid->variables.at("proj_div" ), coarsened_grid );
        /*
        demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );
        */
        fprintf( stdout, "Final mean values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

        // Remove grid-scale oscillations along lines of constant latitude
        // NEED TO ACCOUNT FOR NTIME AND NDEPTH
        for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
            long double alternating_sum_Psi = 0.,
                        alternating_sum_Phi = 0.;

            // Get the alternating mode magnitude
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                int sign = (Ilon % 2 == 0) ? 1 : -1;
                alternating_sum_Psi += sign * soln_Psi[Ipt] / Nlon_coarse;
                alternating_sum_Phi += sign * soln_Phi[Ipt] / Nlon_coarse;
            }

            // Remove the mode
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                size_t Ipt = coarsened_grid->local_index( 0, 0, Ilat, Ilon ); 
                int sign = (Ilon % 2 == 0) ? 1 : -1;
                soln_Psi[Ipt] -= sign * alternating_sum_Psi;
                soln_Phi[Ipt] -= sign * alternating_sum_Phi;
            }
        }
        get_vels_vort_and_div_Eigen( coarsened_grid, *coarsened_grid, false );

        // Need to remove the cos(lat) prefactor on vels and (cos(lat))^2 prefactor on vort/div.
        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                int Itime, Idepth, Ilat, Ilon;
                coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

                double cos_lat = 
                    coarsened_grid->sin_cos_lat_pre_comp ?
                    coarsened_grid->cos_lat[Ilat] :
                    cos( coarsened_grid->latitude[Ilat] );

                coarsened_grid->variables.at("u_lon")[Ipt] =
                    coarsened_grid->variables.at("u_lon")[Ipt] / cos_lat;
                coarsened_grid->variables.at("u_lat")[Ipt] =
                    coarsened_grid->variables.at("u_lat")[Ipt] / cos_lat;

                /*
                coarsened_grid->variables.at("proj_uo_tor")[Ipt] =
                    coarsened_grid->variables.at("proj_uo_tor")[Ipt] / cos_lat;
                coarsened_grid->variables.at("proj_vo_tor")[Ipt] =
                    coarsened_grid->variables.at("proj_vo_tor")[Ipt] / cos_lat;

                coarsened_grid->variables.at("proj_uo_pot")[Ipt] =
                    coarsened_grid->variables.at("proj_uo_pot")[Ipt] / cos_lat;
                coarsened_grid->variables.at("proj_vo_pot")[Ipt] =
                    coarsened_grid->variables.at("proj_vo_pot")[Ipt] / cos_lat;
                */

                const double cos_sqr = cos_lat * cos_lat;

                coarsened_grid->variables.at("vort")[Ipt] = 
                    coarsened_grid->variables.at("vort")[Ipt] / cos_sqr;
                coarsened_grid->variables.at("div")[Ipt] = 
                    coarsened_grid->variables.at("div")[Ipt] / cos_sqr;

                /*
                coarsened_grid->variables.at("proj_vort")[Ipt] = 
                    coarsened_grid->variables.at("proj_vort")[Ipt] / cos_sqr;
                coarsened_grid->variables.at("proj_div")[Ipt] = 
                    coarsened_grid->variables.at("proj_div")[Ipt] / cos_sqr;
                */
            }
        }

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                      coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                      coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "ProcessingSolution"); }
        
        // Write the output
        std::ostringstream name_stream_v2;
        name_stream_v2 << "projection_R" << (refine_level) << "_NoCos.nc";
        const std::string output_filename_v2 = name_stream_v2.str();
        write_Helmholtz_output( output_filename_v2, Helmholtz_data, *coarsened_grid, solution_grid );
        //std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );


        if ( refine_level > 0 ) {
            // Smooth before passing up to the next grid level
            //   otherwise the interpolator causes all kinds of noisy mess
            const double filter_ell = 4 * constants::R_earth * ( M_PI / Nlat_coarse );
            coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
            coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");
            filter_vort_and_div( coarsened_grid, filter_ell, true );
            solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
            solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
        }
        
        //
        //// Also write a version that uses a minor smoothing filter at the end
        //

        /*
        get_vels_vort_and_div_Eigen( coarsened_grid, *coarsened_grid );

        demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );
        fprintf( stdout, "Final mean values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                      coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                      coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }

        // Write the output
        std::ostringstream name_stream_filtered;
        name_stream_filtered << "projection_FILT_R" << (refine_level) << ".nc";
        const std::string output_filename_filtered = name_stream_filtered.str();
        write_Helmholtz_output( output_filename_filtered, Helmholtz_data, *coarsened_grid, solution_grid );
        */

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
