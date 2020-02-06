#ifndef FUNCTIONS_SW_HPP
#define FUNCTIONS_SW_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>

void filtering_sw_2L(
        const std::vector<double> & full_u,
        const std::vector<double> & full_v,
        const std::vector<double> & full_h,
        const std::vector<double> & rhos,
        const double nu,
        const std::vector<double> & scales,
        const std::vector<double> & dAreas,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD);

void filtering_sw(
        const std::vector<double> & u, 
        const std::vector<double> & v, 
        const std::vector<double> & h,
        const std::vector<double> & scales, 
        const std::vector<double> & dAreas, 
        const std::vector<double> & time, 
        const std::vector<double> & depth,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD);

void Compute_pressure_2L(
        std::vector<double> & pressure, 
        const std::vector<double> & h, 
        const std::vector<double> & rho,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_gradient(
        std::vector<double> & field_x, 
        std::vector<double> & field_y, 
        const std::vector<double> & field,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon,
        const std::vector<double> & mask);

void Compute_Laplacians(
        std::vector<double> & h_lap_u, 
        std::vector<double> & h_lap_v, 
        std::vector<double> &   lap_h, 
        const std::vector<double> & full_u, 
        const std::vector<double> & full_v, 
        const std::vector<double> & full_h, 
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon,
        const std::vector<double> & mask);

void Compute_Viscous_Loss_KE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & lap_u_tilde, 
        const std::vector<double> & lap_v_tilde, 
        const double nu,
        const std::vector<double> & rho, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Viscous_Loss_PE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & lap_h_bar, 
        const double nu,
        const std::vector<double> & rho, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Baroclinic_Transfer_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & p_bar,
        const double alpha,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Traditional_Transfer_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & PE,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Pi_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & uu_tilde,
        const std::vector<double> & uv_tilde,
        const std::vector<double> & vv_tilde,
        const std::vector<double> & rho,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Gamma_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar,
        const std::vector<double> & v_bar,
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & rho,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Transport_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u,
        const std::vector<double> & v,
        const std::vector<double> & field,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_Alt_PE_Transport(
        std::vector<double> & out_array, 
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & h_bar,
        const std::vector<double> & rho,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_KE_Transport_smallscales_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & uu_tilde,
        const std::vector<double> & uv_tilde,
        const std::vector<double> & vv_tilde,
        const std::vector<double> & rho,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_misc_conversion_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos,
        const double alpha,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_misc1_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_misc2_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & tau_hpx, 
        const std::vector<double> & tau_hpy, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_misc3_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos,
        const double alpha,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon);

void Compute_full_bc_term_KE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & p_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        );

void Compute_full_bc_term_PE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        );

void Compute_full_bc_term_PE_parts(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        );

#endif
