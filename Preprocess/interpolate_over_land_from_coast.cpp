#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/interpolation.h"
#include <vector>
#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"

void interpolate_over_land_from_coast(
        std::vector<double> & interp_field,
        const std::vector<double> & field,
        const int nlayers,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<bool>   & mask,
        const std::vector<int>    & myCounts,
        const MPI_Comm comm
        ){

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const bool cast_to_sphere = false;

    interp_field.resize(field.size());

    //
    // Step 1: RBF model creation
    //
    alglib::rbfmodel model;
    if (cast_to_sphere) { 
        alglib::rbfcreate(3, 1, model);
    } else {
        alglib::rbfcreate(2, 1, model);
    }

    //
    // Step 2: we add dataset of known points
    //

    std::vector<double> xyzf_vec;

    int cntr, num_coast=0;
    double R, val;
    int Ilat, Ilon;
    size_t index = 0, num_interp_points = 0;
    #if DEBUG >= 0
    const double perc_base = 5;
    double perc;
    int tid, perc_count = 0;
    #endif

    double x, y, z, lat, lon;

    double rbase;
    // largest interpolation scale is scale * dx
    //   smallest is largest / 2^nlayers
    // Current setting sweeps dx -> 2^nlayers * dx
    const double scale = pow(2., (double)nlayers);

    std::vector<double> lon_coast;
    std::vector<double> lat_coast;
    std::vector<double> field_coast;

    if (cast_to_sphere) { 
        rbase = scale * distance(0, 0, 0, latitude.at(1) - latitude.at(0));
        #if DEBUG >= 2
        fprintf(stdout, " Using rbase = %.3g km and nlayers = %d\n", rbase, nlayers);
        #endif
    } else {
        rbase = scale * fabs(latitude.at(1) - latitude.at(0));
        #if DEBUG >= 2
        fprintf(stdout, " Using rbase = %.3g degrees and nlayers = %d\n", rbase, nlayers);
        #endif
    }

    for (int Itime = 0; Itime < Ntime; Itime++) {
    #if DEBUG >= 1
    fprintf(stdout, "  Itime = %d of %d\n", Itime + 1, Ntime);
    #endif
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #if DEBUG >= 1
            fprintf(stdout, "    Idepth = %d of %d\n", Idepth + 1, Ndepth);
            #endif
    
            get_coast(lon_coast, lat_coast, field_coast,
                      longitude, latitude, field, mask,
                      Itime, Idepth,
                      Ntime, Ndepth, Nlat, Nlon);

            num_coast = lon_coast.size();

            if (cast_to_sphere) { 
                xyzf_vec.resize(4 * num_coast);
            } else {
                xyzf_vec.resize(3 * num_coast);
            }

            R = constants::R_earth - depth.at(Idepth);
            #if DEBUG >= 1
            fprintf(stdout, "      Adding seed data to the interpolator object...");
            fflush(stdout);
            perc = perc_base;
            #endif
            cntr = 0;

            for (int Icoast = 0; Icoast < num_coast; Icoast++) {
                lat = lat_coast.at(Icoast);
                lon = lon_coast.at(Icoast);
                val = field_coast.at(Icoast);

                if (cast_to_sphere) { 
                    xyzf_vec.at(4*cntr + 0) = R * cos(lat) * cos(lon);
                    xyzf_vec.at(4*cntr + 1) = R * cos(lat) * sin(lon);
                    xyzf_vec.at(4*cntr + 2) = R * sin(lat);
                    xyzf_vec.at(4*cntr + 3) = val;
                } else {
                    xyzf_vec.at(3*cntr + 0) = lon;
                    xyzf_vec.at(3*cntr + 1) = lat;
                    xyzf_vec.at(3*cntr + 2) = val;
                }

                cntr++;

            }
            #if DEBUG >= 1
            fprintf(stdout, " done\n");
            fflush(stdout);
            #endif

            #if DEBUG >= 1
            fprintf(stdout, "      Setting points for the model...");
            fflush(stdout);
            #endif
            alglib::real_2d_array xyzf;
            if (cast_to_sphere) { 
                xyzf.setlength(num_coast, 4);
                xyzf.setcontent(num_coast, 4, &xyzf_vec[0]);
            } else {
                xyzf.setlength(num_coast, 3);
                xyzf.setcontent(num_coast, 3, &xyzf_vec[0]);
            }

            #if DEBUG >= 1
            fprintf(stdout, " ...");
            fflush(stdout);
            #endif
        
            alglib::rbfsetpoints(model, xyzf);
            #if DEBUG >= 1
            fprintf(stdout, " done\n");
            fflush(stdout);
            #endif
        
            //
            // Step 3: rebuild model
            //
            // After we've configured model, we should rebuild it -
            // it will change coefficients stored internally in the
            // rbfmodel structure.
            //
            // We use hierarchical RBF algorithm with following parameters:
            // * RBase - set to 1.0
            // * NLayers - three layers are used
            // * LambdaReg - is set to zero value, no smoothing is required
            //
            #if DEBUG >= 1
            fprintf(stdout, "      Building the interpolator object.\n");
            fflush(stdout);
            #endif
            alglib::rbfreport rep;
            alglib::rbfsetalgohierarchical(model, rbase, nlayers, 0.0);
            alglib::rbfbuildmodel(model, rep);
        
            //
            // Step 4: model was built
            //         now we use the interpolator at each land point
            //
            #if DEBUG >= 1
            fprintf(stdout, "      Filling land points with the interpolator object.\n");
            perc = perc_base;
            #endif

            for ( Ilat = 0; Ilat < Nlat; Ilat++) {
                for ( Ilon = 0; Ilon < Nlon; Ilon++) {
                    lat = latitude.at(Ilat);
                    lon = longitude.at(Ilon);

                    index = Index(Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon);

                    if (not(mask.at(index))) {
                        // If we're on a land cell, then use the interpolator
                        if (cast_to_sphere) { 
                            x = R * cos(lat) * cos(lon);
                            y = R * cos(lat) * sin(lon);
                            z = R * sin(lat);
                            interp_field.at(index) = alglib::rbfcalc3(model, x, y, z);
                        } else {
                            interp_field.at(index) = alglib::rbfcalc2(model, lon, lat);
                        }

                        num_interp_points++;
                    } else {
                        // Otherwise use the true data
                        interp_field.at(index) = field.at(index);
                    }

                    #if DEBUG >= 0
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every 5 percent, print a dot
                        if ( ((double)(Ilat*Nlon + Ilon) / (Nlat*Nlon)) * 100 >= perc ) {
                            if (perc == perc_base) { fprintf(stdout, "          "); }
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif
                }
            }
        }
        #if DEBUG >= 0
        fprintf(stdout, "\n");
        #endif
    }

    #if DEBUG >= 1
    fprintf(stdout, "  %'d seed points used to determine %'zu interpolated points.\n", num_coast, num_interp_points);
    #endif

}

