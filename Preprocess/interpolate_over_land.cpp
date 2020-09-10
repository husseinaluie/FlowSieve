#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/interpolation.h"
#include <vector>
#include "../constants.hpp"
#include "../functions.hpp"

void interpolate_over_land(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<bool>   &mask)
{
    // Count how mayn water points there are (since these are the data points
    //   that are known / valid)
    int num_water_pts = 0;
    for (size_t II = 0; II < mask.size(); II++) { 
        if (mask.at(II)) { num_water_pts++; }
    }

    int Ntime  = (int)time.size();
    int Ndepth = (int)depth.size();
    int Nlat   = (int)latitude.size();
    int Nlon   = (int)longitude.size();

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
    if (cast_to_sphere) { 
        xyzf_vec.resize(4*num_water_pts);
    } else {
        xyzf_vec.resize(3*num_water_pts);
    }

    int cntr;
    double R;
    int index, mask_index, num_seed_points = 0, num_interp_points = 0;
    #if DEBUG >= 1
    double perc_base = 10;
    double perc;
    #endif

    double x, y, z, lat, lon;
    const double D2R = M_PI / 180;

    double rbase;
    const int nlayers = 4;
    if (cast_to_sphere) { 
        rbase = 10 * distance(0, 0, 0, latitude.at(1) - latitude.at(0));
        #if DEBUG >= 2
        fprintf(stdout, " Using rbase = %.3gkm and nlayers = %d\n", rbase, nlayers);
        #endif
    } else {
        rbase = 10 * fabs(latitude.at(1) - latitude.at(0));
        #if DEBUG >= 2
        fprintf(stdout, " Using rbase = %.3grad and nlayers = %d\n", rbase, nlayers);
        #endif
    }

    //for (int Itime = 0; Itime < 1; Itime++) {
    for (int Itime = 0; Itime < Ntime; Itime++) {
    #if DEBUG >= 1
    fprintf(stdout, "  Itime = %d of %d\n", Itime + 1, Ntime);
    #endif
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #if DEBUG >= 1
            fprintf(stdout, "    Idepth = %d of %d\n", Idepth + 1, Ndepth);
            #endif
            R = constants::R_earth - depth.at(Idepth);
            #if DEBUG >= 1
            fprintf(stdout, "      Adding seed data to the interpolator object.\n");
            perc = perc_base;
            #endif
            cntr = 0;
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                lat = latitude.at(Ilat);
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                    lon = longitude.at(Ilon);

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask.at(mask_index)) {
                        if (cast_to_sphere) { 
                            xyzf_vec.at(4*cntr + 0) = R * cos(lat * D2R) * cos(lon * D2R);
                            xyzf_vec.at(4*cntr + 1) = R * cos(lat * D2R) * sin(lon * D2R);
                            xyzf_vec.at(4*cntr + 2) = R * sin(lat * D2R);
                            xyzf_vec.at(4*cntr + 3) = field.at(index);
                        } else {
                            xyzf_vec.at(3*cntr + 0) = lon;
                            xyzf_vec.at(3*cntr + 1) = lat;
                            xyzf_vec.at(3*cntr + 2) = field.at(index);
                        }

                        num_seed_points++;
                        cntr++;
                    }

                    #if DEBUG >= 1
                    // Every 5 percent, print a dot
                    if ( ((double)(mask_index) / (Nlat*Nlon)) * 100 >= perc ) {
                        if (perc == perc_base) { fprintf(stdout, "          "); }
                        fprintf(stdout, ".");
                        fflush(stdout);
                        perc += perc_base;
                    }
                    #endif
                }
            }
            #if DEBUG >= 1
            fprintf(stdout, "\n");
            #endif

            alglib::real_2d_array xyzf;
            if (cast_to_sphere) { 
                xyzf.setlength(num_water_pts, 4);
                xyzf.setcontent(num_water_pts, 4, &xyzf_vec[0]);
            } else {
                xyzf.setlength(num_water_pts, 3);
                xyzf.setcontent(num_water_pts, 3, &xyzf_vec[0]);
            }
        
            alglib::rbfsetpoints(model, xyzf);
        
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
            #endif
            alglib::rbfreport rep;
            fflush(stdout);
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

            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                lat = latitude.at(Ilat);
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                    lon = longitude.at(Ilon);

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (not(mask.at(mask_index))) {
                        // If we're on a land cell, then use the interpolator
                        if (cast_to_sphere) { 
                            x = R * cos(lat * D2R) * cos(lon * D2R);
                            y = R * cos(lat * D2R) * sin(lon * D2R);
                            z = R * sin(lat * D2R);
                            interp_field.at(index) = alglib::rbfcalc3(model, x, y, z);
                        } else {
                            interp_field.at(index) = alglib::rbfcalc2(model, lon, lat);
                        }

                        num_interp_points++;
                    } else {
                        // Otherwise use the true data
                        interp_field.at(index) = field.at(index);
                    }

                    #if DEBUG >= 1
                    // Every 5 percent, print a dot
                    if ( ((double)(mask_index) / (Nlat*Nlon)) * 100 >= perc ) {
                        if (perc == perc_base) { fprintf(stdout, "          "); }
                        fprintf(stdout, ".");
                        fflush(stdout);
                        perc += perc_base;
                    }
                    #endif

                }
            }
        }
        #if DEBUG >= 1
        fprintf(stdout, "\n");
        #endif
    }

    #if DEBUG >= 1
    fprintf(stdout, "  %d seed points used to determine %d interpolated points.\n", num_seed_points, num_interp_points);
    #endif

}

