
#include <algorithm>
#include <math.h>
#include "../functions.hpp"
#include <vector>
#include "../preprocess.hpp"


void get_coast(
        std::vector<double> & lon_coast,
        std::vector<double> & lat_coast,
        std::vector<double> & field_coast,
        const std::vector<double> & lon_full,
        const std::vector<double> & lat_full,
        const std::vector<double> & field_full,
        const std::vector<bool>   & mask,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon) {

    int num_coast = 0;
    bool is_coast;

    const int coast_width = 20; // width in grid points

    int curr_lon, curr_lat, curr_index, index;

    // First, we do a large loop to simply
    //   count how many coast points we have.
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {

            is_coast = false;

            index = Index(Itime, Idepth, Ilat, Ilon,
                          Ntime, Ndepth, Nlat, Nlon);

            if (mask.at(index)) {
                for (int II = 0; II < coast_width+1; II++) {

                    // Look northward
                    curr_lat = std::min(Nlat-1, Ilat+II);
                    for (int JJ = 0; JJ < coast_width+1-II; JJ++) {
                        
                        // Look westward
                        curr_lon = Ilon + JJ;
                        if (curr_lon >= Nlon) { curr_lon -= Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }

                        // Look eastward
                        curr_lon = Ilon - JJ;
                        if (curr_lon < 0) { curr_lon += Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }
                    }  

                    // Look southward
                    curr_lat = std::max(0, Ilat-II);
                    for (int JJ = 0; JJ < coast_width+1-II; JJ++) {
                        
                        // Look westward
                        curr_lon = Ilon + JJ;
                        if (curr_lon >= Nlon) { curr_lon -= Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }

                        // Look eastward
                        curr_lon = Ilon - JJ;
                        if (curr_lon < 0) { curr_lon += Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }
                    }  
                }
            }
            if (is_coast) { num_coast++; }
        }
    }

    fprintf(stdout, "      ... found %d coast points (%.4g percent)\n", num_coast, (100.*num_coast) / (Nlat*Nlon));

    // Now that we know how many points there are, we can resize the appropriate arrays
    lon_coast.resize(num_coast);
    lat_coast.resize(num_coast);
    field_coast.resize(num_coast);

    // Now that we know how many points there are,
    //   we can actually store them.
    int cnt = 0;
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {

            index = Index(Itime, Idepth, Ilat, Ilon,
                          Ntime, Ndepth, Nlat, Nlon);

            is_coast = false;

            if (mask.at(index)) {
                for (int II = 0; II < coast_width+1; II++) {

                    // Look northward
                    curr_lat = std::min(Nlat-1, Ilat+II);
                    for (int JJ = 0; JJ < coast_width+1-II; JJ++) {
                        
                        // Look westward
                        curr_lon = Ilon + JJ;
                        if (curr_lon >= Nlon) { curr_lon -= Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }

                        // Look eastward
                        curr_lon = Ilon - JJ;
                        if (curr_lon < 0) { curr_lon += Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }
                    }  

                    // Look southward
                    curr_lat = std::max(0, Ilat-II);
                    for (int JJ = 0; JJ < coast_width+1-II; JJ++) {
                        
                        // Look westward
                        curr_lon = Ilon + JJ;
                        if (curr_lon >= Nlon) { curr_lon -= Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }

                        // Look eastward
                        curr_lon = Ilon - JJ;
                        if (curr_lon < 0) { curr_lon += Nlon; } 
                        curr_index = Index(Itime, Idepth, curr_lat, curr_lon,
                                           Ntime, Ndepth, Nlat,     Nlon);
                        if (not(mask.at(curr_index))) { is_coast = true; }
                    }  
                }
            }
            if (is_coast) {
                lon_coast.at(cnt)   = lon_full.at(Ilon);
                lat_coast.at(cnt)   = lat_full.at(Ilat);
                field_coast.at(cnt) = field_full.at(index);
                cnt++;
            }
        }
    }

}
