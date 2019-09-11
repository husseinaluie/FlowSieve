#include <math.h>
#include <string>
#include <vector>
#include "../postprocess.hpp"

const std::vector< bool(*)(double, double) > RegionTest::all_regions = {
    RegionTest::Global, 
    RegionTest::GulfofMexico,
    RegionTest::GulfStream, 
    RegionTest::Equator, 
    RegionTest::NorthofEquator, 
    RegionTest::SouthofEquator,
    RegionTest::Kuroshio,
    RegionTest::AntarcticCircumpolar
};

// MUST BE IN THE SAME ORDER AS ABOVE!!!
const std::vector< std::string > RegionTest::region_names = {
    "Global", 
    "Gulf_of_Mexico",
    "Gulf_Stream", 
    "Equator",
    "North_of_Equator", 
    "South_of_Equator",
    "Kuroshio",
    "ACC"
};


// Global: everywhere, always returns true
bool RegionTest::Global(double latitude, double longitude) {
    return true;
}


// Gulf of Mexico
bool RegionTest::GulfofMexico(double latitude, double longitude){
    double lon = longitude * D2R;
    double lat = latitude * D2R;
    if (     (lon > -100)
         and (lon < -83) 
         and (lat > 15) 
         and (lat < 32) 
         and (lat - 21 > -0.68 * (lon + 100))
       ) { return true;  }
    else { return false; }
}


// Gulf Stream
bool RegionTest::GulfStream(double latitude, double longitude){
    // Bounded on east by -35
    // Bounded on west by -80
    // Bounded above and below by abs( lat - 0.4 * lon - 62) > 5
    double lon = longitude * D2R;
    double lat = latitude * D2R;
    if (     ( lon <= -35 )
         and ( lon >= -80.75 )
         and ( fabs( lat - 0.4 * lon - 62 ) <= 6 ) 
       ) { return true;  }
    else { return false; }
}


// Equator
bool RegionTest::Equator(double latitude, double longitude){

    if ( fabs(latitude) < 13 * D2R ) { return true;  }
    else                             { return false; }

}


// Northern Hemisphere (less equator)
bool RegionTest::NorthofEquator(double latitude, double longitude){

    if ( latitude >= 13 * D2R ) { return true;  }
    else                        { return false; }
}


// Souther Hemisphere (less equator)
bool RegionTest::SouthofEquator(double latitude, double longitude){

    if ( latitude <= -14 * D2R ) { return true;  }
    else                         { return false; }
}


// Kuroshio
bool RegionTest::Kuroshio(double latitude, double longitude) {
    double lon = longitude * D2R;
    double lat = latitude * D2R;
    if (     (lon > 120)
         and (lon < 170)
         and (lat > 17)
         and (lat < 45)
         and (1 - (lat - 30 > (lon - 120) * ( 15/20 )) )
         and (1 - (lat - 17 < (lon - 120) * (  8/20 )) * (lon <= 140) )
         and (1 - (lat < 25) * (lon >= 140))
       ) { return true;  }
    else { return false; }
}


// Antarctic Circumpolar Current
bool RegionTest::AntarcticCircumpolar(double latitude, double longitude){
    double lon = longitude * D2R;
    double lat = latitude * D2R;
    if (     (lat < -33) 
         and (lat > -70) 
         and (1 - (lat + 45 > (lon + 180) * (- 5/108)) * (lon < -72) ) 
         and (1 - (lat + 33 > (lon -  20) * (-12/160)) * (lon >  20) )
       ) { return true;  }
    else { return false; }
}
