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
//   name length should not exceed 20 characters
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
    const double lon = longitude * R2D;
    const double lat = latitude  * R2D;
    if (     (lon > -100.)
         and (lon < -83.) 
         and (lat > 15.) 
         and (lat < 32.) 
         and (lat - 21. > -0.68 * (lon + 100.))
       ) { return true;  }
    else { return false; }
}


// Gulf Stream
bool RegionTest::GulfStream(double latitude, double longitude){
    const double lon = longitude * R2D;
    const double lat = latitude  * R2D;
    if (     ( lon <= -35. )
         and ( lon >= -80.75 )
         and ( fabs( lat - 0.4 * lon - 62. ) <= 6. ) 
       ) { return true;  }
    else { return false; }
}


// Equator
bool RegionTest::Equator(double latitude, double longitude){
    const double lat = latitude * R2D;

    if ( fabs(lat) < 15. ) { return true;  }
    else                   { return false; }

}


// Northern Hemisphere (less equator)
bool RegionTest::NorthofEquator(double latitude, double longitude){
    const double lat = latitude * R2D;

    if ( lat >= 15. ) { return true;  }
    else              { return false; }
}


// Souther Hemisphere (less equator)
bool RegionTest::SouthofEquator(double latitude, double longitude){
    const double lat = latitude * R2D;

    if ( lat <= -15. ) { return true;  }
    else               { return false; }
}


// Kuroshio
bool RegionTest::Kuroshio(double latitude, double longitude) {
    const double lon = longitude * R2D;
    const double lat = latitude  * R2D;
    if (     (lon > 120.)
         and (lon < 170.)
         and (lat > 17.)
         and (lat < 45.)
         and not( (  (lat - 30.) > ( (lon - 120.) * ( 15./20. ) ) ) )
         and not( (  (lat - 17.) < ( (lon - 120.) * (  8./20. ) ) ) and (lon <= 140.) )
         and not( (lat < 25.) and (lon >= 140.))
       ) { return true;  }
    else { return false; }
}


// Antarctic Circumpolar Current
bool RegionTest::AntarcticCircumpolar(double latitude, double longitude){
    const double lon = longitude * R2D;
    const double lat = latitude  * R2D;
    if (     (lat < -33.) 
         and (lat > -70.) 
         and not( ( (lat + 45.) > ( (lon + 180.) * (- 5./108.) ) ) and (lon < -72.) ) 
         and not( ( (lat + 33.) > ( (lon -  20.) * (-12./160.) ) ) and (lon >  20.) )
       ) { return true;  }
    else { return false; }
}

