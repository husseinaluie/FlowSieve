//
// This functionality was copied from
//   https://stackoverflow.com/questions/3613284/c-stdstring-to-boolean
//

#include <sstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cctype>

bool string_to_bool(std::string str) {

    bool b;

    // Convert the input string to lower case, for user convenience
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

    std::istringstream is(str);

    is >> std::boolalpha >> b;

    return b;
}
