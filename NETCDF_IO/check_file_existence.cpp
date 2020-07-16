#include <string>
#include <stdio.h>

bool check_file_existence (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        fprintf(stderr, " %s does not exist.\n", name.c_str() );
        return false;
    }   
}

