#include <algorithm>
#include <mpi.h>

#include "../../constants.hpp"
#include "../../functions.hpp"

// This parser tool was provided by StackOverflow user 'iain' at the following post.
//    mild modifications have been made to adapt the code to our purposes.
// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i) {
        this->tokens.push_back(std::string(argv[i]));
    }
}

const std::string InputParser::getCmdOption(
        const std::string &option,
        const std::string &default_value
        ) const{

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    std::vector<std::string>::const_iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        #if DEBUG >= 1
        if (wRank == 0) {
            fprintf(stdout, " Commandline flag \"%s\" got value \"%s\"\n", option.c_str(), itr->c_str());
        }
        #endif
        return *itr;
    }
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, " Commandline flag \"%s\" received no value - will use default \"%s\"\n", 
                option.c_str(), default_value.c_str());
    }
    #endif
    return default_value;
}

bool InputParser::cmdOptionExists(const std::string &option) const{
    return ( std::find(this->tokens.begin(), this->tokens.end(), option)  !=  this->tokens.end() );
}
