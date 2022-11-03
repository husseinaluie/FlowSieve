#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <mpi.h>
#include <cassert>

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
        const std::string &default_value,
        const bool help,
        const std::string &description
        ) const{

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    if (help) {
        const int paddedGoal   = 30;
        const int paddedLength = option.size() >= paddedGoal ? 0 : paddedGoal - option.size();
        std::string opt_name = option;
        std::string padded_description = description;
        std::string dot_padding;

        // pad spaces after the option name to make things format nicely
        //opt_name.insert(opt_name.end(), paddedLength, '.');
        dot_padding.insert(dot_padding.end(), paddedLength, '.');

        // Also add padding before any newlines and the beginning of the description
        std::vector<int> newline_locs( 1, 0 ); 

        // find the newlines
        for( int ii = 0; ii < description.size(); ++ii ) {
            if (description[ii] == '\n') { newline_locs.push_back( ii+1 ); }
        }

        // put the padding in
        for( int ii = newline_locs.size()-1; ii >= 0; --ii ) {
            padded_description.insert( padded_description.begin() + newline_locs[ii], paddedGoal+4, ' ');
        }

        // print out the infor
        fprintf(stdout, "   \033[1m%s\033[0m%s [ %s ]\n", opt_name.c_str(), dot_padding.c_str(), default_value.c_str());
        if ( not( description == "" ) ) {
            fprintf(stdout, "%s\n", padded_description.c_str() );
        }
        return default_value;
    } else {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);

        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            #if DEBUG >= 0
            if (wRank == 0) {
                fprintf(stdout, " Commandline flag \"%s\" got value \"%s\"\n", option.c_str(), itr->c_str());
            }
            #endif
            return *itr;
        }
        #if DEBUG >= 0
        if (wRank == 0) {
            fprintf(stdout, " Commandline flag \"%s\" received no value - will use default \"%s\"\n", 
                    option.c_str(), default_value.c_str());
        }
        #endif
        return default_value;
    }
}

bool InputParser::cmdOptionExists(const std::string &option) const{
    return ( std::find(this->tokens.begin(), this->tokens.end(), option)  !=  this->tokens.end() );
}

void InputParser::getFilterScales( 
        std::vector<double> &filter_scales, 
        const std::string &argname,
        const bool help
        ) const{

    int wRank=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );

    //using namespace std;

    const std::string string_of_scales = getCmdOption( argname, " ", help,
            "Space-separated list of filter scales in metres.\ne.g. '100.e3 250.e3 1000.e3'");
    if (help) {
        return;
    }
    assert( string_of_scales.size() > 0 );

    std::istringstream iss( string_of_scales );

    // Split up the list of inputs based on white space into separate strings
    std::vector< std::string > scales_as_strings;
    copy( std::istream_iterator< std::string >(iss), 
          std::istream_iterator< std::string >(), 
          std::back_inserter( scales_as_strings ));

    // Convert the strings into doubles
    #if DEBUG >= 0
    if (wRank == 0) { fprintf(stdout, "Filter scales (%zu) are: ", scales_as_strings.size()); }
    #endif
    filter_scales.resize( scales_as_strings.size() );
    for (size_t ii = 0; ii < scales_as_strings.size(); ++ii) {
        filter_scales.at(ii) = strtod( scales_as_strings.at(ii).c_str(), NULL );
        if (filter_scales.at(ii) <= 0) {
            fprintf(stderr, 
                    "\nReceived bad filter scale (%s). Input must be of form '1.3e4 678e6' etc (i.e. space-separated list of numbers).\n", 
                    scales_as_strings.at(ii).c_str());
            assert(false);
        }
        #if DEBUG >= 0
        if (wRank == 0) { 
            double curr_scale = filter_scales.at(ii);
            if ( curr_scale >= 1000. ) {
                fprintf(stdout, " %'gkm, ", curr_scale / 1e3 ); 
            } else {
                fprintf(stdout, " %'gm, ",  curr_scale ); 
            }
        }
        #endif
    }
    #if DEBUG >= 0
    if (wRank == 0) { fprintf(stdout, "\n\n"); }
    #endif

}

void InputParser::getListofStrings( 
        std::vector<std::string> &list_of_strings, 
        const std::string &argname,
        const bool help,
        const std::string &description
        ) const{

    int wRank=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );

    //using namespace std;
    const std::string raw_input_string = getCmdOption( argname, "", help, description );
    if (help) { return; }
    assert( raw_input_string.size() > 0 );

    std::istringstream iss( raw_input_string );

    // Split up the list of inputs based on white space into separate strings
    copy( std::istream_iterator< std::string >(iss), 
          std::istream_iterator< std::string >(), 
          std::back_inserter( list_of_strings ));

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "String arguments for %s are: ", argname.c_str()); }
    for ( size_t II = 0; II < list_of_strings.size(); II++ ) {
        if (wRank == 0) { fprintf(stdout, "  %s", list_of_strings.at(II).c_str()); }
    }
    if (wRank == 0) { fprintf(stdout, "\n\n"); }
    #endif

}
