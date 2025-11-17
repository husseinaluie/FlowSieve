#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

void HelmholtzDataClass::VerifyRowNorms(
        ) {

    // Compute row and column norms
    std::vector<double> col_norms(2*Ncol,0), row_norms( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow,0);
    size_t col = 0, row = 0;
    for ( int k = 0; k < LHS.outerSize(); ++k ) {
        for ( Eigen::SparseMatrix<double>::InnerIterator it(LHS,k); it; ++it ) {
            double val = it.value();

            row = it.row();   // row index
            row_norms[row] += pow(val,2) / (2*Ncol);

            col = it.col();   // col index (here it is equal to k)
            col_norms[col] += pow(val,2) / ( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow);
        }
    }
    double min_row_norm = 1e10, max_row_norm = 0,
           min_col_norm = 1e10, max_col_norm = 0;
    size_t num_zero_rows = 0,
           num_zero_cols = 0;
    for ( row = 0; row < ( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow); row++ ) {
        row_norms[row] = sqrt(row_norms[row]);
        if ( row_norms[row] == 0 ) {
            num_zero_rows++;
            #if DEBUG >= 1
            fprintf( stdout, "  Row %zu has zero norm. Setting accompanying RHS to zero.\n", row );
            #endif
            RHS[row] = 0;
        } else {
            min_row_norm = std::fmin( min_row_norm, row_norms[row] );
            max_row_norm = std::fmax( max_row_norm, row_norms[row] );
        }
    }
    for ( col = 0; col < 2*Ncol; col++ ) {
        col_norms[col] = sqrt(col_norms[col]);
        if ( col_norms[col] == 0 ) {
            num_zero_cols++;
        } else {
            min_col_norm = std::fmin( min_col_norm, col_norms[col] );
            max_col_norm = std::fmax( max_col_norm, col_norms[col] );
        }
    }

    #if DEBUG >= 0
    fprintf(stdout, "Column norms were bounded between %e and %e.\n", min_col_norm, max_col_norm);
    fprintf(stdout, "Row norms were bounded between %e and %e.\n", min_row_norm, max_row_norm);
    if ( num_zero_rows > 0 ) { fprintf(stdout, "  %'zu of the rows were zero\n", num_zero_rows); }
    if ( num_zero_cols > 0 ) { fprintf(stdout, "  %'zu of the cols were zero\n", num_zero_cols); }
    #endif

}
