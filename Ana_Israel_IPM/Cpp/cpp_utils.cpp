#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// convert a matrix to vector. Slightly faster than just using as.vector for large
// objects

// [[Rcpp::export]]

NumericVector mat_to_df(NumericMatrix mat) {
  
  int rows = mat.nrow();
  int cols = mat.ncol();
  
  int out_size = rows * cols;
  
  NumericVector out(out_size);
  
  int it = 0;
  
  for(int i = 0; i < rows; i++) {
    
    for(int j = 0; j < cols; j++) {
      
      out(it) = mat(i, j);
      it = it + 1;
      
    }
    
  }
  
  return out;
  
}

