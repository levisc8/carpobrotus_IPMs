#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// determine who lived and died based on presence-absence of id value at T+1

// [[Rcpp::export]]

Rcpp::IntegerVector update_survival(Rcpp::CharacterVector pop2,
                                    Rcpp::IntegerVector survival) {
  
  int N = pop2.size();
  Rcpp::IntegerVector out(N);
  
  for(int i = 0; i < N; i++) {
    if(survival[i] == 1) { 
      
      out[i] = 1L;
      
    } else if(!Rcpp::CharacterVector::is_na(pop2[i]) && 
      Rcpp::IntegerVector::is_na(survival[i])) {
      
      out[i] = 1L;
      
    } else {
      
      out[i] = 0L;
      
    }
  }
  
  return out;
}

// merge ramet IDs together if one has subsumed the other

// [[Rcpp::export]]
Rcpp::IntegerVector merge_ramets(Rcpp::IntegerVector id,
                                 Rcpp::IntegerVector absorbed,
                                 Rcpp::IntegerVector to_sub) {
  
  int M = absorbed.size();

  for(int j = 0; j < M; j++) {
    
    Rcpp::IntegerVector::iterator i = std::find(id.begin(), id.end(), absorbed[j]);
    
    int index = std::distance(id.begin(), i);
    
    // convert ramet id if it's been absorbed by another ramet
    id[index] = to_sub[j];
    
  } 
  
  return id;
}