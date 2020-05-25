#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector negative_subset(NumericVector vec, NumericVector filter) {

  LogicalVector subsets (vec.length(), 0);

  for(int i = 0; i < filter.length(); i++){
    for(int j = 0; j< filter.length(); j++){
      if( filter[j] == i){
        subsets[ filter[i] ] = 1;
      }
    }
  }
 return vec[ !subsets ];
}
