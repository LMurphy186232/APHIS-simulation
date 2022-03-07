#include <Rcpp.h>

#include <numeric>
using namespace Rcpp;

// This counts neighbor trees within a certain distance.
Rcpp::NumericVector countacers(NumericVector targetx, NumericVector targety, double maxdist) {
  int ntrees = targetx.length();
  Rcpp::NumericVector count(ntrees);
  float xdiff, ydiff, d;
    
    try {
        
        for (size_t i=0; i < ntrees; i++) {
          
            count[i] = 0;
            
            for (size_t j = 0; j < ntrees; j++) {
                if (j != i) {
                    
                    // Quick test to see if it's worth calculating distance
                    xdiff = abs(targetx[i] - targetx[j]);
                    ydiff = abs(targety[i] - targety[j]);
                    
                    if (xdiff <= maxdist && ydiff <= maxdist) {
                        //---------------------------------------------------//
                        // Calculate distance 
                        //---------------------------------------------------//
                        d = sqrt(xdiff*xdiff + ydiff*ydiff);                
                        if (d <= maxdist) {
                            count[i] = count[i] + 1;                                            
                        }
                    }
                }    
            }
           // printf("Row %d of %d\n", i, nrow);
        }    
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
    return count;
}


// This updates neighbor counts. Given a set of trees that died and coordinates,
// this will tell you which trees to subtract from
Rcpp::NumericVector updateacers(NumericVector targetx, NumericVector targety,
        NumericVector deadx, NumericVector deady, double maxdist) {
  
    int ntrees = targetx.length(), ndead = deadx.length();
    Rcpp::NumericVector count(ntrees);
    float xdiff, ydiff, d;
  
    // Initialize the return with 0s
    for (size_t i=0; i < ntrees; i++) count[i] = 0;
    
    try {
        
        for (size_t i=0; i < ndead; i++) {
          
            for (size_t j = 0; j < ntrees; j++) {
                // Quick test to see if it's worth calculating distance
                xdiff = abs(deadx[i] - targetx[j]);
                ydiff = abs(deady[i] - targety[j]);
                if (xdiff <= maxdist && ydiff <= maxdist) {
                
                    //-------------------------------------------------------//
                    // Calculate distance 
                    //-------------------------------------------------------//
                    d = sqrt(xdiff*xdiff + ydiff*ydiff);               
                    if (d <= maxdist) {
                        count[j] = count[j] + 1;                                            
                    }
                }    
            }
        }    
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
    return count;
}
   
    
// Dispatch function for distance calculations
//
// [[Rcpp::export]]
Rcpp::NumericVector count_acer_neighbors(NumericVector targetx, NumericVector targety, double maxdist) {
    
    return countacers(targetx, targety, maxdist);
}
// [[Rcpp::export]]
Rcpp::NumericVector update_acer_neighbors(NumericVector targetx, NumericVector targety, NumericVector deadx, NumericVector deady, double maxdist) {
    
    return updateacers(targetx, targety, deadx, deady, maxdist);
}
