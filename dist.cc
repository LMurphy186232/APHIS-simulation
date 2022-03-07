#include <Rcpp.h>
// This calculates distance matrixes for a previously declared big memory 
// matrix. Given vectors of target coordinates and source coordinates, this
// will populate the matrix with the integer distance between the two. This
// respects max matrix size and will skip distances over the limit. 
//
// This also calculates direction from target to source, and bins it into bins
// of 30 degrees. Azimuth is 0 to 2PI, 0 = north, increasing clockwise. 
// 
// This approach will be the most straightforward method - calculate direction, 
// reduce to bin. Other approaches will be tested for speed.
//
// This is largely following an example in the Rcpp gallery.


// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>
using namespace Rcpp;

//This line relates to the MatrixAccessor term - handles conversion
//from R to C++. You need to have a T for each type of matrix (float, char, etc)
//This one uses integer and char types
template <typename T1, typename T2>
void dodist(XPtr<BigMatrix> pMat,   //Distance matrix 
            XPtr<BigMatrix> bMat,   //Beta index matrix
            XPtr<BigMatrix> ddMat,  //Direction index matrix
            MatrixAccessor<T1> distmat, //Distance matrix 
            MatrixAccessor<T2> betamat, //Beta index matrix
            MatrixAccessor<T2> dirmat,  //Direction index matrix
            NumericVector targetx, NumericVector targety, 
            NumericVector neighx, NumericVector neighy, 
            NumericVector tyear, NumericVector syear, int maxdist) {

    float d, distx, disty, azimuth,
          azichunk = (2*M_PI)/12; //Size of azimuth bins in radians
    int nrow = pMat->nrow(), ncol = pMat->ncol(), nsource = neighx.length(), col = 0;
    try {
        if (targetx.length() > nrow) {
            stop("Too many targets");
        }
        if (neighx.length() != neighy.length()) stop("Neighbor coordinates don't match");
        if (neighx.length() != syear.length()) stop("syear doesn't match");
        if (targetx.length() != targety.length()) stop("Target coordinates don't match");
        if (targetx.length() != tyear.length()) stop("Invalid tyear");
        
        for (size_t i=0; i < nrow; i++) {
            col = 0;
            for (size_t j = 0; j < nsource; j++) {  
                //-----------------------------------------------------------//
                // Calculate distance 
                //-----------------------------------------------------------//
                // Start with this in a float
                d = sqrt(pow((targetx[i] - neighx[j]), 2) + pow((targety[i] - neighy[j]), 2));                
                if (d <= maxdist) {
                    if (col == ncol) {
                        stop("Invalid matrix index");
                    }
                    if (d == 0) {
                        distmat[col][i] = NA_INTEGER;
                        betamat[col][i] = NA_REAL;
                    } else {
                        distmat[col][i] = round(d);
                        
                        //----------------------------------------------------//
                        // Get the beta index based on year of infestation
                        // relative to year of observation
                        //----------------------------------------------------//
                        // This should be 0-indexed
                        betamat[col][i] = tyear[i] - syear[j] - 1;
                        // It is possible for this to go higher than 5 (4 when
                        // zero-indexing) if the year gap is larger than 5, so:
                        if (betamat[col][i] > 4) {
                            betamat[col][i] = 4;
                        }
                    }
                    
                    //--------------------------------------------------------//
                    // Since it was a valid distance, calculate direction
                    //--------------------------------------------------------//
                    //Calculate the azimuth - correct for quadrant of "to" 
                    // relative to "from" (counting clockwise from upper right)
                    distx = neighx[j] - targetx[i];
                    disty = neighy[j] - targety[i];
                    azimuth = 0;
                    if ((disty > 0) && (distx >= 0)) //first quadrant
                        azimuth = atan((distx) / (disty));
                    if ((distx > 0) && (disty <= 0)) //second quadrant
                        azimuth = (M_PI / 2.0) + atan((-1.0 * disty) / (distx));
                    if ((distx <= 0) && (disty < 0)) //third quadrant
                        azimuth = M_PI + atan((-1.0 * distx) / (-1.0 * disty));
                    if ((distx < 0) && (disty >= 0)) //fourth quadrant
                        azimuth = (1.5 * M_PI) + atan((disty) / (-1.0 * distx));
                    
                    dirmat[col][i] = trunc(azimuth / azichunk);
                    
                    // Advance the column counter to the next spot
                    col = col + 1;
                }
            }
           // printf("Row %d of %d\n", i, nrow);
        }    
    } catch(std::exception &ex) {	
	    forward_exception_to_r(ex);
    } catch(...) { 
	    ::Rf_error("c++ exception (unknown reason)"); 
    }
}

   
    
// Dispatch function for distance calculations
//
// [[Rcpp::export]]
void dodist(SEXP pdistMat, SEXP pbetaMat, SEXP pdirMat,
    NumericVector targetx, NumericVector targety, NumericVector neighx, 
    NumericVector neighy, NumericVector tyear, NumericVector syear, int maxdist) {
    // First we have to tell Rcpp what class to use for big.matrix objects.
    // This object stores the attributes of the big.matrix object passed to it
    // by R.
    XPtr<BigMatrix> xpdistMat(pdistMat);
    XPtr<BigMatrix> xpbetaMat(pbetaMat);
    XPtr<BigMatrix> xpdirMat(pdirMat);
    
    if (xpdistMat->matrix_type() != 2) {
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
   // if (xpbetaMat->matrix_type() != 1) {
   //     throw Rcpp::exception("unknown type detected for big.matrix object!");
   // }

    return dodist(xpdistMat, xpbetaMat, xpdirMat, MatrixAccessor<short>(*xpdistMat), 
        MatrixAccessor<char>(*xpbetaMat), MatrixAccessor<char>(*xpdirMat),
        targetx, targety, neighx, neighy, tyear, syear, maxdist);
}