////////////////////////////////////////////////////////////////////////////////
//// Helper Functions to Simulate Dispersal
////////////////////////////////////////////////////////////////////////////////
#include "Rcpp.h"
#include "math.h"
using namespace Rcpp;
using namespace std;
//' Calculate New Absolute Turning Angle
//'
//' Function to calculate new absolute turning angle
//' @export
//' @param absta Absolute turning angle
//' @param ta Relative turning angle
//' @return numeric value
// [[Rcpp::export]]
NumericVector getAbsNew(double absta, NumericVector ta) {

  // Prepare loop counter based on vector length
  int n = ta.length();

  // Calculate new absolute turning angle
  NumericVector absta_new(n);
  for (int i = 0; i < n; i++){
    absta_new[i] = absta + ta[i];

    // We need to make sure the turning angle is between 0 and 2 Pi
    if (absta_new[i] > 2 * M_PI){
      absta_new[i] = absta_new[i] - 2 * M_PI;
    } else if (absta_new[i] < 0){
      absta_new[i] = absta_new[i] + 2 * M_PI;
    }
  }

  // Return it
  return absta_new;
}

//' Calculate New Endpoints
//'
//' Function to calculate new endpoints
//' @export
//' @param xy \code{matrix} containing xy coordinates
//' @param absta absolute turning angle
//' @param sl step lengths
//' @return numeric value
// [[Rcpp::export]]
NumericMatrix calcEndpoints(
    NumericMatrix xy    // Coordinates of start point
  , NumericVector absta // Absolute turning angle
  , NumericVector sl    // Step length
){

  // Prepare loop counter based on vector length
  int n = sl.length();

  // Make sure the step length covers at least a meter and convert to degrees
  NumericVector sl_new(n);
  for (int i = 0; i < n; i++){
    if (sl[i] < 1) {
      sl_new[i] = 1 / 111000;
    } else {
      sl_new[i] = sl[i] / 111000;
    }
  }

  // Prepare matrix into which the new coordinates will go
  NumericMatrix xy_new(n, 2);
  for (int i = 0; i < n; i++) {
    xy_new(i, 0) = xy(0, 0) + sin(absta[i]) * sl_new[i];
    xy_new(i, 1) = xy(0, 1) + cos(absta[i]) * sl_new[i];
  }

  // Return the results
  return xy_new;
}

//' Interpolate Between two Points
//'
//' Function to interpolate coordinates between two points
//' @export
//' @param x1 numeric containing x coordiante of first point
//' @param x2 numeric containing x coordiante of second point
//' @param y1 numeric containing y coordiante of first point
//' @param y2 numeric containing y coordiante of second point
//' @param by numeric Approximate distance at which points should be interpolated
//' @return matrix
// [[Rcpp::export]]
NumericMatrix interpolatePoints(
      const double x1
    , const double x2
    , const double y1
    , const double y2
    , const double by
  ){

  // Calculate length of line
  const double length = sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1));

  // Calculate number segments required for "by"
  int nsegs = ceil(length / by);

  // Need at least 1 segment for a line
  nsegs = max(nsegs, 1);

  // Calculate delta x and delta y for a single segment
  double del_x = (x2 - x1) / nsegs;
  double del_y = (y2 - y1) / nsegs;

  // Calculate coordinates of interpolated points
  NumericMatrix xy_new(nsegs + 1, 2);
  for (int i = 0; i <= nsegs; i++){
    xy_new(i, 0) = x1 + i * del_x;
    xy_new(i, 1) = y1 + i * del_y;
  }

  // return the new coordinates
  return xy_new;
}
