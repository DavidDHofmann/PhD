#include "Rcpp.h"
#include "math.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix interpolatePointsC(
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
