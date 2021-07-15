#include "drawer.hpp"

#include <fstream>
#include <utility>

using namespace std;
namespace bulk {
namespace draw {

SliceContour::SliceContour
  (Gri2hd<double> copy)
  : data(move(copy))
{}


const char* contourHalfTheScreenPreamble
    /* my screen resolution */
  = "set term png size 1920, 1080\n"
    /* contour takes half of the screen */
    "set size 0.5, 1\n"
    "set contour\n"
    "unset surface\n"
    "set view map\n"
    /* don't quite get all the details */
    "splot '-' using 1:2:3 with lines\n"
    ;

Gnuplot& operator<< (Gnuplot& gp, const SliceContour& sc)
  /* output file must be already set at this point */
{
  gp << contourHalfTheScreenPreamble;
  //cerr << contourHalfTheScreenPreamble;
  int I = sc.data.size();
  int J = sc.data[0].size();
  for (int i = 0; i < I; ++i){
    for (int j = 0; j < J; ++j){
      gp << i << " " << j << " " << sc.data[i][j] << "\n";
      //cerr << i << " " << j << " " << data[i][j] << "\n";
    }
    gp << "\n";
    //cerr << "\n";
  }
  //TODO: slices, refinements
  return gp;
}

}
}
