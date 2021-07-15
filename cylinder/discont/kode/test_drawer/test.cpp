//#include "drawer.hpp"
#include <string>
#include <vector>
#include <tuple>
#include <cstdio>

#include "gnuplot-iostream.h"


using namespace std;

string plotname1 = "plot1.png";

Gnuplot& setFileName(Gnuplot& gp){
  gp << "set output '" << plotname1 << "'\n";
  return gp;
}

string preamble =
  "set term png size 1920, 1080\n"
  /* contour takes half of the screen */
  "set size 0.5, 1\n"
  "set contour\n"
  "unset surface\n"
  "set view map\n"
  /* don't quite get all the detail */
  "splot '-' using 1:2:3 with lines\n"
  ;

string toPlot =
  "0 0 0\n"
  "0 1 0\n"
  "0 2 0\n"
  "\n"
  "1 0 2\n"
  "1 1 2\n"
  "1 2 2\n"
  "\n"
  "2 0 1\n"
  "2 1 1\n"
  "2 2 1\n"
  "\n"
  ;

vector<vector<tuple<double, double, double>>>
  formatteData = { {{ 0, 0, 0 }, { 0, 1, 2 }, { 0, 0, 0 }}
                 , {{ 1, 1, 1 }, { 0, 1, 2 }, { 2, 2, 2 }}
                 , {{ 2, 2, 2 }, { 0, 1, 2 }, { 1, 1, 1 }}
                 };
int main() {
  Gnuplot gp;
  setFileName(gp) << preamble;
  //gp << preamble;
  //gp << toPlot;
  gp.send2d(formatteData);

  return 0;
}
