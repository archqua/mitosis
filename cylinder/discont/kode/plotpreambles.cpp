#include "plotpreambles.hpp"

const char* leftHalfPreamble
  = "set size 0.5, 0.89\n"
    "set origin 0, 0.11\n"
    ;
const char* contourPreamble
  = "set contour\n"
    "set title 'contour'\n"
    "unset surface\n"
    "set view map\n"
    /* don't quite get all the detail */
    "splot '-' using 1:2:3 with lines\n"
    ;
const char* contourLabelsPreamble
  = "set contour\n"
    "set title 'Contour'\n"
    "unset surface\n"
    "set view map\n"
    /* don't quite get all the detail */
    "splot '-' using 1:2:3 with lines notitle, '-' using 1:2:3 with labels rotate by 45 notitle\n"
    ;

const char* northEastPreamble
  = "set size 0.5, 0.445\n"
    "set origin 0.5, 0.555\n"
    ;
const char* radialSlicePreamble
  = "set title 'radial central slice'\n"
    "set grid\n"
    "plot '-' using 1:2 with lines notitle\n"
    ;

const char* southEastPreamble
  = "set size 0.5, 0.445\n"
    "set origin 0.5, 0.11\n"
    ;
const char* axialSlicePreamble
  = "set title 'axial central slice'\n"
    "set grid\n"
    "plot '-' using 1:2 with lines notitle\n"
    ;

