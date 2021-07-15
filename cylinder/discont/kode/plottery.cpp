#include "plottery.hpp"

using namespace std;

namespace bulk {
namespace plottery {
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

//const char* annotativeLabelNames[]
//  = { "time"
//    , "ht"
//    , "hr"
//    , "hz"
//    , "omega"
//    , "epsilon"
//    };
#define LABEL_FRACTION_OF_SCREEN 0.11
void setAnnotativeLabels(stringstream& gp, const AnnotationInfo& ainfo) {
  const vector<string>& lNames = ainfo.annotativeLabelNames;
  const vector<double>& lValues = ainfo.annotativeLabelValues;
  for (size_t i = 0; i < ainfo.size; ++i) {
    gp << "set label " << (i+1) << " " <<
          "'" << lNames[i] << ": " << lValues[i] << "' " <<
          "at screen " << (i+1)*1. / (ainfo.size+1) << ", " <<
             "screen " << 0.5 * LABEL_FRACTION_OF_SCREEN << "\n";
  }
}

void unsetAnnotativeLabels(stringstream& gp, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    gp << "unset label " << (i+1) << "\n";
  }
}

stringstream plotUsual(Gri2hd<double> data, AnnotationInfo ainfo){
  /* my screen resolution */
  stringstream gp;
  gp << "set term png size 1920, 1080 enhanced\n"
        "set encoding utf8\n"
        "set multiplot\n";
  {//main contour at the left
    gp << "set xlabel 'R'\n" << "set ylabel 'Z'\n";
    gp << leftHalfPreamble << contourLabelsPreamble;
    gp << data << "e" << "\n";
    /* labels require repeating */
    gp << data << "e" << "\n";
  }
  {//central axial slice at the top right
    gp << "set xlabel 'Z'\n" << "set ylabel 'active AurB'\n";
    const auto& row = data.row(0);
    gp << northEastPreamble << axialSlicePreamble;
    row.send(gp) << "e" << "\n";
  }
  setAnnotativeLabels(gp, ainfo);
  {//central radial slice at the bottom right
    gp << "set xlabel 'R'\n" << "set ylabel 'active AurB'\n";
    const auto& col = data.col(data[0].size() / 2);
    gp << southEastPreamble << radialSlicePreamble;
    col.send(gp) << "e" << "\n";
  }
  unsetAnnotativeLabels(gp, ainfo.size);
  //TODO
  return gp;
}

}//namespace plottery
}//namespace bulk
