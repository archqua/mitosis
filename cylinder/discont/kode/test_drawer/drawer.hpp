#pragma once
#ifndef BDRAWER_HPP
#define BDRAWER_HPP

#include "gri2hd.hpp"

#include "gnuplot-iostream.h"//-I/home/tor/gitlibs/gnuplot-iostream

#include <string>
#include <memory>
#include <utility>

namespace bulk {
namespace draw {

class SliceContour {
  friend Gnuplot& operator<< (Gnuplot&, const SliceContour&);
    Gri2hd<double> data;
  public:
    SliceContour(Gri2hd<double> copy);
};
Gnuplot& operator<< (Gnuplot& gp, const SliceContour& sc);

}
}

#endif
