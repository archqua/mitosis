#pragma once
#ifndef BDRAWER_HPP
#define BDRAWER_HPP

#include "gri2hd.hpp"

#include "gnuplot-iostream.h"//-I/home/tor/gitlibs/gnuplot-iostream

#include <string>
#include <memory>
#include <utility>

template <class Num>
Gnuplot& operator<< (Gnuplot& gp, const Gri2hd<Num>& grid) {
  const auto& grpar = grid.params();
  double xinit = grpar.xinit;
  double yinit = grpar.yinit;
  double hx = grpar.hx;
  double hy = grpar.hy;
  for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j = 0; j < grid.size(); ++j) {
      gp << (xinit + i * hx) << " " << (yinit + j * hy) << " " << grid[i][j] << "\n";
    }
    gp << "\n";
  }
  return gp;
}

template <class Data>
class Draw {
  protected:
    std::vector<const char*> preambles;
    std::vector<Data> datas;
  public:
    virtual ~Draw() = default;
    Draw(std::vector<const char*> pcopies, std::vector<Data> dcopies)
      : preambles(move(pcopies))
      , datas(move(dcopies))
    {}
    /* output file must be already set at this point */
    virtual Gnuplot& gprint(Gnuplot& gp) {
      for (std::size_t i = 0; i < preambles.size() && i < datas.size(); ++i) {
        gp << preambles[i];
        gp << datas[i];
      }
      return gp;
    }
};
/* no need for friend declaration??? */
template <class Data>
Gnuplot& operator<< (Gnuplot& gp, const Draw<Data>& drawer) {
  return drawer.gprint(gp);
}

namespace draw {

template <class Num>
Gnuplot& operator<< (Gnuplot& gp, const Gri2hd<Num>& grid) {
  int I = grid.size();
  int J = grid[0].size();
  for (int i = 0; i < I; ++i){
    for (int j = 0; j < J; ++j){
      gp << i << " " << j << " " << grid[i][j] << "\n";
      //cerr << i << " " << j << " " << data[i][j] << "\n";
    }
    gp << "\n";
    //cerr << "\n";
  }
  return gp;
}

class CylContourWithSlices : public Draw<Gri2hd<double>> {
  public:
    CylContourWithSlices(Gri2hd<double> copy)
      : Draw<Gri2hd<double>>( { contourHalfTheScreenPreamble }
                            , { std::move(copy) }) 
      {}
    
};
//Gnuplot& operator<< (Gnuplot& gp, const CylContourWithSlices&);

}

#endif
