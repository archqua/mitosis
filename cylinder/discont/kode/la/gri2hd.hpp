#ifndef GRI2HD_HPP
#define GRI2HD_HPP
#include <vector>
#include <iostream>

/* 2.5 dimension grid
 * that is actually just 2D grid */
struct GridParams {
  double xinit;
  double yinit;
  double hx;
  double hy;
};
template<class Value>
class Gri2hd{
  using size_t = std::size_t;
  using Line = std::vector<Value>;
  std::vector<Line> grid;
  /* i is x index, j is y index */
  struct Params {
    double xinit = 0;
    double yinit = 0;
    double hx = 1;
    double hy = 1;
  } _params;

public:
  Gri2hd() = default;
  Gri2hd(const Gri2hd& other) = default;
  Gri2hd(Gri2hd&& other) = default;
  Gri2hd& operator= (const Gri2hd& other) = default;
  Gri2hd& operator= (Gri2hd&& other) = default;
  Gri2hd(size_t I, size_t J){
    grid.reserve(I);
    for (size_t i = 0; i < I; ++i){
      grid.emplace_back(J);
    }
  }
  Gri2hd( size_t I, size_t J
        , double xinit, double yinit
        , double hx, double hy
        )
    : _params({xinit, yinit, hx, hy})
  {
    grid.reserve(I);
    for (size_t i = 0; i < I; ++i){
      grid.emplace_back(J);
    }
  }
  Line& operator[] (size_t i){
    return grid[i];
  }
  const Line& operator[] (size_t i) const {
    return grid[i];
  }
  /* rows and cols */
  class Col {
      Gri2hd<Value>& _grid;
      size_t j;
      struct Params {
        double xinit;
        double hx;
      };
    public:
      Col(Gri2hd<Value>& grid_, size_t j_)
        : _grid(grid_)
        , j(j_)
      {}
      Value& operator[] (size_t i) {
        return _grid[i][j];
      }
      const Value& operator[] (size_t i) const {
        return _grid[i][j];
      }
      Params params() const {
        return {_grid.params().xinit, _grid.params().hx};
      }
      size_t size() const {
        return _grid[0].size();
      }
      std::ostream& send(std::ostream& gp) const {
        const Params cpar = params();
        double xinit = cpar.xinit;
        double hx = cpar.hx;
        for (size_t i = 0; i < size(); ++i) {
          gp << (xinit + i * hx) << " " << (*this)[i] << "\n";
        }
        //gp << "\n";
        return gp;
      }
  };
//  class constCol {
//      const Gri2hd<Value>& _grid;
//      size_t j;
//    public:
//      const Value& operator[] (size_t i) const {
//        return _grid[i][j];
//      }
//  };
  class Row {
      Gri2hd<Value>& _grid;
      size_t i;
      struct Params {
        double yinit;
        double hy;
      };
    public:
      Row(Gri2hd<Value>& grid_, size_t i_)
        : _grid(grid_)
        , i(i_)
      {}
      Value& operator[] (size_t j) {
        return _grid[i][j];
      }
      const Value& operator[] (size_t j) const {
        return _grid[i][j];
      }
      Params params() const {
        return {_grid.params().xinit, _grid.params().hx};
      }
      size_t size() const {
        return _grid.size();
      }
      std::ostream& send(std::ostream& gp) const {
        const Params rpar = params();
        double yinit = rpar.yinit;
        double hy = rpar.hy;
        for (size_t j = 0; j < size(); ++j) {
          gp << (yinit + j * hy) << " " << (*this)[j] << "\n";
        }
        //gp << "\n";
        return gp;
      }
  };
//  class constRow {
//      const Gri2hd<Value>& _grid;
//      size_t i;
//    public:
//      const Value& operator[] (size_t j) const {
//        return _grid[i][j];
//      }
//  };
  Col col(size_t j) {
    return {*this, j};
  }
//  constCol col(size_t j) const {
//    return {*this, j};
//  }
  Row row(size_t i) {
    return {*this, i};
  }
//  constRow row(size_t i) const {
//    return {*this, i};
//  }

  size_t size() const {
    return grid.size();
  }

  const std::vector<Line>& operator* () const {
    return grid;
  }
  const Params& params() const {
    return _params;
  }

}; //class Gri2hd<Value>

template <class Num>
std::ostream& operator<< (std::ostream& gp, const Gri2hd<Num>& grid) {
//std::ostream& operator<< (std::ostream& gp, const Gri2hd<Num>& grid) {
  const auto& grpar = grid.params();
  double xinit = grpar.xinit;
  double yinit = grpar.yinit;
  double hx = grpar.hx;
  double hy = grpar.hy;
  for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j = 0; j < grid[0].size(); ++j) {
      gp << (xinit + i * hx) << " " << (yinit + j * hy) << " " << grid[i][j] << "\n";
    }
    gp << "\n";
  }
  return gp;
}
//template <class Num>
//gnuplotio::Gnuplot& operator<< (gnuplotio::Gnuplot& gp, const typename Gri2hd<Num>::Col& col) {
////std::ostream& operator<< (std::ostream& gp, const typename Gri2hd<Num>::Col& col) {
//  const auto& cpar = col.params();
//  double yinit = cpar.yinit;
//  double hy = cpar.hy;
//  for (size_t j = 0; j < col.size(); ++j) {
//    gp << (yinit + j * hy) << " " << col[j] << "\n";
//  }
//  gp << "\n";
//  return gp;
//}
//template <class Num>
//gnuplotio::Gnuplot& operator<< (gnuplotio::Gnuplot& gp, const typename Gri2hd<Num>::Row& row) {
////std::ostream& operator<< (std::ostream& gp, const typename Gri2hd<Num>::Row& row) {
//  const auto& rpar = row.params();
//  double xinit = rpar.xinit;
//  double hx = rpar.hx;
//  for (size_t i = 0; i < row.size(); ++i) {
//    gp << (xinit + i * hx) << " " << row[i] << "\n";
//  }
//  gp << "\n";
//  return gp;
//}

#endif
