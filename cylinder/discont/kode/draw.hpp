#pragma once
#ifndef BDRAWER_HPP
#define BDRAWER_HPP

#include "gri2hd.hpp"

#include <string>
#include <memory>
#include <utility>
#include <iostream>

#include "gnuplot-iostream.h"//-I/home/tor/gitlibs/gnuplot-iostream

template <class Num>
//Gnuplot& operator<< (Gnuplot& gp, const Gri2hd<Num>& grid) {
std::ostream& operator<< (std::ostream& gp, const Gri2hd<Num>& grid) {
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
template <class Num>
//gnuplotio::Gnuplot& operator<< (gnuplotio::Gnuplot& gp, const typename Gri2hd<Num>::Col& col) {
std::ostream& operator<< (std::ostream& gp, const typename Gri2hd<Num>::Col& col) {
  const auto& cpar = col.params();
  double yinit = cpar.yinit;
  double hy = cpar.hy;
  for (size_t j = 0; j < col.size(); ++j) {
    gp << (yinit + j * hy) << " " << col[j] << "\n";
  }
  gp << "\n";
  return gp;
}
template <class Num>
//gnuplotio::Gnuplot& operator<< (gnuplotio::Gnuplot& gp, const typename Gri2hd<Num>::Row& row) {
std::ostream& operator<< (std::ostream& gp, const typename Gri2hd<Num>::Row& row) {
  const auto& rpar = row.params();
  double xinit = rpar.xinit;
  double hx = rpar.hx;
  for (size_t i = 0; i < row.size(); ++i) {
    gp << (xinit + i * hx) << " " << row[i] << "\n";
  }
  gp << "\n";
  return gp;
}


#endif
