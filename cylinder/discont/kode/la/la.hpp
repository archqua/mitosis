#pragma once
#ifndef _LIN_ALG_HPP
#define _LIN_ALG_HPP
#include <array>
#include <utility>
#include <initializer_list>
#include <armadillo>
#include <iostream>

#define _LIN_ALG_ABS_(num) ((num) > 0 ? (num) : -(num))

namespace la {

template <std::size_t length>
class Vector {
    std::array<double, length> data = {0};
  public:
    double* begin() {
      return data.begin();
    }
    const double* begin() const {
      return data.begin();
    }
    double* end() {
      return data.end();
    }
    const double* end() const {
      return data.end();
    }
    double& operator[] (std::size_t ind) {
      return data[ind];
    }
    const double& operator[] (std::size_t ind) const {
      return data[ind];
    }
    Vector<length>& operator+= (const Vector<length>& other) {
      auto thit = begin();
      for (auto it = other.begin(); it != other.end(); ++it) {
        *thit++ += *it;
      }
      return *this;
    }
    Vector<length> operator+ (const Vector<length>& other) const {
      Vector<length> temp = *this;
      return temp += other;
    }
    Vector<length>& operator*= (double coeff) {
      for (auto it = begin(); it != end(); ++it){
        *it *= coeff;
      }
      return *this;
    }
    Vector<length> operator* (double coeff) const {
      Vector<length> temp = *this;
      return temp *= coeff;
    }
    Vector<length>& operator-= (const Vector<length>& other) {
      auto thit = begin();
      for (auto it = other.begin(); it != other.end(); ++it){
        *thit++ -= *it;
      }
      return *this;
    }
    Vector<length> operator- (const Vector<length>& other) const {
      Vector<length> temp = *this;
      return temp -= other;
    }
    Vector<length> operator- () const {
      Vector<length> temp = *this;
      return temp *= -1;
    }
    double dot(const Vector<length>& other) const {
      double res = 0;
      auto thit = begin();
      for (auto it = other.begin(); it != other.end(); ++it){
        res += (*thit++) * *it;
      }
      return res;
    }
    double norm() const {
      double res = 0;
      for (auto it = begin(); it != end(); ++it) {
        res = res > _LIN_ALG_ABS_(*it) ? res : _LIN_ALG_ABS_(*it);
      }
      return res;
    }
    std::array<double, length>& operator* () {
      return data;
    }
    const std::array<double, length>& operator* () const {
      return data;
    }
};
template <std::size_t length>
Vector<length> operator* (double coeff, Vector<length> vec) {
  return std::move(vec)*coeff;
}
/* avoid writing constructors */
template <std::size_t length>
Vector<length> make_vec(std::initializer_list<double> l){
  Vector<length> res;
  auto thit = res.begin();
  auto it = l.begin();
  while (thit != res.end() && it != l.end()){
    *thit++ = *it++;
  }
  return res;
}
template <std::size_t length, class Range>
Vector<length> make_vec(Range l){
  Vector<length> res;
  auto thit = res.begin();
  auto it = l.begin();
  while (thit != res.end() && it != l.end()){
    *thit++ = *it++;
  }
  return res;
}
template<std::size_t length>
std::ostream& operator<< (std::ostream& os, const Vector<length>& v) {
  bool first = true;
  for (double e : v) {
    if (first) {
      first = false;
    } else {
      os << " ";
    }
    os << e;
  }
  return os;
}

template <std::size_t dim>
class Matrix {
    std::array<Vector<dim>, dim> data;
  public:
    Vector<dim>* begin() {
      return data.begin();
    }
    const Vector<dim>* begin() const {
      return data.begin();
    }
    Vector<dim>* end() {
      return data.end();
    }
    const Vector<dim>* end() const {
      return data.end();
    }
    Vector<dim>& operator[] (std::size_t ind) {
      return data[ind];
    }
    const Vector<dim>& operator[] (std::size_t ind) const {
      return data[ind];
    }
    Vector<dim> operator* (const Vector<dim>& vec) const {
      Vector<dim> res;
      for (std::size_t i = 0; i != dim; ++i) {
        res[i] = vec.dot(data[i]);
      }
      return res;
    }
    Matrix<dim>& operator*= (double coeff) {
      for (auto it = begin(); it != end(); ++it) {
        *it *= coeff;
      }
      return *this;
    }
    Matrix<dim> operator* (double coeff) const {
      Matrix<dim> temp = *this;
      return temp *= coeff;
    }
    Matrix<dim>& operator+= (const Matrix<dim>& other) {
      auto thit = begin();
      for (auto it = other.begin(); it != other.end(); ++it) {
        *thit++ += *it;
      }
      return *this;
    }
    Matrix<dim> operator+ (const Matrix<dim>& other) const {
      Matrix<dim> temp = *this;
      return temp += other;
    }
    Matrix<dim>& operator-= (const Matrix<dim>& other) {
      auto thit = begin();
      for (auto it = other.begin(); it != other.end(); ++it){
        *thit++ -= *it;
      }
      return *this;
    }
    Matrix<dim> operator- (const Matrix<dim>& other) const {
      Matrix<dim> temp = *this;
      return temp -= other;
    }
    Matrix<dim> operator- () const {
      Matrix<dim> temp = *this;
      return temp *= -1;
    }
    std::array<Vector<dim>, dim>& operator* () {
      return data;
    }
    const std::array<Vector<dim>, dim>& operator* () const {
      return data;
    }

    Matrix<dim> i() const {
      arma::mat::fixed<dim, dim> inter; inter.fill(0.);
      for (std::size_t j = 0; j < dim; ++j) {
        for (std::size_t k = 0; k < dim; ++k){
          inter(j,k) = (*this)[j][k];
        }
      }
      inter = inter.i();
      Matrix<dim> res;
      for (std::size_t j = 0; j < dim; ++j) {
        for (std::size_t k = 0; k < dim; ++k){
          res[j][k] = inter(j,k);
        }
      }
      return res;
    }
};
template <std::size_t dim>
Matrix<dim> operator* (double coeff, Matrix<dim> mat) {
  return std::move(mat) *= coeff;
}

template <std::size_t dim>
Matrix<dim> eye() {
  Matrix<dim> res;
  for (std::size_t i = 0; i < dim; ++i) {
    res[i][i] = 1;
  }
  return res;
}
template <std::size_t dim>
Matrix<dim> make_mat(std::initializer_list<Vector<dim>> l) {
  Matrix<dim> res;
  auto thit = res.begin();
  auto it = l.begin();
  while (thit != res.end() && it != l.end()){
    *thit++ = *it++;
  }
  return res;
}
template <std::size_t dim, class Range>
Matrix<dim> make_mat(Range l) {
  Matrix<dim> res;
  auto thit = res.begin();
  auto it = l.begin();
  while (thit != res.end() && it != l.end()){
    *thit++ = *it++;
  }
  return res;
}
template <std::size_t dim>
Matrix<dim> inv(const Matrix<dim>& mat) {
  return mat.i();
}

template<std::size_t dim>
std::ostream& operator<< (std::ostream& os, const Matrix<dim>& m) {
  bool first = true;
  for (const Vector<dim>& v : m) {
    if (first) {
      first = false;
    } else {
      os << "\n";
    }
    os << v;
  }
  return os;
}

}//namespace la


#endif
