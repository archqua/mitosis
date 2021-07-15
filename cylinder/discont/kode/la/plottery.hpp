#pragma once
#ifndef BPLOTTERY_HPP
#define BPLOTTERY_HPP
#include "gri2hd.hpp"
#include <sstream>
#include <vector>

namespace bulk {
namespace plottery {

struct AnnotationInfo {
  std::size_t size;
  const std::vector<std::string> annotativeLabelNames;
  const std::vector<double> annotativeLabelValues;
};

std::stringstream plotUsual(Gri2hd<double> data, AnnotationInfo ainfo);

}//namespace plottery
}//namespace bulk

#endif//BPLOTTERY_HPP guard
