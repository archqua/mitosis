#include "simulation.hpp"
#include <cassert>
#include <iostream>

#define DIFFRC(i, j) ((arma::mat::fixed<6,6>)diffr_coeff(i, j))(0, 0)
#define DIFFZC(i, j) ((arma::mat::fixed<6,6>)diffz_coeff(i, j))(0, 0)
#define PURE_DIFFUSION_TEST \
{\
    assert(DIFFRC(0, 0) == 0);\
    assert(DIFFRC(0, 4) == 0);\
    assert(DIFFRC(0, 9) == 0);\
    assert(DIFFRC(10, 0) == 0);\
    assert(DIFFRC(10, 4) == 0);\
    assert(DIFFRC(10, 9) == 0);\
    assert(DIFFZC(0, 0) == 0);\
    assert(DIFFZC(4, 0) == 0);\
    assert(DIFFZC(9, 0) == 0);\
    assert(DIFFZC(0, 10) == 0);\
    assert(DIFFZC(4, 10) == 0);\
    assert(DIFFZC(9, 10) == 0);\
}

#define AA(i,j) A_ij(i, j)(0, 0)
#define BB(i,j) B_ij(i, j)(0, 0)
#define BJL(i,j) BJL_ij[i][j](0, 0)
#define CC(i,j) C_ij(i, j)(0, 0)
#define GG(i,j) G_ij(i, j)(0, 0)
#define FF(i,j) F_ij(i, j)(0, 0)
#define ABCD_TEST \
{\
    assert(AA(0, 0) == 0);\
    assert(AA(0, 4) == 0);\
    assert(AA(0, 9) == 0);\
    assert(CC(9, 0) == 0);\
    assert(CC(9, 4) == 0);\
    assert(CC(9, 9) == 0);\
    assert(GG(0, 0) == 0);\
    assert(GG(4, 0) == 0);\
    assert(GG(9, 0) == 0);\
    assert(FF(0, 9) == 0);\
    assert(FF(4, 9) == 0);\
    assert(FF(9, 9) == 0);\
}

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define EPS 1e-05
#define DISTORSION(i, j) (BJL_ij[i][j](0,0) - (i == j) - AA(i,j) - CC(i,j) - FF(i,j) - GG(i,j))
#define UNIFORMITY_TEST \
{\
  assert(ABS(DISTORSION(4,4)) == 0);\
  assert(ABS(DISTORSION(9,4)) == 0);\
  assert(ABS(DISTORSION(9,9)) == 0);\
  assert(ABS(DISTORSION(4,9)) == 0);\
  assert(ABS(DISTORSION(4,0)) == 0);\
  assert(ABS(DISTORSION(9,0)) == 0);\
  assert(ABS(DISTORSION(0,0)) == 0);\
  assert(ABS(DISTORSION(0,4)) == 0);\
  assert(ABS(DISTORSION(0,9)) == 0);\
}

namespace bulk {

void printCoeffs(std::ostream& os) {
  /* flow along r */
}

class TestSimulation : public Simulation {
  public:
    TestSimulation(): Simulation( GridSpec{.ht = 0.1, .hr = 1, .hz = 1}
                                , SizeSpec{.T = 1, .R = 4, .Z = 4}
                                , 1, 1e-03
                                , 2
                                , 1
                                , PathSpec{ .in = "stupid4x4ic.scv", .out = "stupid4x4_out.csv"
                                          , .log = "", .plots = "./stupid4x4plots/"}
                                , mode::Jacobian::off
                                , mode::Reactions::off
                                )
    {
//      PURE_DIFFUSION_TEST
//      std::cerr << "PureDiffusionTest OK" << std::endl;
//      ABCD_TEST
//      std::cerr << "ABCDTest OK" << std::endl;
//      UNIFORMITY_TEST
//      std::cerr << "UniformityTest OK" << std::endl;
      std::cerr << "I: " << I << ", J: " << J << "\n";
      for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 4; ++i) {
          std::cerr << DIFFZC(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "-----------\n";
      for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 5; ++i) {
          std::cerr << DIFFRC(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "A-----------\n";
      for (int j = 0; j < 4; ++j){
        for (int i = 0; i < 4; ++i){
          std::cerr << AA(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "C-----------\n";
      for (int j = 0; j < 4; ++j){
        for (int i = 0; i < 4; ++i){
          std::cerr << CC(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "F-----------\n";
      for (int j = 0; j < 4; ++j){
        for (int i = 0; i < 4; ++i){
          std::cerr << FF(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "G-----------\n";
      for (int j = 0; j < 4; ++j){
        for (int i = 0; i < 4; ++i){
          std::cerr << GG(i,j) << " ";
        }
        std::cerr << "\n";
      }
      std::cerr << "B-----------\n";
      for (int j = 0; j < 4; ++j){
        for (int i = 0; i < 4; ++i){
          std::cerr << BJL(i,j) << " ";
        }
        std::cerr << "\n";
      }
      //run();
    }
};

}

int main() {
  bulk::TestSimulation();
  return 0;
}
