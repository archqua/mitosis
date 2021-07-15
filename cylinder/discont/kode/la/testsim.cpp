#include "simulation.hpp"
#include <cassert>
#include <iostream>

#define DIFFRC(i, j) diffr_coeff(i, j)[0][0]
#define DIFFZC(i, j) diffz_coeff(i, j)[0][0]
#define PURE_DIFFUSION_TEST \
{\
    assert(DIFFRC(0, 0) == 0);\
    assert(DIFFRC(0, 2) == 0);\
    assert(DIFFRC(0, 3) == 0);\
    assert(DIFFRC(4, 0) == 0);\
    assert(DIFFRC(4, 2) == 0);\
    assert(DIFFRC(4, 3) == 0);\
    assert(DIFFZC(0, 0) == 0);\
    assert(DIFFZC(2, 0) == 0);\
    assert(DIFFZC(3, 0) == 0);\
    assert(DIFFZC(0, 4) == 0);\
    assert(DIFFZC(2, 4) == 0);\
    assert(DIFFZC(3, 4) == 0);\
}

#define AA(i,j) A_ij(i, j)[0][0]
#define BB(i,j) B_ij(i, j)[0][0]
#define CC(i,j) C_ij(i, j)[0][0]
#define GG(i,j) G_ij(i, j)[0][0]
#define FF(i,j) F_ij(i, j)[0][0]
#define ABCD_TEST \
{\
    assert(AA(0, 0) == 0);\
    assert(AA(0, 2) == 0);\
    assert(AA(0, 3) == 0);\
    assert(CC(3, 0) == 0);\
    assert(CC(3, 2) == 0);\
    assert(CC(3, 3) == 0);\
    assert(GG(0, 0) == 0);\
    assert(GG(2, 0) == 0);\
    assert(GG(3, 0) == 0);\
    assert(FF(0, 3) == 0);\
    assert(FF(2, 3) == 0);\
    assert(FF(3, 3) == 0);\
}

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define EPS 1e-05
#define DISTORSION(i, j) (BJL_ij[i][j][0][0] - 1 - AA(i,j) - CC(i,j) - FF(i,j) - GG(i,j))
#define UNIFORMITY_TEST \
{\
  std::cerr << DISTORSION(2,2) << std::endl;\
  std::cerr << DISTORSION(3,2) << std::endl;\
  std::cerr << DISTORSION(3,3) << std::endl;\
  std::cerr << DISTORSION(2,3) << std::endl;\
  std::cerr << DISTORSION(2,0) << std::endl;\
  std::cerr << DISTORSION(2,0) << std::endl;\
  std::cerr << DISTORSION(0,0) << std::endl;\
  std::cerr << DISTORSION(0,2) << std::endl;\
  std::cerr << DISTORSION(0,3) << std::endl;\
  assert(ABS(DISTORSION(2,2)) < 1e-05);\
  assert(ABS(DISTORSION(3,2)) < 1e-04);\
  assert(ABS(DISTORSION(3,3)) < 1e-04);\
  assert(ABS(DISTORSION(2,3)) < 1e-04);\
  assert(ABS(DISTORSION(2,0)) < 1e-04);\
  assert(ABS(DISTORSION(3,0)) < 1e-04);\
  assert(ABS(DISTORSION(0,0)) < 1e-04);\
  assert(ABS(DISTORSION(0,2)) < 1e-04);\
  assert(ABS(DISTORSION(0,3)) < 1e-04);\
}

namespace bulk {

class TestSimulation : public Simulation {
  public:
    TestSimulation(): Simulation( GridSpec{.ht = 0.1, .hr = 1, .hz = 1}
                                , SizeSpec{.T = 1, .R = 4, .Z = 4}
                                , 1, 1e-03
                                , 0
                                , 1
                                , PathSpec{ .in = "stupid4x4ic.scv", .out = "stupid4x4_out.csv"
                                          , .log = "", .plots = "./stupid4x4plots/"}
                                , mode::Jacobian::off
                                , mode::Reactions::off
                                )
    {
      PURE_DIFFUSION_TEST
      std::cerr << "PureDiffusionTest OK" << std::endl;
      ABCD_TEST
      std::cerr << "ABCDTest OK" << std::endl;
      UNIFORMITY_TEST
      std::cerr << "UniformityTest OK" << std::endl;
      run();
    }
};

}

int main() {
  bulk::TestSimulation();
  return 0;
}
