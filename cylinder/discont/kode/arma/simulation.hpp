#ifndef BSIMULATION_HPP
#define BSIMULATION_HPP
#include "gri2hd.hpp"
#include "parameters.hpp"
#include <armadillo>
#include <string>
//#include <iostream>
#include <iosfwd>

namespace bulk {

extern const arma::mat::fixed<6,6> E_mx;
extern const arma::mat::fixed<6,6> diffusionMask;
extern const arma::mat::fixed<6,6> dbulk_mx;
extern const arma::mat::fixed<6,6> dmt_mx;
extern const arma::mat::fixed<6,6> dav_mx;
extern const arma::mat::fixed<6,6> dzer_mx;

struct GridSpec{
  double ht, hr, hz;
};
struct SizeSpec{
  double T, R, Z;
};
struct PathSpec{
  std::string in, out, log, plots;
};

class Simulation{
  /* parameters 
   * should be explicitly initialized at construction */
  /*       timestep  gridsteps along R and Z */
  const GridSpec grid { .ht = dflt::ht, .hr = dflt::hr, .hz = dflt::hz };
  /*     simulation time    cell radius and length */
  const SizeSpec size { .T = dflt::T, .R = dflt::R, .Z = dflt::Z };
  /* SOR parameter  relaxation stop criterion */
  const double omg = dflt::omg, eps = dflt::eps;
  /* MT radius */
  const double RMT = dflt::RMT;
  const int each;
  const PathSpec path { .in = "", .out = "", .log = "" , .plots = ""};
  /* enums are keept for logging */
  const mode::Jacobian jmode = mode::Jacobian::on;
  const mode::Reactions rmode = mode::Reactions::on;

  /* these are initialized via enums at construction */
  const Jacobian jac;
  const RHS chem;
  /* Jacobian == function<mat<6,6>(const vec<6>&, int i, int j)> 
   * RHS == function<vec<6>(const vec<6>)> */

  /* derived from parameters 
   * mustn't be initialized explicitly */
  const int timesteps = size.T / grid.ht, I = size.R / grid.hr, J = size.Z / grid.hz; //will probably work fine
  const int rmti = RMT / grid.hr;
  //const double invomg = 1./omg;

public:
  Simulation(const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int each
      , const PathSpec& ps
      , mode::Jacobian// = mode::Jacobian::on
      , mode::Reactions// = mode::Reactions::on
      );
  void run();

protected:
  /* full timestep */
  void step();
  /* print current active aurB concentration field to file */
  void stamp(std::ostream& os);
  /* log simulation parameters and max relaxation steps per timestep */
  void log(std::ostream& os);

private:
  template <std::size_t n>
  using vec = arma::vec::fixed<n>;
  template <std::size_t i, std::size_t j>
  using mat = arma::mat::fixed<i, j>;

  /* relaxation iteration
   * return value is stop_relaxations (satisfied) */
  bool relax();
  /* particular node-cell relaxation
   * returns norm of change during update */
  double update(int i, int j, double worstshift);

  /* computation participants */
  Gri2hd<vec<6>> u_ij = Gri2hd<vec<6>>(I+2, J+2); //u_ij[i][j] stands for u_i-0.5,j-0.5
  /* time-independent */
  //Gri2hd<sp_mat> BJL_ij = Gri2hd<sp_mat>(I, J);
protected:
  Gri2hd<mat<6,6>> BJL_ij = Gri2hd<mat<6,6>>(I, J);
private:
  /* time-dependent */
  Gri2hd<mat<6,6>> BJF_ij = Gri2hd<mat<6,6>>(I, J);
  Gri2hd<mat<6,6>> B_ij = Gri2hd<mat<6,6>>(I, J);
  Gri2hd<mat<6,6>>Binv_ij = Gri2hd<mat<6,6>>(I, J);
  /* time-independent */
  std::vector<mat<6,6>> C_i = std::vector<mat<6,6>>(I);
  std::vector<mat<6,6>> A_i = std::vector<mat<6,6>>(I);
  mat<6,6> F_i[3];
  mat<6,6> G_i[3];
  //1D->2D
protected:
  const mat<6,6>& C_ij(int i, int j);
  const mat<6,6>& A_ij(int i, int j);
  const mat<6,6>& F_ij(int i, int j);
  const mat<6,6>& G_ij(int i, int j);
  bool r_in_mt(int i, int j);
  bool r_in_bulk(int i, int j);
  //bool r_in_av(int i, int j);//is incorrect and not called
  bool r_in_vanish(int i, int j);
  const mat<6,6>& diffr_coeff(int i, int j); // diffr_coeff(i,j) stands for d_i,j+0.5
  bool z_in_mt(int i, int j);
  bool z_in_bulk(int i, int j);
  //bool z_in_av(int i, int j);//not called
  bool z_in_vanish(int i, int j);
  const mat<6,6>& diffz_coeff(int i, int j); // diffz_coeff(i,j) stands for d_i+0.5,j
private:
  Gri2hd<vec<6>> dt_phi = Gri2hd<vec<6>>(I, J);
  
  /* logging */
  int maxRelaxationsPerTimestep = 0;

  /* initialization supplements */
  void set_init_cond_from_file(const std::string& path);
  void iniTimeIndep();

  /* automatic plotting */
protected:
  void draw(int t, int n_len);
  /* missing feature:
   * various drawing modes*/
};

void runSimulation(const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int each
      , const PathSpec& ps
      , mode::Jacobian
      , mode::Reactions
      );
void simParametersNotify(std::ostream& os
      , const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int each
      , const PathSpec& ps
      , mode::Jacobian
      , mode::Reactions
      );


}

#endif
