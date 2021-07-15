#include "simulation.hpp"
#include "plottery.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "gnuplot-iostream.h"//-I/home/tor/gitlibs/gnuplot-iostream


using namespace std;
namespace bulk {

Simulation::Simulation(const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int _each
      , const PathSpec& ps
      , mode::Jacobian _jmode = mode::Jacobian::on
      , mode::Reactions _rmode = mode::Reactions::on
      )
  : grid(gs)
    , size(ss)
    , omg(_omg), eps(_eps)
    , RMT(_RMT)
    , each(_each)
    , path(ps) 
    , jmode(_jmode)
    , rmode(_rmode)
    , jac(jacs[static_cast<int>(_jmode)])
    , chem(rhss[static_cast<int>(_rmode)])
{
  if (timesteps % each)
    cerr << "each " << each << "doesn't divide timesteps " << timesteps << endl;
  try {
    set_init_cond_from_file(path.in);
  } catch (const char* msg) {
    cerr << "couldn't set initial conditions" << endl;
    cerr << msg << endl;
    exit(1);
  }
  iniTimeIndep();
}

void Simulation::run(){
  ofstream out(path.out);
  if (!out.is_open()){
    cerr << "unable to save results to " << path.out << endl;
    exit(2);
  }
  //stamp(out);
  for (int t = 0; t < timesteps + 1; ++t){
    if (t % each == 0){
      stamp(out);
      if (path.plots != ""){
        try {
          draw(t, (int) to_string(timesteps).size());
        } catch(...) {
          cerr << "unable to draw plot to " << path.plots << " directory" << endl;
        }
      }
    }
    /* one extra step not recorded */
    step();
  }
  out.close();
  if (path.log != ""){
    ofstream logStream(path.log);
    if (logStream)
      log(logStream);
    else
      cerr << "unable to log to " << path.log << endl;
      /* too late to throw or exit) */
    logStream.close();
  }
}

/* full timestep */
void Simulation::step(){
  const double ht = grid.ht;//, hr = grid.hr, hz = grid.hz;
  for (int i=0; i<I; ++i){
    for (int j=0; j<J; ++j){
      const arma::vec::fixed<6> uij = u_ij[i+1][j+1];
      BJF_ij[i][j] = - ht*jac(uij,i,j); //mind the sign
      B_ij[i][j] = BJL_ij[i][j] + BJF_ij[i][j];
      Binv_ij[i][j] = B_ij[i][j].i();
      /* dangerous place */
      dt_phi[i][j] = ht*chem(uij,i,j) + BJF_ij[i][j]*uij + uij; //SIGN!!!
    }
  }
  int relaxations = 0;
  do{
    ++relaxations;
  } while (!relax());
  maxRelaxationsPerTimestep = std::max(relaxations, maxRelaxationsPerTimestep);
}

/* relaxation iteration */
bool Simulation::relax(){
  int jinit=0;
  double worstshift=0;
  for (int i=0; i<I; ++i){
    for (int j=jinit; j<J; j+=2){
      worstshift = update(i, j, worstshift);
    }
    jinit = (jinit + 1) % 2;
  }
  jinit=1;
  for (int i=0; i<I; ++i){
    for (int j=jinit; j<J; j+=2){
      worstshift = update(i, j, worstshift);
    }
    jinit = (jinit + 1) % 2;
  }
  return worstshift < eps; //relaxations stop at !relax()
}

/* particular node-cell relaxation
 * returns norm of change during update */
double Simulation::update(int i, int j, double worstshift){
  vec<6> uij_prev = u_ij[i+1][j+1];
  u_ij[i+1][j+1] = (1 - omg)*u_ij[i+1][j+1] 
    + omg*Binv_ij[i][j]
       * (dt_phi[i][j]
          + C_ij(i,j)*u_ij[i+2][j+1]
          + A_ij(i,j)*u_ij[i][j+1]
          + F_ij(i,j)*u_ij[i+1][j+2]
          + G_ij(i,j)*u_ij[i+1][j]);
  double shift = norm(u_ij[i+1][j+1] - uij_prev, "inf")
    / std::max(1., norm(uij_prev, "inf"));
  return std::max(worstshift, shift);
}

void Simulation::stamp(ostream& os){
  for (int j=1; j <= J; ++j){
    for (int i=1; i<I; ++i){
      os << u_ij[i][j][1] << ",";
    }
    os << u_ij[I][j][1] << "\n";
  }
}

void Simulation::log(ostream& os){
  os << "ht = " << grid.ht << ", hr = " << grid.hr << ", = " << grid.hz << "\n"
    << "simulation time = " << size.T << ", radius = " << size.R << ", length = " << size.Z << "\n"
    << "MT radius = " << RMT << "\n"
    << "each " << each << "'th step printed to " << path.out << "\n"
    << "SOR omega = " << omg << ", relaxation stop epsilon = " << eps << "\n"
    << "maximum relaxations per timestep: " << maxRelaxationsPerTimestep << "\n"
    << "initial conditions file: " << path.in << "\n"
    << "jacobian mode: " << jmode
    << ", reacions mode: " << rmode << "\n";
  if (path.plots != ""){
    os << "directory with plots: " << path.plots << "\n";
  }
  os.flush();
}

void Simulation::set_init_cond_from_file(const std::string& path){ //JxI.csv
  std::ifstream istream(path);
  std::string line;
  int i=0;
  int j=0;
  int k=0;
  if (istream.is_open()){
    while (getline(istream,line)){
      if (j>=J){throw "J value exceeded\n";}
      std::stringstream ss(line);
      std::vector<std::vector<double>> us;
      std::vector<double> u;
      std::string su;
      std::string _su;
      while (getline(ss, su, ';')){
        if (i>=I){throw "I value exceeded\n";}
        std::stringstream _ss(su);
        while (getline(_ss, _su, ',')){
          if (k>=6){throw "vector of dimension gt 6\n";}
          u.push_back(stod(_su));
          ++k;
        }
        arma::vec _u(u);
        u_ij[i+1][j+1] = _u;
        u.clear();
        k=0;
        ++i;
      }
      i=0;
      ++j;
    }
  } else {
    throw "couldn't open init cond file\n";
  }
}
void Simulation::iniTimeIndep(){
  const double ht = grid.ht, hr = grid.hr, hz = grid.hz;
  for (int i=0; i<I; ++i){
    C_i[i] = (ht*(i + 1.) / (hr*hr*(i + 0.5))) * diffr_coeff(i+1,0);
    A_i[i] = (ht*i / (hr*hr*(i+0.5))) * diffr_coeff(i,0);
  }
  F_i[0] = dzer_mx;
  F_i[1] =(ht / (hz*hz)) * dbulk_mx;
  F_i[2] =(ht / (hz*hz)) * dmt_mx;
  G_i[0] = dzer_mx;
  G_i[1] =(ht / (hz*hz)) * dbulk_mx;
  G_i[2] =(ht / (hz*hz)) * dmt_mx;
  for (int i=0; i<I; ++i){
    for (int j=0; j<J; ++j){
      //BJL_ij[i][j] = E_mx + (ht / (hr*hr*(i+0.5))) * ((i + 1.)*diffr_coeff(i+1,j) + i*diffr_coeff(i,j))
      //         + (ht / (hz*hz)) * (diffz_coeff(i,j+1) + diffz_coeff(i,j));
      BJL_ij[i][j] = E_mx + C_ij(i,j) + A_ij(i,j) + F_ij(i,j) + G_ij(i,j);
    }
  }
}

//1D->2D
bool Simulation::r_in_mt(int i, int j){
  return i < rmti && i > 0;
}
bool Simulation::r_in_bulk(int i, int j){
  return i > rmti && i < I; //was <=
}
//bool Simulation::r_in_av(int i, int j){//is incorrect and not called
//  return i == rmti;
//}
bool Simulation::r_in_vanish(int i, int j){
  return i==0 || i==I; // was I+1
}
const arma::mat::fixed<6,6>& Simulation::diffr_coeff(int i, int j){
  if (r_in_bulk(i,j)) {return dbulk_mx;}
  if (r_in_mt(i,j)) {return dmt_mx;}
  if (r_in_vanish(i,j)) {return dzer_mx;}
  return dav_mx;
} // diffr_coeff(i,j) stands for d_i,j+0.5
bool Simulation::z_in_mt(int i, int j){
  return i < rmti && j > 0 && j < J; //was <=
}
bool Simulation::z_in_bulk(int i, int j){
  return i >= rmti && j > 0 && j < J; //was <=
}
//bool Simulation::z_in_av(int i, int j);
bool Simulation::z_in_vanish(int i, int j){
  return j == 0 || j == J; // was J+1
}
const arma::mat::fixed<6,6>& Simulation::diffz_coeff(int i, int j){//??? 
  if (z_in_bulk(i,j)) {return dbulk_mx;}
  if (z_in_mt(i,j)) {return dmt_mx;}
  if (z_in_vanish(i,j)) {return dzer_mx;}
  return dav_mx;
} // diffz_coeff(i,j) stands for d_i+0.5,j

const arma::mat::fixed<6,6>& Simulation::C_ij(int i, int j){
  return C_i[i];
}
const arma::mat::fixed<6,6>& Simulation::A_ij(int i, int j){
  return A_i[i];
}
const arma::mat::fixed<6,6>& Simulation::F_ij(int i, int j){
  /* question is what return type should be (to & or not to &)?
   * another question is who optimizes better -- me or compiler? */
  if (j < J-1){
    if (i > rmti) {
      return F_i[1];
    }
    return F_i[2];
  }
  return F_i[0];
  //return (j < J-1)*((i >= rmti) ? F_i[1] : F_i[2]) + dzer_mx; //was i>
}
const arma::mat::fixed<6,6>& Simulation::G_ij(int i, int j){
  /* same questions here */
  if (j > 0){
    if (i > rmti){
      return G_i[1];
    }
    return G_i[2];
  }
  return G_i[0];
  //return (j > 0)*((i >= rmti) ? G_i[1] : G_i[2]) + dzer_mx; //was i>
}

void runSimulation(const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int each
      , const PathSpec& ps
      , mode::Jacobian jmode = mode::Jacobian::on
      , mode::Reactions rmode = mode::Reactions::on
      ){
  //Simulation(gs, ss, _omg, _eps, _RMT, each, ps, jmode, rmode).run();
  Simulation sim(gs, ss, _omg, _eps, _RMT, each, ps, jmode, rmode);
  sim.run();

}
void simParametersNotify(ostream& os
      , const GridSpec& gs
      , const SizeSpec& ss
      , double _omg, double _eps
      , double _RMT
      , int each
      , const PathSpec& ps
      , mode::Jacobian jmode = mode::Jacobian::on
      , mode::Reactions rmode = mode::Reactions::on
      ){
  os << "starting simulation\n";
  os << "ht = " << gs.ht << ", hr = " << gs.hr << ", = " << gs.hz << "\n"
    << "simulation time = " << ss.T << ", radius = " << ss.R << ", length = " << ss.Z << "\n"
    << "MT radius = " << _RMT << "\n"
    << "each " << each << "'th step will be printed to " << ps.out << "\n"
    << "SOR omega = " << _omg << ", relaxation stop epsilon = " << _eps << "\n"
    << "initial conditions file: " << ps.in << "\n"
    << "jacobian mode: " << jmode << ", reactions mode: " << rmode << "\n";
  if (ps.plots != ""){
    os << "directory with plots: " << ps.plots << "\n";
  }
  if (ps.log != "")
    os << "log file: " << ps.log << "\n";
  os.flush();
}

const arma::mat::fixed<6,6> E_mx =
           { { 1, 0, 0, 0, 0, 0}
           , { 0, 1, 0, 0, 0, 0}
           , { 0, 0, 1, 0, 0, 0}
           , { 0, 0, 0, 1, 0, 0}
           , { 0, 0, 0, 0, 1, 0}
           , { 0, 0, 0, 0, 0, 1}
};
const arma::mat::fixed<6,6> diffusionMask = 
                 { { 1, 0, 0, 0, 0, 0}
                 , { 0, 1, 0, 0, 0, 0}
                 , { 0, 0, 1, 0, 0, 0}
                 , { 0, 0, 0, 0, 0, 0}
                 , { 0, 0, 0, 0, 0, 0}
                 , { 0, 0, 0, 0, 0, 0}
};
const arma::mat::fixed<6,6> dbulk_mx = dbulk * diffusionMask;
const arma::mat::fixed<6,6> dmt_mx = dmt * diffusionMask;
const arma::mat::fixed<6,6> dav_mx = dav * diffusionMask;
const arma::mat::fixed<6,6> dzer_mx = 0 * diffusionMask;

/* plotting */
void Simulation::draw(int t, int n_len){
  stringstream filename;
  filename << path.plots << "contour" << setfill('0') << setw(n_len) << to_string(t) << ".png";
  Gnuplot gp;
  gp << "set output '" + filename.str() + "'\n";
  Gri2hd<double> aaConcentration(I, J
                                , 0, 0 //xinit, yinit ??
                                , grid.hr, grid.hz //hx, hy??
                                );
  for (int i = 0; i < I; ++i){
    for (int j = 0; j < J; ++j){
      aaConcentration[i][j] = u_ij[i+1][j+1][1];
    }
  }
  vector<string> annotationNames = {"time", "ht", "hr", "hz", "ω", "ε"};//"\\omega", "\\epsilon"};
  vector<double> annotationValues = {t * grid.ht, grid.ht, grid.hr, grid.hz, omg, eps};
  auto ainfoGenerate = [](vector<string> an, vector<double> av){
                               if (an.size() != av.size()) {
                                 cerr << "annotation names and values size mismatch" << endl;
                                 return plottery::AnnotationInfo({std::min(an.size(), av.size()), an, av});
                               } else {
                                 return plottery::AnnotationInfo({an.size(), an, av});
                               }
                             };

  gp << plottery::plotUsual( move(aaConcentration)
                           , ainfoGenerate(move(annotationNames), move(annotationValues))
                           )
                           .str();
  //TODO
}

}
