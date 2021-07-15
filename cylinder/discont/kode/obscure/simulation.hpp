#ifndef SIMULATION_HPP
#define SIMULATION_HPP
#include "gri2hd.hpp"
#include "parameters.hpp"
#include <armadillo>
#include <string>

namespace bulk {

const arma::umat unitLocations = { { 0, 1, 2, 3, 4, 5 }
			   	 , { 0, 1, 2, 3, 4, 5 }
}; //2xN
const arma::sp_mat E_mx(unitLocations, arma::vec({1, 1, 1, 1, 1, 1}));// 6x6 identity

const arma::umat diffusionLocations = { { 0, 1, 2 }
				      , { 0, 1, 2 }};
const arma::sp_mat dbulk_mx=arma::sp_mat(diffusionLocations, dbulk * arma::vec({1,1,1}));//(6,6);
const arma::sp_mat dmt_mx=arma::sp_mat(diffusionLocations, dmt * arma::vec({1,1,1}));
const arma::sp_mat dav_mx=arma::sp_mat(diffusionLocations, dav * arma::vec({1,1,1}));
const arma::sp_mat dzer_mx=arma::sp_mat(6, 6);

class Simulation{
	/* parameters 
	 * should be initialized at construction */
	/* 	timestep	gridsteps along R and Z */
	const double ht = dflt::ht, hr = dflt::hr, hz = dflt::hz;
	/* simulation time	cell radius and length */
	const double T = dflt::T, R = dflt::R, Z = dflt::Z;
	/* SOR parameter	relaxation stop criterion */
	const double omg = dflt::omg, eps = dflt::eps;
	/* MT radius */
	const double RMT = dflt::RMT;
	const int each = dflt::each;
	const std::string initCondPath = "";
	const std::string outPath = "";
	const std::string logPath = "";

	/* derivedfrom parameters 
	 * mustn't be initialized explicitly */
	const int timesteps = T / ht, I = R / hr, J = Z / hz; //will probably work fine
	const int rmti = RMT / hr;

	template <std::size_t n>
	using vec = arma::vec::fixed<n>;
	template <std::size_t i, std::size_t j>
	using mat = arma::mat::fixed<i, j>;
	using sp_mat = arma::sp_mat;

	Gri2hd<vec<6>> u_ij; //u_ij[i][j] stands for u_i-0.5,j-0.5
	/* time-independent */
	Gri2hd<sp_mat> BJL_ij;
	/* time-dependent */
	Gri2hd<mat<6,6>> BJF_ij, B_IJ, Binv_IJ;
	/* time-independent */
	std::vector<sp_mat> C_i;
	std::vector<sp_mat> A_i;
	sp_mat F_i[3];
	sp_mat G_i[3];

	void set_init_cond_from_file(const std::string& path);
	void initialize();

	public:
	Simulation(double _ht, double _hr, double _hz
			, double _T, double _R, double _Z
			, double _omg, double _eps
			, double _RMT
			, int each
			, const std::string& in, const std::string& out, const std::string& log);
	void step();
	void run();
};

void runSimulation(double _ht, double _hr, double _hz
			, double _T, double _R, double _Z
			, double _omg, double _eps
			, double _RMT
			, int each
			, const std::string& in, const std::string& out, const std::string& log);

}

#endif
