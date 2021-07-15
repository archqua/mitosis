#ifndef BPARAMETERS_HPP
#define BPARAMETERS_HPP

#include <armadillo>
#include <cmath>
#include <functional>
//#include <iostream>
#include <iosfwd>

namespace bulk {

namespace dflt {
	extern double ht, hr, hz;
	extern double T, R, Z;
	extern double omg, eps;
	extern double RMT;
	extern int each;
}

double fgus(int i, int j, double sg, double mur, double muz);
using Chromosome = std::function<double(int, int)>;
double n0(int i, int j);

extern const double dbulk;
extern const double dmt;
extern const double dav;

namespace mode{
enum class Jacobian{
	off = 0,
	on = 1
};
enum class Reactions{
	off = 0,
	on = 1
};
std::istream& operator>> (std::istream& is, Jacobian& j);
std::ostream& operator<< (std::ostream& os, const Jacobian& j);
std::istream& operator>> (std::istream& is, Reactions& r);
std::ostream& operator<< (std::ostream& os, const Reactions& r);
}

using Jacobian 
	= std::function<
		arma::mat::fixed<6,6>( const arma::vec::fixed<6>&
				     , int , int 
				     )
	>;
using RHS 
	= std::function<
		arma::vec::fixed<6>( const arma::vec::fixed<6>&
				   , int, int )
	>;
//jacobian
arma::mat::fixed<6,6> jac0(const arma::vec::fixed<6>& u, int i, int j);

//double _knfan456(const arma::vec::fixed<6>& u, int i, int j){
	//return knfn0(i,j) - knfal0*(u[3] + u[4] + u[5]);
//}
arma::mat::fixed<6,6> jac1(const arma::vec::fixed<6>& u, int i, int j);

extern const Jacobian jacs[2];

//right hand side (chemical reactions)
arma::vec::fixed<6> rhs0(const arma::vec::fixed<6>& u, int i, int j);

arma::vec::fixed<6> rhs1(const arma::vec::fixed<6>& u, int i, int j);

extern const RHS rhss[2];

}

#endif
