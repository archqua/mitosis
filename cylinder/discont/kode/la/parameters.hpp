#ifndef BPARAMETERS_HPP
#define BPARAMETERS_HPP

#include "la.hpp"
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
		la::Matrix<6>( const la::Vector<6>&
				     , int , int 
				     )
	>;
using RHS 
	= std::function<
		la::Vector<6>( const la::Vector<6>&
				   , int, int )
	>;
//jacobian
la::Matrix<6> jac0(const la::Vector<6>& u, int i, int j);

//double _knfan456(const la::Vector<6>& u, int i, int j){
	//return knfn0(i,j) - knfal0*(u[3] + u[4] + u[5]);
//}
la::Matrix<6> jac1(const la::Vector<6>& u, int i, int j);

extern const Jacobian jacs[2];

//right hand side (chemical reactions)
la::Vector<6> rhs0(const la::Vector<6>& u, int i, int j);

la::Vector<6> rhs1(const la::Vector<6>& u, int i, int j);

extern const RHS rhss[2];

}

#endif
