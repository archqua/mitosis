#ifndef BOPTIONS_HPP
#define BOPTIONS_HPP

#include "parameters.hpp"

#include <boost/program_options.hpp>
#include <vector>
#include <string>
//#include <iostream>
#include <iosfwd>

namespace bulk{
namespace options{

extern std::vector<double> hts, hrs, hzs;
extern std::vector<double> Ts, Rs, Zs;
extern std::vector<std::string> ins, outs, logs, plotDirs;
extern std::vector<double> omgs, epss;
extern std::vector<int> eachs;
extern size_t jobs;
extern std::vector<double> rmts;
extern std::vector<mode::Jacobian> jmodes;
extern std::vector<mode::Reactions> rmodes;
extern boost::program_options::variables_map vm;

void set_from_cmd(int argc, char** argv);
void print_help(std::ostream& os);
bool check();

}
}

#endif
