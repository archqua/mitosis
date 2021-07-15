#include "options.hpp"
#include <iostream>

using namespace std;
namespace bulk{
namespace options{

vector<double> hts = {dflt::ht}, hrs = {dflt::hr}, hzs = {dflt::hz};
vector<double> Ts = {dflt::T}, Rs = {dflt::R}, Zs = {dflt::Z};
vector<string> ins, outs, logs = {""}, plotDirs = {""};
vector<double> omgs = {dflt::omg}, epss = {dflt::eps};
vector<int> eachs = {dflt::each};
size_t jobs = 1;
vector<double> rmts = {dflt::RMT};
vector<mode::Jacobian> jmodes = {mode::Jacobian::on};
vector<mode::Reactions> rmodes = {mode::Reactions::on};
boost::program_options::variables_map vm;
boost::program_options::options_description desc("valid options");
template <class T>
void append(vector<T>& paramList, T value, size_t finalSize){
	while (paramList.size() < finalSize){
		paramList.push_back(value);
	}
}
bool check(){
	if (ins.size() < 1) {
		cerr << "yout must specify at least one file with initial conditions" << endl;
		return false;
	} else if (outs.size() < 1) {
		cerr << "yout must specify at least one ouput file path" << endl;
		return false;
	} else if (outs.size() < jobs) {
		cerr << outs.size() << " output files is less than the numbef of simulations: " << jobs 
			<< "\n setting the number of simulations to " << outs.size() << endl;
		jobs = outs.size();
	} else if (outs.size() > jobs){
		cerr << outs.size() << " output files is greater than the number of simulations: " << jobs
			<< "\n setting the number of simulations to " << outs.size() << endl;
		jobs = outs.size();
	}

	append(hts, hts.back(), jobs);
	append(hrs, hrs.back(), jobs);
	append(hzs, hzs.back(), jobs);
	append(Ts, Ts.back(), jobs);
	append(Rs, Rs.back(), jobs);
	append(Zs, Zs.back(), jobs);
	append(ins, ins.back(), jobs);
	/* outs are nor appended for obvious reasons */
	append(logs, string(""), jobs);
	append(plotDirs, string(""), jobs);
	append(omgs, omgs.back(), jobs);
	append(epss, epss.back(), jobs);
	append(eachs, eachs.back(), jobs);
	append(rmts, rmts.back(), jobs);
	append(jmodes, jmodes.back(), jobs);
	append(rmodes, rmodes.back(), jobs);

        /* IILE */
	cerr << [](size_t jobs) {
          return "ready to launch " + to_string(jobs) + " simulations\n"
		 " if n'th simulation wasn't explicitly given a particular parameter,\n"
		 "  it will borrow this parameter from last simulation that has this parameter\n"
		 "  or default it (log file path and directory with plots always default to nothing)\n"
		 " if a parameter was specified more than " + to_string(jobs) + " times.\n"
		 "  extra occurences will be ignored\n";}(jobs);
        cerr.flush();
	return true;
}

void set_from_cmd(int argc, char** argv){
	namespace po = boost::program_options;
	desc.add_options()
		("help,h"
		 	, "show help msg")
		("verbose,v"
		 	, "print info about each simulation to command line")
		("ht",		po::value<vector<double>>(&hts)
			, "set timestep(s)")
		("hr",		po::value<vector<double>>(&hrs)
		 	, "set gridstep(s) along r")
		("hz",		po::value<vector<double>>(&hzs)
		 	, "set gridstep(s) along z")
		("each,e",	po::value<vector<int>>(&eachs)
		 	, "set stamp-to-file period(s)"
		 	  " in terms of timestep count")
		("time,T",	po::value<vector<double>>(&Ts)
		 	, "set simulation time(s)"
		 	  " (not total timestep count!)")
		("radius,R",	po::value<vector<double>>(&Rs)
		 	, "set cell radius(es)")
		("length,Z",	po::value<vector<double>>(&Zs)
		 	, "set cell length(s)")
		("input,in,if,ic,i", po::value<vector<string>>(&ins)
		 	, "set initial conditions file(s)")
		("output,out,of,o", po::value<vector<string>>(&outs)
		 	, "set output file(s)")
                ("plots,p",      po::value<vector<string>>(&plotDirs)
                        , "set directory where plots go")
		("log",		po::value<vector<string>>(&logs)
		 	, "set log file(s)")
		("omega,omg,w", po::value<vector<double>>(&omgs)
		 	, "set omega SOR parameter(s)")
		("eps",		po::value<vector<double>>(&epss)
		 	, "set relaxation stop epsilon(s)")
		("rmt",		po::value<vector<double>>(&rmts)
		 	, "set MT radius(es)")
		("jmode",	po::value<vector<mode::Jacobian>>(&jmodes)
		 	, "set Jacobian mode(s) (on/off ~ 1/0)")
		("rmode",	po::value<vector<mode::Reactions>>(&rmodes)
		 	, "set Reactions mode(s) (on/off ~ 1/0)")
		("jobs,j",	po::value<size_t>(&jobs)
		 	, "set number of simulaitons"
			  ", pleas keep it equal to number of"
			  " copies of each parameter"
			  " and under the number of your CPU's threads")
		("seq"
		 	,"run simulations sequentially (not in parallel)")
		;
	//po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	
}
void print_help(ostream& os){
	os << desc << endl;
}

}
}
