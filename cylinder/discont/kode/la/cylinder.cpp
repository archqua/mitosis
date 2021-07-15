#include "simulation.hpp"
#include "options.hpp"

#include<iostream>
#include<vector>
#include<future>
#include<thread>

using namespace std;
using namespace bulk;

void run();
int main(int argc, char* argv[]){
	ios_base::sync_with_stdio(false);
	options::set_from_cmd(argc, argv);
	if (options::vm.count("help")){
		options::print_help(cout);
		return 0;
	} else if (options::check()){
		run();
		return 0;
	} else {
		cerr << "command line parameters don't make sense" << endl;
		return 3;
	}
}

void run(){
	vector<future<void>> futures;
	using namespace options;
	for (size_t i = 0; i < jobs; ++i){
		GridSpec gs{.ht = hts[i], .hr = hrs[i], .hz = hzs[i]};
		SizeSpec ss{.T = Ts[i], .R = Rs[i], .Z = Zs[i]};
		PathSpec ps{ .in = ins[i], .out = outs[i], .log = logs[i]
                              /* IILE */
                           , .plots = [](const string& dir){
	                             if (dir == ""){
	                               return dir;
	                             } else if (dir.back() != '/'){
	                               return dir + "/";
	                             } else {
	                               return dir;
	                             }
                                   }(plotDirs[i])};
		if (options::vm.count("verbose"))
			simParametersNotify(cerr
					, gs, ss
					, omgs[i], epss[i], rmts[i], eachs[i]
					, ps
					, jmodes[i], rmodes[i]
					);
		if (options::vm.count("seq"))
			runSimulation(gs, ss
				, omgs[i], epss[i], rmts[i], eachs[i]
				, ps
				, jmodes[i], rmodes[i]);
		else
			futures.push_back(async(
				runSimulation
					, gs, ss
					, omgs[i], epss[i], rmts[i], eachs[i]
					, ps
					, jmodes[i], rmodes[i]
			));
		
	}
}
