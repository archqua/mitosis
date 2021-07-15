CXXFLAGS = -O3 --std=c++17 -pthread -Wall -Werror
INCLUDES = -I/home/tor/gitlibs/gnuplot-iostream
LIBS = -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -larmadillo
OBJS = parameters.o options.o simulation.o plottery.o cylinder.o
TOBJS = parameters.o options.o simulation.o plottery.o testsim.o
MAIN = a.out
TMAIN = t.out

COMMAND = g++ $(CXXFLAGS) $(INCLUDES) $(LIBS)

all: $(MAIN)

$(MAIN): $(OBJS)
	$(COMMAND) $(OBJS) -o $(MAIN)
parameters.o: parameters.cpp parameters.hpp la.hpp
	$(COMMAND) parameters.cpp -c

options.o: options.cpp options.hpp parameters.hpp la.hpp
	$(COMMAND) options.cpp -c

plottery.o: plottery.cpp plottery.hpp gri2hd.hpp
	$(COMMAND) plottery.cpp -c

simulation.o: simulation.cpp simulation.hpp gri2hd.hpp parameters.hpp plottery.hpp la.hpp
	$(COMMAND) simulation.cpp -c

cylinder.o: cylinder.cpp simulation.hpp gri2hd.hpp parameters.hpp plottery.hpp options.hpp la.hpp
	$(COMMAND) cylinder.cpp -c
