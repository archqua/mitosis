#! /bin/sh
# compile.sh
g++ --std=c++17 ./*cpp -I/home/tor/gitlibs/gnuplot-iostream -pthread -lboost_program_options -larmadillo -lboost_iostreams -lboost_system -lboost_filesystem -Wall -Werror -O3
