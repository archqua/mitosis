#include <vector>
#include <tuple>

#include "gnuplot-iostream.h"

using namespace std;

vector<tuple<double, double>>
  vt = { {0, 0}
       , {1, 1}
       , {2, 1}
       , {3, 0}
       };
tuple<vector<double>, vector<double>>
  tv = { {0, 1, 2, 3}
       , {0, 1, 1, 0}
       };

int main() {
  Gnuplot gp;
  gp << "set term png size 1920,1080\n";
  gp << "set output vt.png\n";
  gp << "plot '-' with lines\n";
  gp.send1d(vt);
}
