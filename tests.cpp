#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/tools/roots.hpp>

#include "particles.h"

using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

int main(){
    long double T = 1e6;
    long double MS = 100.0;
    long double lHS = 1e-11;

    particle p(91.1876,2,1);

    cout<<p.ColissionTerm(T,MS,lHS)<<endl;

}
