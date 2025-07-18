#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/tools/roots.hpp>

using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;
using namespace boost::math::tools;

#define vev 246.0L
#define Gammah 4.07e-3L
#define branchRatio 1.0L - 5.792e-1L -6.240e-2L - 2.165e-4L - 2.876e-4L
#define geffIni 106.75L
#define planckMass 2.4e18L
#define T0 2.3482233139345615e-13L

class particle{
    private:
        long double mass;
        int dof;
        int nc;

    public:
        particle(long double m,int g,int n): mass(m),dof(g),nc(n) {};
        long double get_mass(){return mass;}
        long double get_g(){return dof;}
        long double get_nc(){return nc;}

        long double ColissionTerm(long double T, long double MS, long double lHS){
            double m_H = 125.2;


        auto integrand = [=] (long double s){    
            long double sigma;
            switch (dof){
                case 1:  // scalar boson (Higgs)
                if ((s-4*pow(MS,2))/(s-4*pow(mass,2))>0 && (s-4*pow(mass,2))!=0){
                    long double A = pow(dof,-2)*pow(2*lHS,2)/(32*M_PI*s);
                    long double B = sqrt((s-4*pow(MS,2))/(s-4*pow(mass,2)));
                    std::complex<long double> aux (s-pow(mass,2),mass*Gammah);
                    std::complex<long double> C;
                    C = 2*lHS+3*mass*1.0/aux-4*pow(vev,2)*2*lHS/(s-2*pow(mass,2));
                    sigma = A*B*norm(C);
                }
                else{
                    sigma = 0;
                }
                break;
                case 2:  // fermion
                if ((s-4*pow(MS,2))*(s-4*pow(mass,2))>0){
                    long double A = pow(dof,-2)*pow(2*lHS,2)*nc/(16*M_PI*s);
                    long double B = sqrt((s-4*pow(MS,2))*(s-4*pow(mass,2)));
                    long double C = 1/(pow((s-pow(m_H,2)),2)+pow(m_H*Gammah,2));
                    sigma = A*B*C;
                }
                else {
                    sigma = 0;
                }
                break;
                case 3:   // gauge boson
                if ((s-4*pow(MS,2))/(s-4*pow(mass,2))>0){
                    long double A = pow(dof,-2)*pow(2*lHS,2)/(32*M_PI*s);
                    long double B = sqrt((s-4*pow(MS,2))/(s-4*pow(mass,2)));
                    long double C = (pow(s,2)-4*s*pow(mass,2)+12*pow(mass,4))/(pow(s-pow(m_H,2),2)+pow(m_H*Gammah,2));
                    sigma = A*B*C;
                }
                else {
                    sigma = 0;
                }
                break;
            }
            return sigma*(s-4*mass*mass)*sqrt(s)*boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };

        return (dof*dof*T)/(32*M_PI*M_PI)*exp_sinh<long double>().integrate(integrand,
                                            4.0L*mass*mass,
                                            INFINITY);
    }

};