#ifndef ICC_HPP
#define ICC_HPP

//c libraries
#include <cstdlib>
//c++ libraries
#include <iostream>
//Eigen
#include <Eigen/Dense>
//boost
//#include <boost/math/special_functions/expint.hpp>
//ame
#include "math_const.hpp"

#ifndef DEBUG_ICC
#define DEBUG_ICC 0
#endif

struct ICC{
	// Constants
	static const double erf_const;
	static const int ICC_NUM=4;
	// Form
	enum Form{
		IDEAL=0,
		LINEAR=1,
		EXP=2,
		ERF=3,
		UNKNOWN=-1
	};
	static Form read(const char* str);
	// Interaction matrices
	static double itensor_ideal(double dr, double a=1.0);
	static double itensor_exp(double dr, double a=1.0);
	static double itensor_linear(double dr, double a=1.0);
	static double itensor_erf(double dr, double a=1.0);
	// Scaling factors
	static double scale_ideal(double alpha1, double alpha2);
	static double scale_exp(double alpha1, double alpha2);
	static double scale_linear(double alpha1, double alpha2);
	static double scale_erf(double alpha1, double alpha2);
	// Functions - Interaction matrices
	static const std::function<double (double dr, double a)> iTensor[ICC_NUM];
	// Functions - Scaling Factors
	static const std::function<double (double a1, double a2)> scale[ICC_NUM];
};
std::ostream& operator<<(std::ostream& out, const ICC::Form& t);
	
#endif