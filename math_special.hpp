#ifndef MATH_SPECIAL_HPP
#define MATH_SPECIAL_HPP

//c libraries
#include <cmath>
//c++ libraries
#include <iostream>
#include <vector>
//local
#include "math_const.hpp"
#include "math_func.hpp"

#ifndef DEBUG_MATH_SPECIAL
#define DEBUG_MATH_SPECIAL 0
#endif 

namespace special{
	
	//**************************************************************
	//Typedefs
	//**************************************************************
	
	typedef unsigned int uint;
	
	//**************************************************************
	//Constants
	//**************************************************************
	
	static const double prec=1E-8;
	
	//**************************************************************
	//Sign functions
	//**************************************************************

	template <class T> inline int sign(T x){return (x>0)-(x<0);}
	
	//**************************************************************
	//Modulus functions
	//**************************************************************
	
	template <class T> inline T mod(T n, T z){return (n%z+z)%z;}
	template <> inline int mod<int>(int n, int z){return (n%z+z)%z;}
	template <> inline uint mod<uint>(uint n, uint z){return (n%z+z)%z;}
	template <> inline float mod<float>(float n, float z){return fmod(fmod(n,z)+z,z);}
	template <> inline double mod<double>(double n, double z){return fmod(fmod(n,z)+z,z);}
	template <class T> inline T mod(T n, T lLim, T uLim){return mod<T>(n-lLim,uLim-lLim)+lLim;}
		
	//**************************************************************
	//Step functions
	//**************************************************************
	
	template <class T> inline T step(T x, T c){return x>c;}
	template <class T> inline T rect(T x, T b, T e){return (x>b)*(x<e);}
	
	//**************************************************************
	//Sigmoid function
	//**************************************************************
	
	template <class T> inline T sigmoid(T x){return 1.0/(1.0+std::exp(-x));}
	
	//**************************************************************
	//Delta functions
	//**************************************************************
	
	inline int delta(int x1, int x2){return (x1==x2);}
	inline uint delta(uint x1, uint x2){return (x1==x2);}
	inline float delta(float x, float zero=num_const::ZERO){return std::fabs(x)<zero;}
	inline double delta(double x, double zero=num_const::ZERO){return std::fabs(x)<zero;}
	
	//**************************************************************
	//Complementary Error Function - Approximations
	//**************************************************************
	
	struct erfca_const{
		static const double a1[5];
		static const double a2[5];
		static const double a3[7];
		static const double a4[7];
	};
	
	double erfca1(double x);//max error: 5e-4
	double erfca2(double x);//max error: 2.5e-5
	double erfca3(double x);//max error: 3e-7
	double erfca4(double x);//max error: 1.5e-7
	
	//**************************************************************
	//Kummer's (confluent hypergeometric) function 
	//**************************************************************
	
	double M(double a, double b, double z, double prec=1e-8);
	
	//**************************************************************
	//Legendre Poylnomials
	//**************************************************************
	
	std::vector<double>& legendre(unsigned int n, std::vector<double>& c);
	
	//**************************************************************
	//Chebyshev Polynomials
	//**************************************************************
	
	double chebyshev1r(unsigned int n, double x);//chebyshev - first kind  - recursive
	double chebyshev1i(unsigned int n, double x);//chebyshev - first kind  - iterative
	double chebyshev2r(unsigned int n, double x);//chebyshev - second kind - recursive
	double chebyshev2i(unsigned int n, double x);//chebyshev - second kind - iterative
	std::vector<double>& chebyshev1i(unsigned int n, double x, std::vector<double>& r);//chebyshev - first kind  - iterative - all
	std::vector<double>& chebyshev2i(unsigned int n, double x, std::vector<double>& r);//chebyshev - second kind - iterative - all
	std::vector<double>& chebyshev1_root(unsigned int n, std::vector<double>& r);//chebyshev - first kind  - roots
	std::vector<double>& chebyshev2_root(unsigned int n, std::vector<double>& r);//chebyshev - second kind - roots
	
	//**************************************************************
	//Jacobi Polynomials
	//**************************************************************
	double jacobi(unsigned int n, double a, double b, double x);
	std::vector<double>& jacobi(unsigned int n, double a, double b, std::vector<double>& c);
	
	//**************************************************************
	//Laguerre Polynomials
	//**************************************************************
	std::vector<double>& laguerre(unsigned int n, std::vector<double>& c);
}

#endif
