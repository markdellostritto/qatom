#ifndef MATH_FUNC_HPP
#define MATH_FUNC_HPP

#include <cmath>
#include <vector>

#ifndef DEBUG_MATH_SPECIAL
#define DEBUG_MATH_SPECIAL 0
#endif 

namespace function{
	
	//**************************************************************
	//polynomial evaluation
	//**************************************************************
	
	double poly(double x, const std::vector<double>& a);
	double poly(double x, const double* a, unsigned int s);
	template <unsigned int N> double poly(double x, const double* a)noexcept{
		unsigned int s=N;
		double result=a[--s];
		while(s>0) result=x*result+a[--s];
		return result;
	}
	template <> inline double poly<0>(double x, const double* a)noexcept{
		return a[0];
	}
	template <> inline double poly<1>(double x, const double* a)noexcept{
		return x*a[1]+a[0];
	}
	template <> inline double poly<2>(double x, const double* a)noexcept{
		return x*(x*a[2]+a[1])+a[0];
	}
	template <> inline double poly<3>(double x, const double* a)noexcept{
		return x*(x*(x*a[3]+a[2])+a[1])+a[0];
	}
	template <> inline double poly<4>(double x, const double* a)noexcept{
		return x*(x*(x*(x*a[4]+a[3])+a[2])+a[1])+a[0];
	}
	template <> inline double poly<5>(double x, const double* a)noexcept{
		return x*(x*(x*(x*(x*a[5]+a[4])+a[3])+a[2])+a[1])+a[0];
	}
	template <> inline double poly<6>(double x, const double* a)noexcept{
		return x*(x*(x*(x*(x*(x*a[6]+a[5])+a[4])+a[3])+a[2])+a[1])+a[0];
	}
	
}

#endif