#ifndef SIGNAL_HPP
#define SIGNAL_HPP

//c libraries
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
//c++ libraries
#include <iostream>
//math
#include "math_const.hpp"

namespace filter{
	
	struct Gauss{
	public:
		double N,c,s;
		Gauss():s(1),c(0),N(1.0/std::sqrt(2*num_const::PI)){};
		Gauss(double s, double c=0.0):s(s),c(c),N(1.0/(s*std::sqrt(2*num_const::PI))){};
		~Gauss(){};
		double operator()(double x)const{return N*std::exp(-0.5*(x-c)*(x-c)/(s*s));}
		static std::vector<double>& gen(double s, std::vector<double>& f);
	};
	
	typedef std::vector<double> ArrayT;
	
	ArrayT& filter(ArrayT& out, const ArrayT& in, const ArrayT& f, bool pad=true);
}

namespace window{
	
	struct WINDOW_FUNC{
		enum type{
			IDENTITY,
			GAUSSIAN,
			KAISERBESSEL,
			BLACKMANHARRIS,
			UNKNOWN
		};
		static WINDOW_FUNC::type read(const char* str);
	};
	
	std::ostream& operator<<(std::ostream& out, const WINDOW_FUNC::type t);
	
	class Functor{
	public:
		unsigned int N;
		Functor():N(1){};
		Functor(unsigned int n):N(n){};
		virtual ~Functor(){};
		virtual double operator()(unsigned int t)const=0;
	};
	
	class Identity: public Functor{
	public:
		double operator()(unsigned int t)const{return 1.0;};
	};
	
	/*class KaiserBessel: public Functor{
	public:
		double alpha;
		KaiserBessel():alpha(1.0){};
		KaiserBessel(unsigned int n, double a):Functor(n),alpha(a){};
		double operator()(unsigned int t)const{
			return boost::math::cyl_bessel_i(0,num_const::PI*alpha*std::sqrt(1.0-(2.0*t/(N-1.0)-1.0)*(2.0*t/(N-1.0)-1.0)))
			/boost::math::cyl_bessel_i(0,num_const::PI*alpha);
		};
	};*/
	
	class BlackmanHarris: public Functor{
	public:
		BlackmanHarris(){};
		BlackmanHarris(unsigned int n):Functor(n){};
		double operator()(unsigned int t)const{
			/*return 0.35875
				-0.48829*std::cos(2.0*num_const::PI*t/(N-1.0))
				+0.14128*std::cos(4.0*num_const::PI*t/(N-1.0))
				-0.01168*std::cos(6.0*num_const::PI*t/(N-1.0));*/
			double p=2.0*num_const::PI/(N-1.0);
			return 0.35875
				-0.48829*std::cos(t*p)
				+0.14128*std::cos(t*p*2.0)
				-0.01168*std::cos(t*p*3.0);
		};
	};
	
	class Gaussian: public Functor{
	public:
		double sigma;
		Gaussian():sigma(1.0){};
		Gaussian(unsigned int n, double s):Functor(n),sigma(s){};
		double operator()(unsigned int t)const{
			double x=(t-0.5*(N-1.0))/(sigma*0.5*(N-1.0));
			return std::exp(-0.5*x*x);
		};
	};
}

#endif