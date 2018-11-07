#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cstring>

namespace Interp{
	
	struct INTERP_METHOD{
		enum type{
			LINEAR,
			AKIMA,
			UNKNOWN
		};
		static type load(const char* str);
	};
	std::ostream& operator<<(std::ostream& out, const INTERP_METHOD::type& t);
	
	template <class T>
	unsigned int index(T datum, const std::vector<T>& data){
		unsigned int u=data.size()-1,l=0,m;
		while(u-l>1){
			m=l+(u-l)/2;
			if(data[l]<=datum && datum<=data[m]) u=m;
			else l=m;
		}
		return l;
	}
	
	template <class T>
	unsigned int index(T datum, const T* data, unsigned int size){
		unsigned int u=size-1,l=0,m;
		while(u-l>1){
			m=l+(u-l)/2;
			if(data[l]<=datum && datum<=data[m]) u=m;
			else l=m;
		}
		return l;
	}
	
	class Data{
	private:
		int size_;
		std::vector<double> x_;
		std::vector<double> y_;
	public:
		//constructors/destructors
		Data():size_(0),x_(0),y_(0){};
		Data(int s):size_(s),x_(s),y_(s){};
		~Data(){};
		
		//access
		int size()const{return size_;};
		const std::vector<double>& x()const{return x_;};
		const std::vector<double>& y()const{return y_;};
		double& x(int i){return x_[i];};
		const double& x(int i)const{return x_[i];};
		double& y(int i){return y_[i];};
		const double& y(int i)const{return y_[i];};
		
		//member functions
		void resize(int s){size_=s; x_.resize(s); y_.resize(s);};
	};
	
	//double interpLinear(double x, double x1, double x2, double y1, double y2);
	
	//double interpLinear(double x, const Data& d);
	
	//double interpAkima(double x, const double* xData, const double* yData);
	
	double interpAkima(double x, const Data& d);
	
	//**************************************************************
	//Linear Interpolation
	//**************************************************************
	
	double interpLinear(double x, const double* X, const double* Y, unsigned int size);
	double interpLinear(double x, const std::vector<double>& X, const std::vector<double>& Y);
	
	//**************************************************************
	//Akima Interpolation
	//**************************************************************
	
	double interpAkima(double x, const double* X, const double* Y, unsigned int size);
	double interpAkima(double x, const std::vector<double>& X, const std::vector<double>& Y);
	
	//**************************************************************
	//Poylnomial Interpolation
	//**************************************************************
	
	double interpPoly(double x, unsigned int order, const double* X, const double* Y, unsigned int size);
	double interpPoly(double x, unsigned int order, const std::vector<double>& X, const std::vector<double>& Y);
	double interpPoly(double x, unsigned int loc, unsigned int order, const double* X, const double* Y, unsigned int size);
	double interpPoly(double x, unsigned int loc, unsigned int order, const std::vector<double>& X, const std::vector<double>& Y);
	
};

#endif