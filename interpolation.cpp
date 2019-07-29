#include "interpolation.hpp"

namespace Interp{
	
	INTERP_METHOD::type INTERP_METHOD::read(const char* str){
		if(std::strcmp(str,"LINEAR")==0) return INTERP_METHOD::LINEAR;
		else if(std::strcmp(str,"AKIMA")==0) return INTERP_METHOD::AKIMA;
		else return INTERP_METHOD::UNKNOWN;
	}
	
	std::ostream& operator<<(std::ostream& out, const INTERP_METHOD::type& t){
		if(t==INTERP_METHOD::LINEAR) out<<"LINEAR";
		else if(t==INTERP_METHOD::AKIMA) out<<"AKIMA";
		else out<<"UKNOWN";
		return out;
	}
	
	double interpAkima(double x, const Data& D){
		int i=index(x,D.x());
		
		/*
			find the slopes at i, i+1
		*/
		double t1,t2;
		if(i==D.size()-2){
			t2=(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i));
			t1=0.5*((D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i))+(D.y(i)-D.y(i-1))/(D.x(i)-D.x(i-1)));
		} else if(i==D.size()-3){
			t2=0.5*((D.y(i+2)-D.y(i+1))/(D.x(i+2)-D.x(i+1))+(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i)));
			double m0=(D.y(i-1)-D.y(i-2))/(D.x(i-1)-D.x(i-2));
			double m1=(D.y(i)-D.y(i-1))/(D.x(i)-D.x(i-1));
			double m2=(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i));
			double m3=(D.y(i+2)-D.y(i+1))/(D.x(i+2)-D.x(i+1));
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
		} else if(i==0){
			t1=(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i));
			t2=0.5*((D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i))+(D.y(i+2)-D.y(i+1))/(D.x(i+2)-D.x(i+1)));
		} else if(i==1){
			t1=0.5*((D.y(i)-D.y(i-1))/(D.x(i)-D.x(i-1))+(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i)));
			double m1=(D.y(i)-D.y(i-1))/(D.x(i)-D.x(i-1));
			double m2=(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i));
			double m3=(D.y(i+2)-D.y(i+1))/(D.x(i+2)-D.x(i+1));
			double m4=(D.y(i+3)-D.y(i+2))/(D.x(i+3)-D.x(i+2));
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		} else {
			double m0=(D.y(i-1)-D.y(i-2))/(D.x(i-1)-D.x(i-2));
			double m1=(D.y(i)-D.y(i-1))/(D.x(i)-D.x(i-1));
			double m2=(D.y(i+1)-D.y(i))/(D.x(i+1)-D.x(i));
			double m3=(D.y(i+2)-D.y(i+1))/(D.x(i+2)-D.x(i+1));
			double m4=(D.y(i+3)-D.y(i+2))/(D.x(i+3)-D.x(i+2));
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		}
		
		/*
			calculate the fitting coefficients
		*/
		double deltax=D.x(i+1)-D.x(i);
		double a=(t2+t1)/(deltax*deltax)-2*(D.y(i+1)-D.y(i))/(deltax*deltax*deltax);
		double b=-(t2+2*t1)/deltax+3*(D.y(i+1)-D.y(i))/(deltax*deltax);
		double c=t1;
		double d=D.y(i);
		
		/*
			calculate the interpolated value
		*/
		return a*std::pow(x-D.x(i),3)+b*std::pow(x-D.x(i),2)+c*(x-D.x(i))+d;
	}
	
	//**************************************************************
	//Linear Interpolation
	//**************************************************************
	
	double interpLinear(double x, const double* X, const double* Y, unsigned int size){
		int i=index(x,X,size);
		return Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*x;
	}
	
	double interpLinear(double x, const std::vector<double>& X, const std::vector<double>& Y){
		int i=index(x,X);
		return Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x-X[i]);
	}
	
	//**************************************************************
	//Akima Interpolation
	//**************************************************************
	
	double interpAkima(double x, const double* X, const double* Y, unsigned int size){
		int i=index(x,X,size);
		
		/*
			find the slopes at i, i+1
		*/
		double t1,t2;
		if(i==size-2){
			t2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			t1=0.5*((Y[i+1]-Y[i])/(X[i+1]-X[i])+(Y[i]-Y[i-1])/(X[i]-X[i-1]));
		} else if(i==size-3){
			t2=0.5*((Y[i+2]-Y[i+1])/(X[i+2]-X[i+1])+(Y[i+1]-Y[i])/(X[i+1]-X[i]));
			double m0=(Y[i-1]-Y[i-2])/(X[i-1]-X[i-2]);
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
		} else if(i==0){
			t1=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			t2=0.5*((Y[i+1]-Y[i])/(X[i+1]-X[i])+(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]));
		} else if(i==1){
			t1=0.5*((Y[i]-Y[i-1])/(X[i]-X[i-1])+(Y[i+1]-Y[i])/(X[i+1]-X[i]));
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			double m4=(Y[i+3]-Y[i+2])/(X[i+3]-X[i+2]);
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		} else {
			double m0=(Y[i-1]-Y[i-2])/(X[i-1]-X[i-2]);
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			double m4=(Y[i+3]-Y[i+2])/(X[i+3]-X[i+2]);
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		}
		
		/*
			calculate the fitting coefficients
		*/
		double deltax=X[i+1]-X[i];
		double a=(t2+t1)/(deltax*deltax)-2*(Y[i+1]-Y[i])/(deltax*deltax*deltax);
		double b=-(t2+2*t1)/deltax+3*(Y[i+1]-Y[i])/(deltax*deltax);
		double c=t1;
		double d=Y[i];
		
		/*
			calculate the interpolated value
		*/
		return a*(x-X[i])*(x-X[i])*(x-X[i])+b*(x-X[i])*(x-X[i])+c*(x-X[i])+d;
	}
	
	double interpAkima(double x, const std::vector<double>& X, const std::vector<double>& Y){
		int i=index(x,X);
		
		/*
			find the slopes at i, i+1
		*/
		double t1,t2;
		if(i==X.size()-2){
			t2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			t1=0.5*((Y[i+1]-Y[i])/(X[i+1]-X[i])+(Y[i]-Y[i-1])/(X[i]-X[i-1]));
		} else if(i==X.size()-3){
			t2=0.5*((Y[i+2]-Y[i+1])/(X[i+2]-X[i+1])+(Y[i+1]-Y[i])/(X[i+1]-X[i]));
			double m0=(Y[i-1]-Y[i-2])/(X[i-1]-X[i-2]);
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
		} else if(i==0){
			t1=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			t2=0.5*((Y[i+1]-Y[i])/(X[i+1]-X[i])+(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]));
		} else if(i==1){
			t1=0.5*((Y[i]-Y[i-1])/(X[i]-X[i-1])+(Y[i+1]-Y[i])/(X[i+1]-X[i]));
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			double m4=(Y[i+3]-Y[i+2])/(X[i+3]-X[i+2]);
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		} else {
			double m0=(Y[i-1]-Y[i-2])/(X[i-1]-X[i-2]);
			double m1=(Y[i]-Y[i-1])/(X[i]-X[i-1]);
			double m2=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
			double m3=(Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]);
			double m4=(Y[i+3]-Y[i+2])/(X[i+3]-X[i+2]);
			if(((m0==m1)&&((m2==m3)||(m0==m2)))||((m1==m3)&&((m2==m3)||(m0==m2)))){
				t1=0.5*m1+0.5*m2;
			} else {
				t1=(std::sqrt(std::abs((m2-m0)*(m3-m2)))*m1+std::sqrt(std::abs((m1-m0)*(m3-m1)))*m2)
				/(std::sqrt(std::abs((m2-m0)*(m3-m2)))+std::sqrt(std::abs((m1-m0)*(m3-m1))));
			}
			if(((m1==m2)&&((m3==m4)||(m1==m3)))||((m2==m4)&&((m3==m4)||(m1==m3)))){
				t2=0.5*m2+0.5*m3;
			} else {
				t2=(std::sqrt(std::abs((m3-m1)*(m4-m3)))*m2+std::sqrt(std::abs((m2-m1)*(m4-m2)))*m3)
				/(std::sqrt(std::abs((m3-m1)*(m4-m3)))+std::sqrt(std::abs((m2-m1)*(m4-m2))));
			}
		}
		
		/*
			calculate the fitting coefficients
		*/
		double deltax=X[i+1]-X[i];
		double a=(t2+t1)/(deltax*deltax)-2*(Y[i+1]-Y[i])/(deltax*deltax*deltax);
		double b=-(t2+2*t1)/deltax+3*(Y[i+1]-Y[i])/(deltax*deltax);
		double c=t1;
		double d=Y[i];
		
		/*
			calculate the interpolated value
		*/
		return a*(x-X[i])*(x-X[i])*(x-X[i])+b*(x-X[i])*(x-X[i])+c*(x-X[i])+d;
	}
	
	//**************************************************************
	//Poylnomial Interpolation
	//**************************************************************
	
	double interpPoly(double x, unsigned int order, const double* X, const double* Y, unsigned int size){
		int ii=index(x,X,size);
		unsigned int orderL=order+1;//the "logical" order is one more than the order of approximation
		if(orderL>size) throw std::invalid_argument("Too few points for polynomial interpolation");
		std::vector<double> c(orderL);
		std::vector<double> d(orderL);
		
		int origin;
		if(ii-size/2<0) origin=0;
		else if(ii+size/2>size-1) origin=size-orderL;
		else origin=ii-size/2;
		
		//initialize the poylnomial coefficients (orderL=0 terms)
		for(int i=0; i<orderL; i++){
			c[i]=Y[i]; d[i]=Y[i];
		}
		double y=Y[ii--];//initial approximation
		for(int m=1; m<orderL; m++){
			for(int i=0; i<orderL-m; i++){
				double w=c[i+1]-d[i];
				d[i]=(X[i+m]-x)*w/(X[i]-X[i+m]);
				c[i]=(X[i]-x)*w/(X[i]-X[i+m]);
			}
			y+=(2*(ii+1)<(orderL-m)?c[ii+1]:d[ii--]);
		}
		
		return y;
	}
	
	double interpPoly(double x, unsigned int order, const std::vector<double>& X, const std::vector<double>& Y){
		int ii=index(x,X);
		unsigned int orderL=order+1;//the "logical" order is one more than the order of approximation
		if(orderL>X.size()) throw std::invalid_argument("Too few points for polynomial interpolation");
		std::vector<double> c(orderL);
		std::vector<double> d(orderL);
		
		int origin;
		if(ii-X.size()/2<0) origin=0;
		else if(ii+X.size()/2>X.size()-1) origin=X.size()-orderL;
		else origin=ii-X.size()/2;
		
		//initialize the poylnomial coefficients (orderL=0 terms)
		for(int i=0; i<orderL; i++){
			c[i]=Y[i]; d[i]=Y[i];
		}
		double y=Y[ii--];//initial approximation
		for(int m=1; m<orderL; m++){
			for(int i=0; i<orderL-m; i++){
				double w=c[i+1]-d[i];
				d[i]=(X[i+m]-x)*w/(X[i]-X[i+m]);
				c[i]=(X[i]-x)*w/(X[i]-X[i+m]);
			}
			y+=(2*(ii+1)<(orderL-m)?c[ii+1]:d[ii--]);
		}
		
		return y;
	}
	
	double interpPoly(double x, unsigned int loc, unsigned int order, const double* X, const double* Y, unsigned int size){
		unsigned int orderL=order+1;//the polynomial order is one more than the order of approximation
		if(orderL>size) throw std::invalid_argument("Too few points for polynomial interpolation");
		std::vector<double> c(orderL);
		std::vector<double> d(orderL);
		
		int origin=loc;
		double diff=std::fabs(x-X[origin]);
		unsigned int ii=0;
		for(int i=0; i<orderL; i++){
			if(std::fabs(x-X[origin+i])<diff){
				ii=i;
				diff=std::fabs(x-X[origin+i]);
			}
		}
		
		//initialize the poylnomial coefficients (orderL=0 terms)
		for(int i=0; i<orderL; i++){
			c[i]=Y[origin+i]; d[i]=Y[origin+i];
		}
		double y=Y[origin+ii--];//initial approximation
		for(int m=1; m<orderL; m++){
			for(int i=0; i<orderL-m; i++){
				double w=c[i+1]-d[i];
				d[i]=(X[origin+i+m]-x)*w/(X[origin+i]-X[origin+i+m]);
				c[i]=(X[origin+i]-x)*w/(X[origin+i]-X[origin+i+m]);
			}
			y+=(2*(ii+1)<(orderL-m)?c[ii+1]:d[ii--]);
		}
		
		return y;
	}
	
	double interpPoly(double x, unsigned int loc, unsigned int order, const std::vector<double>& X, const std::vector<double>& Y){
		unsigned int orderL=order+1;//the polynomial order is one more than the order of approximation
		if(orderL>X.size()) throw std::invalid_argument("Too few points for polynomial interpolation");
		std::vector<double> c(orderL);
		std::vector<double> d(orderL);
		
		int origin=loc;
		double diff=std::fabs(x-X[origin]);
		unsigned int ii=0;
		for(int i=0; i<orderL; i++){
			if(std::fabs(x-X[origin+i])<diff){
				ii=i;
				diff=std::fabs(x-X[origin+i]);
			}
		}
		
		//initialize the poylnomial coefficients (orderL=0 terms)
		for(int i=0; i<orderL; i++){
			c[i]=Y[origin+i]; d[i]=Y[origin+i];
		}
		double y=Y[origin+ii--];//initial approximation
		for(int m=1; m<orderL; m++){
			for(int i=0; i<orderL-m; i++){
				double w=c[i+1]-d[i];
				d[i]=(X[origin+i+m]-x)*w/(X[origin+i]-X[origin+i+m]);
				c[i]=(X[origin+i]-x)*w/(X[origin+i]-X[origin+i+m]);
			}
			y+=(2*(ii+1)<(orderL-m)?c[ii+1]:d[ii--]);
		}
		
		return y;
	}
	
}