#include "icc.hpp"

//*******************************************************
// Constants
//*******************************************************

const double ICC::erf_const=std::sqrt(2.0/num_const::PI);

//*******************************************************
// Form
//*******************************************************
	
ICC::Form ICC::read(const char* str){
	//std::cout<<"str = "<<str<<"\n";
	if(std::strcmp(str,"IDEAL")==0) return ICC::Form::IDEAL;
	else if(std::strcmp(str,"LINEAR")==0) return ICC::Form::LINEAR;
	else if(std::strcmp(str,"EXP")==0) return ICC::Form::EXP;
	else if(std::strcmp(str,"ERF")==0) return ICC::Form::ERF;
	else return ICC::Form::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const ICC::Form& t){
	if(t==ICC::Form::IDEAL) out<<"IDEAL";
	else if(t==ICC::Form::LINEAR) out<<"LINEAR";
	else if(t==ICC::Form::EXP) out<<"EXP";
	else if(t==ICC::Form::ERF) out<<"ERF";
	else out<<"UNKNOWN";
	return out;
}

//*******************************************************
// Interaction matrices
//*******************************************************

double ICC::itensor_ideal(double dr, double a){
	if(DEBUG_ICC>0) std::cout<<"ICC::itensor_ideal(const Eigen::Vector3d&,double):\n";
	return 1.0/dr;
}

double ICC::itensor_exp(double dr, double a){
	if(DEBUG_ICC>0) std::cout<<"ICC::itensor_exp(const Eigen::Vector3d&,double):\n";
	//return boost::math::expint(-a*dr)-0.5*std::exp(-a*dr)*(3.0+a*dr)-std::log(dr);
	return 0;
}

double ICC::itensor_linear(double dr, double a){
	if(DEBUG_ICC>0) std::cout<<"ICC::itensor_linear(const Eigen::Vector3d&,double):\n";
	return 0.0;//NOT YET IMPLEMENTED
}

double ICC::itensor_erf(double dr, double a){
	if(DEBUG_ICC>0) std::cout<<"ICC::itensor_erf(const Eigen::Vector3d&,double):\n";
	return std::erf(dr*a)/dr;
}

//*******************************************************
// Scaling factors
//*******************************************************

double ICC::scale_ideal(double a1, double a2){
	if(DEBUG_ICC>0) std::cout<<"ICC::scale_ideal(const Eigen::Vector3d&,double):\n";
	return 1.0;
}

double ICC::scale_exp(double c1, double c2){
	if(DEBUG_ICC>0) std::cout<<"ICC::scale_exp(const Eigen::Vector3d&,double):\n";
	return 1.0/(std::sqrt(c1*c1+c2*c2));
}

double ICC::scale_linear(double c1, double c2){
	if(DEBUG_ICC>0) std::cout<<"ICC::scale_linear(const Eigen::Vector3d&,double):\n";
	return 1.0;
}

double ICC::scale_erf(double c1, double c2){
	if(DEBUG_ICC>0) std::cout<<"ICC::scale_erf(const Eigen::Vector3d&,double):\n";
	return 1.0/(std::sqrt(2.0)*erf_const*std::sqrt(c1*c1+c2*c2));
}

//*******************************************************
// Function - Interaction Tensors
//*******************************************************

const std::function<double (double dr, double a)> ICC::iTensor[ICC_NUM]={
	itensor_ideal,
	itensor_linear,
	itensor_exp,
	itensor_erf
};

//*******************************************************
// Function - Scaling factors
//*******************************************************

const std::function<double (double a1, double a2)> ICC::scale[ICC_NUM]={
	scale_ideal,
	scale_linear,
	scale_exp,
	scale_erf
};
