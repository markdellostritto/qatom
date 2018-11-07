#include "idd.hpp"
	
//*******************************************************
// Constants
//*******************************************************

const double IDD::erf_const=std::pow(2.0/(9.0*num_const::PI),1.0/6.0);

//*******************************************************
// Form
//*******************************************************

IDD::Form IDD::load(const char* str){
	if(std::strcmp(str,"IDEAL")==0) return IDD::Form::IDEAL;
	else if(std::strcmp(str,"LINEAR")==0) return IDD::Form::LINEAR;
	else if(std::strcmp(str,"EXP")==0) return IDD::Form::EXP;
	else if(std::strcmp(str,"ERF")==0) return IDD::Form::ERF;
	else return IDD::Form::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const IDD::Form& t){
	if(t==IDD::Form::IDEAL) out<<"IDEAL";
	else if(t==IDD::Form::LINEAR) out<<"LINEAR";
	else if(t==IDD::Form::EXP) out<<"EXP";
	else if(t==IDD::Form::ERF) out<<"ERF";
	else out<<"UNKNOWN";
	return out;
}

//*******************************************************
// Interaction matrices
//*******************************************************

Eigen::Matrix3d& IDD::itensor_ideal(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a){
	double dr=r.norm();
	//mat.noalias()=(3.0*r*r.transpose()-Eigen::Matrix3d::Identity()*dr*dr)/(dr*dr*dr*dr*dr);
	mat.noalias()=3.0*r*r.transpose();
	mat.noalias()-=Eigen::Matrix3d::Identity()*dr*dr;
	mat*=1.0/(dr*dr*dr*dr*dr);
	return mat;	
}

Eigen::Matrix3d& IDD::itensor_exp(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a){
	double dr=r.norm();
	double b=a*dr;
	mat.noalias()=r*r.transpose()*3*(1.0-((1.0/6.0)*b*b*b+0.5*b*b+b+1)*std::exp(-b))/(dr*dr*dr*dr*dr)
		-Eigen::Matrix3d::Identity()*(1.0-(0.5*b*b+b+1)*std::exp(-b))/(dr*dr*dr);
	return mat;
}

Eigen::Matrix3d& IDD::itensor_linear(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a){
	double dr=r.norm();
	double c5,c3;
	if(dr<a){c5=3.0*std::pow(dr/a,4.0); c3=4.0*std::pow(dr/a,3.0)-c5;}
	else {c5=3.0; c3=1.0;};
	mat.noalias()=(c5*r*r.transpose()-c3*dr*dr*Eigen::Matrix3d::Identity())/(dr*dr*dr*dr*dr);
	return mat;
};

Eigen::Matrix3d& IDD::itensor_erf(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a){
	double dr=r.norm();
	double exp=2.0/num_const::RadPI*a*std::exp(-dr*dr*a*a);
	mat.noalias()=3.0*r*r.transpose();
	mat.noalias()-=Eigen::Matrix3d::Identity()*dr*dr;
	mat*=(std::erf(dr*a)-dr*exp)/(dr*dr*dr*dr*dr);
	mat.noalias()-=2.0*a*a/(dr*dr)*exp*r*r.transpose();
	return mat;
}

//*******************************************************
// Scaling factors
//*******************************************************

double IDD::scale_ideal(double a1, double a2){
	return 1.0;
}

double IDD::scale_exp(double a1, double a2){
	return std::pow(a1*a2,-1.0/6.0);
}

double IDD::scale_linear(double a1, double a2){
	return std::pow(a1*a2,-1.0/6.0);
}

double IDD::scale_erf(double a1, double a2){
	return 1.0/(erf_const*std::sqrt(std::pow(a1,2.0/3.0)+std::pow(a2,2.0/3.0)));
}

//*******************************************************
// Function - Scaling factors
//*******************************************************

const std::function<Eigen::Matrix3d& (const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a)> IDD::iTensor[IDD_NUM]={
	itensor_ideal,
	itensor_linear,
	itensor_exp,
	itensor_erf
};

//*******************************************************
// Function - Scaling factors
//*******************************************************

const std::function<double (double a1, double a2)> IDD::scale[IDD_NUM]={
	scale_ideal,
	scale_linear,
	scale_exp,
	scale_erf
};
