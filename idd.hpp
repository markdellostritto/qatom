#ifndef IDD_HPP
#define IDD_HPP

//c libraries
#include <cstdlib>
//c++ libraries
#include <iostream>
//Eigen
#include <Eigen/Dense>
//ame
#include "math_const.hpp"

struct IDD{
	// Constants
	static const double erf_const;
	static const int IDD_NUM=4;
	// Form
	enum Form{
		IDEAL=0,
		LINEAR=1,
		EXP=2,
		ERF=3,
		UNKNOWN=-1
	};
	static Form load(const char* str);
	// Interaction matrices
	static Eigen::Matrix3d& itensor_ideal(const Eigen::Vector3d& dr, Eigen::Matrix3d& mat, double a=1.0);
	static Eigen::Matrix3d& itensor_exp(const Eigen::Vector3d& dr, Eigen::Matrix3d& mat, double a=1.0);
	static Eigen::Matrix3d& itensor_linear(const Eigen::Vector3d& dr, Eigen::Matrix3d& mat, double a=1.0);
	static Eigen::Matrix3d& itensor_erf(const Eigen::Vector3d& dr, Eigen::Matrix3d& mat, double a=1.0);
	// Scaling factors
	static double scale_ideal(double alpha1, double alpha2);
	static double scale_exp(double alpha1, double alpha2);
	static double scale_linear(double alpha1, double alpha2);
	static double scale_erf(double alpha1, double alpha2);
	// Functions - Interaction matrices
	static const std::function<Eigen::Matrix3d& (const Eigen::Vector3d& r, Eigen::Matrix3d& mat, double a)> iTensor[IDD_NUM];
	// Functions - Scaling Factors
	static const std::function<double (double a1, double a2)> scale[IDD_NUM];
};

std::ostream& operator<<(std::ostream& out, const IDD::Form& t);

#endif