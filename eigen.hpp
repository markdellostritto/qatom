#ifndef EIGEN_HPP
#define EIGEN_HPP

//c libraries
#include <cstdio>
#include <cstdlib>
//c++ libraries
#include <iostream>
#include <Eigen/Dense>
//local libraries
#include "string.hpp"
#include "serialize.hpp"

namespace eigen{

struct LIN_SOLVER{
	enum type{
		LLT=0,//cholesky decomposition
		LDLT=1,//cholesky-variant
		PPLU=2,//LU - partial pivoting
		FPLU=3,//LU - full pivoting
		HQR=4,//Householder QR
		CPHQR=5,//Householder QR - column pivoting
		UNKNOWN=-1
	};
	static LIN_SOLVER::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const LIN_SOLVER::type& t);

Eigen::Vector3d& load(const char* str, Eigen::Vector3d& vec);

const char* print(char* str, const Eigen::Vector3d& vec);

}

namespace serialize{

//**********************************************
// byte measures
//**********************************************

template <> unsigned int nbytes(const Eigen::Vector3d& obj);
template <> unsigned int nbytes(const Eigen::VectorXd& obj);
template <> unsigned int nbytes(const Eigen::Matrix3d& obj);
template <> unsigned int nbytes(const Eigen::MatrixXd& obj);

//**********************************************
// packing
//**********************************************

template <> void pack(const Eigen::Vector3d& obj, char* arr);
template <> void pack(const Eigen::VectorXd& obj, char* arr);
template <> void pack(const Eigen::Matrix3d& obj, char* arr);
template <> void pack(const Eigen::MatrixXd& obj, char* arr);

//**********************************************
// unpacking
//**********************************************

template <> void unpack(Eigen::Vector3d& obj, const char* arr);
template <> void unpack(Eigen::VectorXd& obj, const char* arr);
template <> void unpack(Eigen::Matrix3d& obj, const char* arr);
template <> void unpack(Eigen::MatrixXd& obj, const char* arr);

}

#endif