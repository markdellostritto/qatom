#ifndef IR_QEQ_HPP
#define IR_QEQ_HPP

//Thole model, bonds

//c libraries
#include <cstdio>
#include <cstdlib>
#include <ctime>
//c++ libraries
#include <iostream>
#include <fstream>
#include <vector>
//eigen libraries
#include <Eigen/Dense>
#include <Eigen/StdVector>
//fftw
#include <fftw3.h>
//simulation
#include "structure.hpp"
//reading
#include "vasp.hpp"
#include "lammps.hpp"
#include "qe.hpp"
#include "string.hpp"
//signal analysis
#include "signal.hpp"
#include "fft.hpp"
//math
#include "math_const.hpp"
#include "math_gradient.hpp"
#include "interpolation.hpp"
//electrostatics
#include "ewald3D.hpp"
#include "qeq.hpp"

#ifndef DEBUG_IR_QEQ
#define DEBUG_IR_QEQ 0
#endif

//***********************************************************************************************************************************
//Profile
//***********************************************************************************************************************************

class Profile{
private:
	double sigma_;
	double b1_,b2_;
public:
	//constructors/destructors
	Profile():sigma_(0),b1_(0),b2_(0){};
	Profile(double s, double b1, double b2):sigma_(s),b1_(b1),b2_(b2){};
	~Profile(){};
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const Profile& profile);
	double operator()(double x){
		if(x<0.5*(b1_+b2_)) return 0.5*std::erfc(-(x-b1_)*sigma_);
		else return 0.5*std::erfc((x-b2_)*sigma_);
	}
	
	//access
	double& sigma(){return sigma_;};
	const double& sigma()const{return sigma_;};
	double& b1(){return b1_;};
	const double& b1()const{return b1_;};
	double& b2(){return b2_;};
	const double& b2()const{return b2_;};
};

//***********************************************************************************************************************************
//Complex class
//***********************************************************************************************************************************

class Complex{
private:
	double data[2];
public:
	Complex(){};
	Complex(double r, double i){data[0]=r; data[1]=i;};
	Complex(const Complex& c){data[0]=c[0]; data[1]=c[1];};
	~Complex(){};
	
	Complex& operator=(const Complex& c){data[0]=c[0]; data[1]=c[1];};
	double& operator[](int i){return data[i];};
	const double& operator[](int i)const{return data[i];};
	friend std::ostream& operator<<(std::ostream& out, const Complex& c){return out<<"("<<c.data[0]<<","<<c.data[1]<<")";};
};

//***********************************************************************************************************************************
//IRQEQ
//***********************************************************************************************************************************

class IRQEQ{
private:
	//FFT parameters
		fourier::FreqUnit::type freqUnit_;//the unit of the frequency
		double freqCut_;//the cutoff for printing the frequency
		unsigned int freqN_;//the integer cutoff for printing the frequency 
		double minFreq_;//the minimum frequency possible, given the number of timesteps
		double maxFreq_;//the maximum frequency possible, given the number of timesteps
		double freqRes_;//the resolution in frequency space in THz
		unsigned int fourierN_;//the total number of data points sampled
	//stride
		unsigned int strideChg_;//the number of timesteps between sampled data points
		unsigned int nStepsChg_;//the total number of calculated steps for charge
	//FFT windowing
		double sigma_;
		std::function<double (int)> window_;
		window::WINDOW_FUNC::type windowType_;
	//temp
		double T_;
	//profile
		Profile profile_;
	//file i/o
		std::string fileSpectrum_;//spectrum file
		std::string fileDipoleT_;//dipole total file
		std::string fileChgT_;//chg total file
		std::string fileChgA_;//chg atom file
	//calculation flags
		bool calcChg_;//whether we will calculate the charges
		bool calcSpectrum_;//whether we will calculate the spectrum
		bool normalize_;//whether to normalize the spectrum
	//i/o flags
		bool writeDipoleT_;
		bool writeChgT_;
		bool writeChgA_;
		bool readChgA_;
	//ir spectrum
		std::vector<double> irSpectrum_;
	//subset
		std::string subsetstr_;
		std::vector<unsigned int> subset_;
public:
	//constants
	static const double mevPerThz;
	static const double cmiPerThz;
	
	//costructors/destructors
	IRQEQ(){defaults();};
	~IRQEQ(){};
	
	//operators
	IRQEQ& operator=(const IRQEQ& ir3d);
	friend std::ostream& operator<<(std::ostream& out, const IRQEQ& ir3d);
	
	//access
	//calculation flags
		bool& calcChg(){return calcChg_;};
		const bool& calcChg()const{return calcChg_;};
		bool& calcSpectrum(){return calcSpectrum_;};
		const bool& calcSpectrum()const{return calcSpectrum_;};
		bool& normalize(){return normalize_;};
		const bool& normalize()const{return normalize_;};
	//i/o flags
		bool& writeDipoleT(){return writeDipoleT_;};
		const bool& writeDipoleT()const{return writeDipoleT_;};
		bool& writeChgT(){return writeChgT_;};
		const bool& writeChgT()const{return writeChgT_;};
		bool& writeChgA(){return writeChgA_;};
		const bool& writeChgA()const{return writeChgA_;};
		bool& readChgA(){return readChgA_;};
		const bool& readChgA()const{return readChgA_;};
	//file i/o
		std::string& fileSpectrum(){return fileSpectrum_;};
		const std::string& fileSpectrum()const{return fileSpectrum_;};
		std::string& fileDipoleT(){return fileDipoleT_;};
		const std::string& fileDipoleT()const{return fileDipoleT_;};
		std::string& fileChgT(){return fileChgT_;};
		const std::string& fileChgT()const{return fileChgT_;};
		std::string& fileChgA(){return fileChgA_;};
		const std::string& fileChgA()const{return fileChgA_;};
	//profile
		Profile& profile(){return profile_;};
		const Profile& profile()const{return profile_;};
	//temp
		double& T(){return T_;}
		const double& T()const{return T_;}
	//subset
		std::vector<unsigned int>& subset(){return subset_;}
		const std::vector<unsigned int>& subset()const{return subset_;}
	
	//member functions
		void defaults();
		void clear(){defaults();};
		void init(Simulation& sim);
	//spectra
		void calcChg(Simulation& sim, const QEQ& qeq, const Ewald3D::Coulomb& ewald);
		void calcSpectrum(Simulation& sim);
	//loading/printing
		void read(const char* file);
		void printSpectrum(const char* file)const;
};

#endif
