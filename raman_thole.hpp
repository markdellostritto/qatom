#ifndef RAMAN3D_THOLEV_HPP
#define RAMAN3D_THOLEV_HPP

//Thole model, bonds

//standard libraries
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
//eigen libraries
#include <Eigen/Dense>
#include <Eigen/StdVector>
//fftw
#include <fftw3.h>
//simulation
#include "sim.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "sim_util.hpp"
//loading
#include "vasp.hpp"
#include "lammps.hpp"
#include "gaussian.hpp"
#include "string.hpp"
//signal analysis
#include "signal.hpp"
#include "fft.hpp"
//math
#include "math_const.hpp"
#include "math_gradient.hpp"
#include "interpolation.hpp"
//electrostatics
#include "electrostatics_util.hpp"
#include "ewald3D.hpp"
#include "thole.hpp"
//units
#include "units.hpp"

#ifndef DEBUG_RAMAN3D_THOLEV
#define DEBUG_RAMAN3D_THOLEV 2
#endif

//***********************************************************************************************************************************
//typedefs
//***********************************************************************************************************************************

typedef Atom<Name,AN,Species,Index,Position,Alpha> AtomT;

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
//Raman3D
//***********************************************************************************************************************************

class Raman3D{
private:
	//FFT parameters
		fourier::FreqUnit::type freqUnit_;//the unit of the frequency
		double freqCut_;//the cutoff for printing the frequency
		unsigned int freqN_;//the integer cutoff for printing the frequency 
		double minFreq_;//the minimum frequency possible, given the number of timesteps
		double maxFreq_;//the maximum frequency possible, given the number of timesteps
		double freqRes_;//the resolution in frequency space in THz
		unsigned int fourierN_;//the total number of data points sampled
		double freqVis_;//the frequency of the visible laser
	//stride
		unsigned int strideAlpha_;//the number of timesteps between sampled data points
		unsigned int nStepsAlpha_;//the total number of calculated steps for charge
	//FFT windowing
		double sigma_;
		std::function<double (int)> window_;
		window::WINDOW_FUNC::type windowType_;
	//profile
		Profile profileCalc_;
		Profile profileLoad_;
	//file i/o
		std::string fileSpectrum_;//file where the spectrum is printed
	//calculation flags
		bool calcAlpha_;//whether we will calculate the charges
		bool calcSpectrum_;//whether we will calculate the spectrum
		bool normalize_;//whether to normalize the spectrum
	//i/o flags
		bool printAlphaT_;
	//ir spectrum
		std::vector<double> ramanp_;//parallel
		std::vector<double> ramans_;//perpendicular (senkrecht)
		std::vector<double> ramant_;//total (sum of squares)
	//temperature
		double T_;
	//logging
		logging::DebugLogger log;
public:
	//constants
	static const double mevPerThz;
	static const double cmiPerThz;
	
	//costructors/destructors
	Raman3D():log("Raman3D"){defaults();};
	~Raman3D(){};
	
	//operators
	Raman3D& operator=(const Raman3D& ir3d);
	friend std::ostream& operator<<(std::ostream& out, const Raman3D& ir3d);
	
	//access
	//FFT parameters
		double& visFreq(){return freqVis_;};
		const double& visFreq()const{return freqVis_;};
	//temperature
		double& T(){return T_;};
		const double& T()const{return T_;};
	//calculation flags
		bool& calcAlpha(){return calcAlpha_;};
		const bool& calcAlpha()const{return calcAlpha_;};
		bool& calcSpectrum(){return calcSpectrum_;};
		const bool& calcSpectrum()const{return calcSpectrum_;};
		bool& normalize(){return normalize_;};
		const bool& normalize()const{return normalize_;};
	//i/o flags
		bool& printAlphaT(){return printAlphaT_;};
		const bool& printAlphaT()const{return printAlphaT_;};
	//file i/o
		std::string& fileSpectrum(){return fileSpectrum_;};
		const std::string& fileSpectrum()const{return fileSpectrum_;};
	//profile
		Profile& profileCalc(){return profileCalc_;};
		const Profile& profileCalc()const{return profileCalc_;};
	
	//member functions
		void defaults();
		void clear(){defaults();};
		void init(SimAtomic<AtomT>& sim);
	//spectra
		void calcAlpha(SimAtomic<AtomT>& sim, const Thole& thole, const Ewald3D::Dipole& ewald);
		void calcSpectrum(SimAtomic<AtomT>& sim);
		void calcSpectrum2(SimAtomic<AtomT>& sim);
		void calcSpectrum3(SimAtomic<AtomT>& sim);
	//loading/printing
		void load(const char* file);
		void printSpectrum(const char* file)const;
};

#endif
