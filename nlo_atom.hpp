#ifndef NLO_ATOM_HPP
#define NLO_ATOM_HPP

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
#include "sim.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "sim_util.hpp"
//loading
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
#include "electrostatics_util.hpp"
#include "ewald3D.hpp"
//polarizability
#include "thole.hpp"
//charge transfer
#include "qeq3.hpp"
//units
#include "units.hpp"

#ifndef DEBUG_NLO_ATOM
#define DEBUG_NLO_ATOM 1
#endif

//***********************************************************************************************************************************
//typedefs
//***********************************************************************************************************************************

typedef Atom<Name,AN,Species,Index,Mass,Position,Charge,Alpha,JZero> AtomT;
typedef Molecule<AtomT,Index,Charge,Mass,Position,Dipole> MolT;

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
//NLO
//***********************************************************************************************************************************

class NLO{
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
		unsigned int nStepsChg_;//the total number of calculated steps 
		unsigned int strideAlpha_;//the number of timesteps between sampled data points
		unsigned int nStepsAlpha_;//the total number of calculated steps 
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
		bool calcChg_;//whether we will calculate atomic charges
		bool calcAlpha_;//whether we will calculate atomic alphas
		bool calcSpectrum_;//whether we will calculate the spectrum
		bool normalize_;//whether to normalize the spectrum
	//i/o flags
		bool printAlphaT_;
		bool printDipoleT_;
		bool printChgT_;
		bool printChg_;
	//ir spectrum
		std::vector<std::vector<std::vector<std::vector<Complex> > > > chi2_;
	//logging
		logging::DebugLogger log;
public:
	//constants
	static const double mevPerThz;
	static const double cmiPerThz;
	
	//costructors/destructors
	NLO():log("NLO"){defaults();};
	~NLO(){};
	
	//operators
	NLO& operator=(const NLO& nlo2d);
	friend std::ostream& operator<<(std::ostream& out, const NLO& nlo2d);
	
	//access
	//calculation flags
		bool& calcChg(){return calcChg_;};
		const bool& calcChg()const{return calcChg_;};
		bool& calcAlpha(){return calcAlpha_;};
		const bool& calcAlpha()const{return calcAlpha_;};
		bool& calcSpectrum(){return calcSpectrum_;};
		const bool& calcSpectrum()const{return calcSpectrum_;};
		bool& normalize(){return normalize_;};
		const bool& normalize()const{return normalize_;};
	//i/o flags
		bool& printAlphaT(){return printAlphaT_;};
		const bool& printAlphaT()const{return printAlphaT_;};
		bool& printDipoleT(){return printDipoleT_;};
		const bool& printDipoleT()const{return printDipoleT_;};
		bool& printChgT(){return printChgT_;};
		const bool& printChgT()const{return printChgT_;};
		bool& printChg(){return printChg_;};
		const bool& printChg()const{return printChg_;};
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
		void calcChg(SimAtomic<AtomT>& sim, const QEQ3& qeq, const Ewald3D::Coulomb& ewald);
		void calcAlpha(SimAtomic<AtomT>& sim, const Thole& thole, const Ewald3D::Dipole& ewald);
		void calcSpectrum(SimAtomic<AtomT>& sim, Bonding& bonding);
	//loading/printing
		void load(const char* file);
		void printSpectrum(const char* file)const;
};

#endif
