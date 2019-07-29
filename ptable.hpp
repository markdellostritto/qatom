#ifndef PERIODIC_TABLE_HPP
#define PERIODIC_TABLE_HPP

#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>

namespace PTable{

//*********************************************
//Periodic Table Data
//*********************************************

static const float inf=std::numeric_limits<float>::max();
static const unsigned int NUM_ELEMENTS=118;

static const char* const ELEMENT_NAME[NUM_ELEMENTS]={
	"H","He","Li","Be","B","C","N","O","F","Ne",
	"Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
	"Sc","Ti","V","Cr","Mn","Fe","Ni","Co","Cu","Zn",
	"Ga","Ge","As","Se","Br","Xe","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
	"Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
	"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	"Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
	"Tl","Pb","Bi","Po","At","Rn","Fr","Ra'","Ac","Th",
	"Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
	"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
	"Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
};

static const float ELEMENT_MASS[NUM_ELEMENTS]={
	1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994,18.9984,20.1797,
	22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453,39.948,39.0983,40.078,
	44.9559,47.867,50.9415,51.9961,54.938,55.845,58.6934,58.9332,63.546,65.39,
	69.723,72.64,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,91.224,
	92.9064,95.94,98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,
	121.76,126.9045,127.6,131.293,132.9055,137.327,138.9055,140.116,140.9077,144.24,
	145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,
	174.967,178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,
	204.3833,207.2,208.9804,209,210,222,223,226,227,231.0359,
	232.0381,237,238.0289,243,244,247,247,251,252,257,
	258,259,261,262,262,264,266,268,272,277,
	278,281,282,285,286,289,290,293
};

static const float ATOMIC_RADII[NUM_ELEMENTS]={
	0.53,	0.31,	1.67,	1.12,	0.87,	0.67,	0.56,	0.48,	0.42,	0.38,
	1.9,	1.45,	1.18,	1.11,	0.98,	0.88,	0.79,	0.71,	2.43,	1.94,
	1.84,	1.76,	1.71,	1.66,	1.61,	1.56,	1.52,	1.49,	1.45,	1.42,
	1.36,	1.25,	1.14,	1.03,	0.94,	0.88,	2.65,	2.19,	2.12,	2.06,
	1.98,	1.9,	1.83,	1.78,	1.73,	1.69,	1.65,	1.61,	1.56,	1.45,
	1.33,	1.23,	1.15,	1.08,	2.98,	2.53,	1.95,	1.85,	2.47,	2.06,
	2.05,	2.38,	2.31,	2.33,	2.25,	2.28,	2.26,	2.26,	2.22,	2.22,
	2.17,	2.08,	2,	1.93,	1.88,	1.85,	1.8,	1.77,	1.74,	1.71,
	1.56,	1.54,	1.43,	1.35,	1.27,	1.2,	0,	0,	1.95,	1.8,
	1.8,	1.75,	1.75,	1.75,	1.75,	0
};
//Angstrom
//E Clementi, D L Raimondi, W P Reinhardt (1963) J Chem Phys. 38:2686.

static const float COVALENT_RADII[NUM_ELEMENTS]={
	0.37,	0.32,	1.34,	0.9,	0.82,	0.77,	0.75,	0.73,	0.71,	0.69,
	1.54,	1.3,	1.18,	1.11,	1.06,	1.02,	0.99,	0.97,	1.96,	1.74,
	1.44,	1.36,	1.25,	1.27,	1.39,	1.25,	1.26,	1.21,	1.38,	1.31,
	1.26,	1.22,	1.19,	1.16,	1.14,	1.1,	2.11,	1.92,	1.62,	1.48,
	1.37,	1.45,	1.56,	1.26,	1.35,	1.31,	1.53,	1.48,	1.44,	1.41,
	1.38,	1.35,	1.33,	1.3,	2.25,	1.98,	1.69,	0,	0,	0,
	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
	1.6,	1.5,	1.38,	1.46,	1.59,	1.28,	1.37,	1.28,	1.44,	1.49,
	1.48,	1.47,	1.46,	0,	0,	1.45,	0,	0,	0,	0,
	0,	0,	0,	0,	0,	0
};
//Angstrom
//WebElements

static const float IONIZATION_ENERGY[NUM_ELEMENTS]={
	13.598434005136,24.587387936,5.391714761,9.322699,8.298019,11.2603,14.53413,13.618054,17.42282,21.56454,
	5.1390767,7.646235,5.985768,8.151683,10.486686,10.36001,12.96763,15.7596112,4.34066354,6.1131552,
	6.56149,6.82812,6.746187,6.76651,7.434038,7.9024678,7.88101,7.639877,7.72638,9.3941968,
	5.9993018,7.899435,9.7886,9.752392,11.81381,13.9996049,4.177128,5.6948672,6.21726,6.6339,
	6.75885,7.09243,7.11938,7.3605,7.4589,8.33686,7.576234,8.99382,5.7863554,7.343917,
	8.608389,9.00966,10.45126,12.1298431,3.893905557,5.211664,5.5769,5.5386,5.47,5.525,
	5.577,5.64371,5.670385,6.1498,5.8638,5.93905,6.0215,6.1077,6.18431,6.254159,
	5.425871,6.825069,7.549571,7.86403,7.83352,8.43823,8.96702,8.95883,9.225553,10.437504,
	6.1082871,7.4166796,7.285516,8.414,9.31751,4.0727409,5.278424,5.380226,6.3067,5.89,
	6.19405,6.2655,6.0258,5.9738,5.9914,6.1978,6.2817,6.3676,6.5,6.58,
	6.65,inf,inf,inf,inf,inf,inf,inf,inf,inf,
	inf,inf,inf,inf,inf,inf,inf,inf
};
//ionization energies are in eV
//Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2014). NIST Atomic Spectra Database (ver. 5.2), [Online]. 
//Available: https://physics.nist.gov/asd [2017, August 9]. National Institute of Standards and Technology, Gaithersburg, MD. 

static const float ELECTRON_AFFINITY[NUM_ELEMENTS]={
	0.754195,-0.52,0.618049,-0.52,0.279723,1.2621226,-0.000725,1.4611136,3.4011898,-1.2,
	0.547926,-0.415,0.43283,1.3895212,0.746607,2.0771042,3.612724,-1,0.501459,0.02455,
	0.188,0.084,0.52766,0.67584,-0.52,0.153236,0.66226,1.15716,1.23578,-0.62,
	0.43,1.2326764,0.8048,2.0206047,3.363588,-0.62,0.485916,0.05206,0.307,0.427,
	0.9174,0.7473,0.55,1.04638,1.14289,0.56214,1.30447,-0.725,0.3,1.11207,
	1.047401,1.970875,3.0590465,-0.83,0.47163,0.14462,0.47,0.65,0.962,1.916,
	0.129,0.162,0.864,0.137,1.165,0.352,0.338,0.312,1.029,-0.02,
	0.346,0.017,0.323,0.81626,0.15,1.1,1.56436,2.1251,2.30861,-0.52,
	0.377,0.356743,0.942362,1.9,2.3,-0.725,0.486,0.1,0.35,1.17,
	0.55,0.53,0.48,-0.5,0.1,0.28,-1.72,-1.01,-0.3,0.35,
	0.98,-2.33,-0.31,-inf,-inf,-inf,-inf,-inf,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf
};
//electron affinities are in eV
//wikipedia

static const float ELECTRONEGATIVITY_PAULING[NUM_ELEMENTS]{
	2.2,-inf,0.98,1.57,2.04,2.55,3.04,3.44,3.98,-inf,
	0.93,1.31,1.61,1.9,2.19,2.58,3.16,-inf,0.82,1,
	1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.9,1.65,
	1.81,2.01,2.18,2.55,2.96,3,0.82,0.95,1.22,1.33,
	1.6,2.16,1.9,2.2,2.28,2.2,1.93,1.69,1.78,1.96,
	2.05,2.1,2.66,2.6,0.79,0.89,1.1,1.12,1.13,1.14,
	-inf,1.17,-inf,1.2,-inf,1.22,1.23,1.24,1.25,-inf,
	1.27,1.3,1.5,2.36,1.9,2.2,2.2,2.28,2.54,2,
	1.62,2.33,2.02,2,2.2,-inf,-inf,0.9,1.1,1.3,
	1.5,1.38,1.36,1.28,1.3,1.3,1.3,1.3,1.3,1.3,
	1.3,1.3,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf
};
//electronegativities are in eV
//wikipedia (web elements)

static const float ELECTRONEGATIVITY_ALLEN[NUM_ELEMENTS]{
	2.3,4.16,0.912,1.576,2.051,2.544,3.066,3.61,4.193,4.787,
	0.869,1.293,1.613,1.916,2.253,2.589,2.869,3.242,0.734,1.034,
	1.19,1.38,1.53,1.65,1.75,1.8,1.84,1.88,1.85,1.59,
	1.756,1.994,2.211,2.424,2.685,2.966,0.706,0.963,1.12,1.32,
	1.41,1.47,1.51,1.54,1.56,1.58,1.87,1.52,1.656,1.824,
	1.984,2.158,2.359,2.582,0.659,0.881,-inf,-inf,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,
	1.09,1.16,1.34,1.47,1.6,1.65,1.68,1.72,1.92,1.76,
	1.789,1.854,2.01,2.19,2.39,2.6,0.67,0.89,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,
	-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf
};

//*********************************************
//Enums
//*********************************************

struct ElectronegativityType{
	enum type{
		MULLIKEN,
		PAULING,
		ALLEN,
		HINZE,
		UNKNOWN
	};
	static ElectronegativityType::type read(const char* str);
};

std::ostream& operator<<(std::ostream& out, const ElectronegativityType::type& t);

//*********************************************
//Function
//*********************************************

//************** NAME ***************
const char* name(unsigned int atomicNumber);
//********** ATOMIC_NUMBER **********
unsigned int an(const char* name);
unsigned int an(double mass);
//************** MASS ***************
double mass(unsigned int atomicNumber);
//************* RADIUS **************
double atomicRadius(unsigned int atomicNumber);
double covalentRadius(unsigned int atomicNumber);
//********** IONIZATION_ENERGY **********
//ionization energy
double ionizationEnergy(unsigned int an);
//ionization energy - Hinze
double ionizationEnergyHinze(unsigned int an, unsigned int vs);
//********** ELECTRON_AFFINITY **********
//electron affinity
double electronAffinity(unsigned int an);
//electron affinity - Hinze
double electronAffinityHinze(unsigned int an, unsigned int vs);
//********** ELECTRONEGATIVITY **********
//electronegativity
double electronegativityMulliken(unsigned int an);
//electronegativity - Pauling
double electronegativityPauling(unsigned int an);
//electronegativity - Allen
double electronegativityAllen(unsigned int an);
//electronegativity - Hinze
double electronegativityHinze(unsigned int an, unsigned int vs);
//********** IDEMPOTENTIAL **********
//idempotential (coulomb self-energy)
double idempotential(unsigned int an);
//idempotential Hinze (coulomb self-energy)
double idempotentialHinze(unsigned int an, unsigned int vs);

}

#endif
