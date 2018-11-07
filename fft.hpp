#ifndef FFT_HPP
#define FFT_HPP

//c libraries
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
//c++ libraries
#include <stdexcept>
#include <iostream>
//fftw3
#include <fftw3.h>

namespace fourier{
	
	//***********************************************************************************************************************************
	//Constants
	//***********************************************************************************************************************************

	static const double mevPerThz=0.24180;
	static const double cmiPerThz=0.02998;

	//***********************************************************************************************************************************
	//FreqUnit struct
	//***********************************************************************************************************************************
	
	struct FreqUnit{
		enum type{
			THZ,//terahertz
			MEV,//millielectronvolts
			CMI,//inverse centimeters
			UNKNOWN
		};
		static FreqUnit::type load(const char* str);
	};
	
	std::ostream& operator<<(std::ostream& out, const FreqUnit::type& t);
	
	//***********************************************************************************************************************************
	//TransformT struct
	//***********************************************************************************************************************************
	
	struct TransformT{
		enum type{
			EXP,//standard discrete fourier transform
			COS,//discrete cosine transform
			SIN,//discete sine transform
			UNKNOWN,
		};
	};
	
	std::ostream& operator<<(std::ostream& out, const TransformT::type& t);
	
	//******************************************************************
	//Class FFT
	//******************************************************************

	class FFT{
	protected:
		unsigned int N_;//the number of points 
		fftw_plan planf;//the fourier transform plan, forward direction
		fftw_plan planr;//the fourier transform plan, reverse direction
	public:
		//constructors/destructors
		FFT():N_(0),planf(NULL),planr(NULL){};
		FFT(unsigned int N):N_(0),planf(NULL),planr(NULL){resize(N);};
		FFT(const FFT& fft):N_(0),planf(NULL),planr(NULL){resize(fft.N());};
		virtual ~FFT();
		
		//operators
		FFT& operator=(const FFT& fft);
		
		//access
		const unsigned int& N()const{return N_;};
		virtual double& inR(unsigned int i)=0;
		virtual double& inI(unsigned int i)=0;
		virtual double& outR(unsigned int i)=0;
		virtual double& outI(unsigned int i)=0;
		
		//member functions
		//resizing
			virtual void resize(unsigned int N);
		//resetting data
			virtual void clear();
		//performing the transform
			virtual void transformf()=0;//forward transform
			virtual void transformr()=0;//backward transform
	};
	
	//******************************************************************
	//Class FFT_C2C
	//******************************************************************
	
	class FFT_C2C: public FFT{
	private:
		fftw_complex* in_;//data - forward transform
		fftw_complex* out_;//data - backward transform
	public:
		FFT_C2C():FFT(),in_(NULL),out_(NULL){};
		FFT_C2C(unsigned int N):FFT(N),in_(NULL),out_(NULL){resize(N);};
		FFT_C2C(const FFT_C2C& fft):FFT(fft.N()),in_(NULL),out_(NULL){resize(fft.N());};
		~FFT_C2C();
		
		//operators
		FFT_C2C& operator=(const FFT_C2C& fft);
		
		//access
		double& inR(unsigned int i){return in_[i][0];};
		double& inI(unsigned int i){return in_[i][1];};
		double& outR(unsigned int i){return out_[i][0];};
		double& outI(unsigned int i){return out_[i][1];};
		const fftw_complex* in()const{return in_;};
		const fftw_complex* out()const{return out_;};
		fftw_complex& in(unsigned int i){return in_[i];};
		const fftw_complex& in(unsigned int i)const{return in_[i];};
		fftw_complex& out(unsigned int i){return out_[i];};
		const fftw_complex& out(unsigned int i)const{return out_[i];};
		
		//member functions
		//resizing
			void resize(unsigned int N);
		//resetting data
			void clear();
		//performing the transform
			void transformf(){fftw_execute(planf);};//forward transform
			void transformr(){fftw_execute(planr);};//backward transform
	};
	
	//******************************************************************
	//Class FFT_R2C
	//******************************************************************
	
	class FFT_R2C: public FFT{
	private:
		double temp_;
		double* in_;//data - forward transform
		fftw_complex* out_;//data - backward transform
	public:
		FFT_R2C():FFT(),in_(NULL),out_(NULL),temp_(0){};
		FFT_R2C(unsigned int N):FFT(N),in_(NULL),out_(NULL),temp_(0){resize(N);};
		FFT_R2C(const FFT_R2C& fft):FFT(fft.N()),in_(NULL),out_(NULL),temp_(0){resize(fft.N());};
		~FFT_R2C();
		
		//operators
		FFT_R2C& operator=(const FFT_R2C& fft);
		
		//access
		double& inR(unsigned int i){return in_[i];};
		double& inI(unsigned int i){return temp_;};
		double& outR(unsigned int i){return out_[i][0];};
		double& outI(unsigned int i){return out_[i][1];};
		const double* in()const{return in_;};
		const fftw_complex* out()const{return out_;};
		double& in(unsigned int i){return in_[i];};
		const double& in(unsigned int i)const{return in_[i];};
		fftw_complex& out(unsigned int i){return out_[i];};
		const fftw_complex& out(unsigned int i)const{return out_[i];};
		
		//member functions
		//resizing
			void resize(unsigned int N);
		//resetting data
			void clear();
		//performing the transform
			void transformf();
			void transformr(){fftw_execute(planr);};//backward transform
	};
	
	//******************************************************************
	//Class FFT_R2R_COS
	//******************************************************************
	
	class FFT_R2R_COS: public FFT{
	private:
		double temp_;
		double* in_;//data - forward transform
		double* out_;//data - backward transform
	public:
		FFT_R2R_COS():FFT(),in_(NULL),out_(NULL){};
		FFT_R2R_COS(unsigned int N):FFT(N),in_(NULL),out_(NULL){resize(N);};
		FFT_R2R_COS(const FFT_R2R_COS& fft):FFT(fft.N()),in_(NULL),out_(NULL){resize(fft.N());};
		~FFT_R2R_COS();
		
		//operators
		FFT_R2R_COS& operator=(const FFT_R2R_COS& fft);
		
		//access
		double& inR(unsigned int i){return in_[i];};
		double& inI(unsigned int i){return temp_;};
		double& outR(unsigned int i){return out_[i];};
		double& outI(unsigned int i){return temp_;};
		const double* in()const{return in_;};
		const double* out()const{return out_;};
		double& in(unsigned int i){return in_[i];};
		const double& in(unsigned int i)const{return in_[i];};
		double& out(unsigned int i){return out_[i];};
		const double& out(unsigned int i)const{return out_[i];};
		
		//member functions
		//resizing
			void resize(unsigned int N);
		//resetting data
			void clear();
		//performing the transform
			void transformf();
			void transformr();
	};
	
	//******************************************************************
	//Class FFT_R2R_SIN
	//******************************************************************
	
	class FFT_R2R_SIN: public FFT{
	private:
		double temp_;
		double* in_;//data - forward transform
		double* out_;//data - backward transform
	public:
		FFT_R2R_SIN():FFT(),in_(NULL),out_(NULL){};
		FFT_R2R_SIN(unsigned int N):FFT(N),in_(NULL),out_(NULL){resize(N);};
		FFT_R2R_SIN(const FFT_R2R_SIN& fft):FFT(fft.N()),in_(NULL),out_(NULL){resize(fft.N());};
		~FFT_R2R_SIN();
		
		//operators
		FFT_R2R_SIN& operator=(const FFT_R2R_SIN& fft);
		
		//access
		double& inR(unsigned int i){return in_[i];};
		double& inI(unsigned int i){return temp_;};
		double& outR(unsigned int i){return out_[i];};
		double& outI(unsigned int i){return temp_;};
		const double* in()const{return in_;};
		const double* out()const{return out_;};
		double& in(unsigned int i){return in_[i];};
		const double& in(unsigned int i)const{return in_[i];};
		double& out(unsigned int i){return out_[i];};
		const double& out(unsigned int i)const{return out_[i];};
		
		//member functions
		//resizing
			void resize(unsigned int N);
		//resetting data
			void clear();
		//performing the transform
			void transformf();
			void transformr();
	};
}

#endif
