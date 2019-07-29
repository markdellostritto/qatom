#include "fft.hpp"

namespace fourier{

//***********************************************************************************************************************************
//FreqUnit struct
//***********************************************************************************************************************************

FreqUnit::type FreqUnit::read(const char* str){
	if(std::strcmp(str,"THZ")==0) return FreqUnit::THZ;
	else if(std::strcmp(str,"MEV")==0) return FreqUnit::MEV;
	else if(std::strcmp(str,"CMI")==0) return FreqUnit::CMI;
	else return FreqUnit::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FreqUnit::type& t){
	if(t==FreqUnit::THZ) return out<<"THZ";
	else if(t==FreqUnit::MEV) return out<<"MEV";
	else if(t==FreqUnit::CMI) return out<<"CMI";
	else return out;
}

//***********************************************************************************************************************************
//TransformT struct
//***********************************************************************************************************************************

std::ostream& operator<<(std::ostream& out, const TransformT::type& t){
	if(t==TransformT::EXP) return out<<"EXP";
	else if(t==TransformT::COS) return out<<"COS";
	else if(t==TransformT::SIN) return out<<"SIN";
	else return out;
}

//******************************************************************
//Class FFT
//******************************************************************

FFT::~FFT(){
	if(planf!=NULL) fftw_destroy_plan(planf);
	if(planr!=NULL) fftw_destroy_plan(planr);
}

FFT& FFT::operator=(const FFT& fft){
	resize(fft.N());
	return *this;
}

void FFT::resize(unsigned int N){
	N_=N;
	if(planf!=NULL){fftw_destroy_plan(planf); planf=NULL;}
	if(planr!=NULL){fftw_destroy_plan(planr); planr=NULL;}
}

void FFT::clear(){
	N_=0;
	if(planf!=NULL){fftw_destroy_plan(planf); planf=NULL;}
	if(planr!=NULL){fftw_destroy_plan(planr); planr=NULL;}
}

//******************************************************************
//Class FFT_C2C
//******************************************************************

FFT_C2C::~FFT_C2C(){
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
}

FFT_C2C& FFT_C2C::operator=(const FFT_C2C& fft){
	FFT::operator=(fft);
	resize(N_);
	return *this;
}

void FFT_C2C::resize(unsigned int N){
	FFT::resize(N);
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
	in_=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_);
	out_=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_);
	for(unsigned int i=0; i<N_; ++i){in_[i][0]=0.0; in_[i][1]=0.0;};
	for(unsigned int i=0; i<N_; ++i){out_[i][0]=0.0; out_[i][1]=0.0;};
	planf=fftw_plan_dft_1d(N_,in_,out_,FFTW_FORWARD,FFTW_ESTIMATE);
	planr=fftw_plan_dft_1d(N_,in_,out_,FFTW_BACKWARD,FFTW_ESTIMATE);
}

void FFT_C2C::clear(){
	FFT::clear();
	if(in_!=NULL){fftw_free(in_);in_=NULL;}
	if(out_!=NULL){fftw_free(out_);out_=NULL;}
}

//******************************************************************
//Class FFT_R2C
//******************************************************************

FFT_R2C::~FFT_R2C(){
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
}

FFT_R2C& FFT_R2C::operator=(const FFT_R2C& fft){
	FFT::operator=(fft);
	resize(N_);
	return *this;
}

void FFT_R2C::resize(unsigned int N){
	FFT::resize(N);
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
	in_=(double*)fftw_malloc(sizeof(double)*N_);
	out_=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N_);
	for(unsigned int i=0; i<N_; ++i) in_[i]=0.0;
	for(unsigned int i=0; i<N_; ++i){out_[i][0]=0.0; out_[i][1]=0.0;};
	planf=fftw_plan_dft_r2c_1d(N_,in_,out_,FFTW_ESTIMATE);
	planr=fftw_plan_dft_c2r_1d(N_,out_,in_,FFTW_ESTIMATE);
}

void FFT_R2C::clear(){
	FFT::clear();
	if(in_!=NULL){fftw_free(in_);in_=NULL;}
	if(out_!=NULL){fftw_free(out_);out_=NULL;}
}

//performing the transform

void FFT_R2C::transformf(){
	fftw_execute(planf);
	for(unsigned int t=N_/2+1; t<N_; ++t){
		out_[t][0]=out_[N_-t][0];
		out_[t][1]=-out_[N_-t][1];
	}
}

//******************************************************************
//Class FFT_R2R_COS
//******************************************************************

FFT_R2R_COS::~FFT_R2R_COS(){
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
}

FFT_R2R_COS& FFT_R2R_COS::operator=(const FFT_R2R_COS& fft){
	FFT::operator=(fft);
	resize(N_);
	return *this;
}

void FFT_R2R_COS::resize(unsigned int N){
	FFT::resize(N);
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
	in_=(double*)fftw_malloc(sizeof(double)*N_);
	out_=(double*)fftw_malloc(sizeof(double)*N_);
	for(unsigned int i=0; i<N_; ++i) in_[i]=0.0;
	for(unsigned int i=0; i<N_; ++i) out_[i]=0.0;
	planf=fftw_plan_r2r_1d(N_,in_,out_,FFTW_REDFT10,FFTW_ESTIMATE);
	planr=fftw_plan_r2r_1d(N_,out_,in_,FFTW_REDFT01,FFTW_ESTIMATE);
}

void FFT_R2R_COS::clear(){
	FFT::clear();
	if(in_!=NULL){fftw_free(in_);in_=NULL;}
	if(out_!=NULL){fftw_free(out_);out_=NULL;}
}

//performing the transform

void FFT_R2R_COS::transformf(){
	fftw_execute(planf);
}

void FFT_R2R_COS::transformr(){
	fftw_execute(planr);
};

//******************************************************************
//Class FFT_R2R_SIN
//******************************************************************

FFT_R2R_SIN::~FFT_R2R_SIN(){
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
}

FFT_R2R_SIN& FFT_R2R_SIN::operator=(const FFT_R2R_SIN& fft){
	FFT::operator=(fft);
	resize(N_);
	return *this;
}

void FFT_R2R_SIN::resize(unsigned int N){
	FFT::resize(N);
	if(in_!=NULL) fftw_free(in_);
	if(out_!=NULL) fftw_free(out_);
	in_=(double*)fftw_malloc(sizeof(double)*N_);
	out_=(double*)fftw_malloc(sizeof(double)*N_);
	for(unsigned int i=0; i<N_; ++i) in_[i]=0.0;
	for(unsigned int i=0; i<N_; ++i) out_[i]=0.0;
	planf=fftw_plan_r2r_1d(N_,in_,out_,FFTW_RODFT10,FFTW_ESTIMATE);
	planr=fftw_plan_r2r_1d(N_,out_,in_,FFTW_RODFT01,FFTW_ESTIMATE);
}

void FFT_R2R_SIN::clear(){
	FFT::clear();
	if(in_!=NULL){fftw_free(in_);in_=NULL;}
	if(out_!=NULL){fftw_free(out_);out_=NULL;}
}

//performing the transform

void FFT_R2R_SIN::transformf(){
	fftw_execute(planf);
}

void FFT_R2R_SIN::transformr(){
	fftw_execute(planr);
};

}
