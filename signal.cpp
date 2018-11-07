#include "signal.hpp"

namespace filter{
	
	std::vector<double>& Gauss::gen(double s, std::vector<double>& f){
		Gauss gauss(s,0.0);
		//int edge=std::ceil(gauss.s*std::sqrt(-2.0*std::log(0.05/gauss.N)));
		int edge=3*s;
		int l=2*edge+1;
		f.resize(l,0);
		for(int i=0; i<l; ++i){
			f[i]=gauss(i-edge);
		}
	}
	
	ArrayT& filter(ArrayT& out, const ArrayT& in, const ArrayT& f, bool pad){
		int edge=f.size()/2;
		double norm=0;
		for(unsigned int i=0; i<f.size(); ++i) norm+=f[i];
		if(pad){
			out.resize(in.size(),0);
			for(int i=0; i<edge; ++i){
				for(int j=edge-i; j<f.size(); ++j){
					out[i]+=in[i+(j-edge)]*f[j];
				}
			}
			for(int i=edge; i<in.size()-edge-f.size()%2; ++i){
				for(int j=0; j<f.size(); ++j){
					out[i]+=in[i+(j-edge)]*f[j];
				}
			}
			for(int i=in.size()-edge-f.size()%2; i<in.size(); ++i){
				for(int j=0; j<f.size()-(i-(in.size()-edge-f.size()%2)); ++j){
					out[i]+=in[i+(j-edge)]*f[j];
				}
			}
		} else {
			
		}
		return out;
	}
}

namespace window{
	
	WINDOW_FUNC::type WINDOW_FUNC::load(const char* str){
		if(std::strcmp(str,"IDENTITY")==0) return WINDOW_FUNC::IDENTITY;
		else if(std::strcmp(str,"GAUSSIAN")==0) return WINDOW_FUNC::GAUSSIAN;
		else if(std::strcmp(str,"KAISER-BESSEL")==0) return WINDOW_FUNC::KAISERBESSEL;
		else if(std::strcmp(str,"BLACKMAN-HARRIS")==0) return WINDOW_FUNC::BLACKMANHARRIS;
		else return WINDOW_FUNC::UNKNOWN;
	}
	
	std::ostream& operator<<(std::ostream& out, const WINDOW_FUNC::type t){
		if(t==WINDOW_FUNC::IDENTITY) out<<"IDENTITY";
		else if(t==WINDOW_FUNC::GAUSSIAN) out<<"GAUSSIAN";
		else if(t==WINDOW_FUNC::KAISERBESSEL) out<<"KAISER-BESSEL";
		else if(t==WINDOW_FUNC::BLACKMANHARRIS) out<<"BLACKMAN-HARRIS";
		else out<<"UNKNOWN";
		return out;
	}
}