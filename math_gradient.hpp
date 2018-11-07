#ifndef MATH_GRADIENT_HPP
#define MATH_GRADIENT_HPP

#include <cstdlib>
#include <vector>
#include <Eigen/Dense>

namespace gradient{
	
	//*********************************************************
	//Vectors - Scalars
	//*********************************************************
	
	/*
		Symmetric first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T dc1o2(const std::vector<T>& v, T step, unsigned int t){
		return 0.5*(v[t+1]-v[t-1])/step;
	}
	template <class T> T dc1o4(const std::vector<T>& v, T step, unsigned int t){
		return (1.0/12.0*v[t-2]-2.0/3.0*v[t-1]+2.0/3.0*v[t+1]-1.0/12.0*v[t+2])/step;
	}
	template <class T> T dc1o6(const std::vector<T>& v, T step, unsigned int t){
		return (-1.0/60.0*v[t-3]+3.0/20.0*v[t-2]-3.0/4.0*v[t-1]+3.0/4.0*v[t+1]-3.0/20.0*v[t+2]+1.0/60.0*v[t+3])/step;
	}
	template <class T> T dc1o8(const std::vector<T>& v, T step, unsigned int t){
		return (1.0/280.0*v[t-4]-4.0/105.0*v[t-3]+1.0/5.0*v[t-2]-4.0/5.0*v[t-1]+4.0/5.0*v[t+1]-1.0/5.0*v[t+2]+4.0/105.0*v[t+3]-1.0/280.0*v[t+4])/step;
	}
	
	/*
		Forward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T df1o1(const std::vector<T>& v, T step, unsigned int t){
		return (-v[t]+v[t+1])/step;
	}
	template <class T> T df1o2(const std::vector<T>& v, T step, unsigned int t){
		return (-3.0/2.0*v[t]+2*v[t+1]-1.0/2.0*v[t+2])/step;
	}
	template <class T> T df1o3(const std::vector<T>& v, T step, unsigned int t){
		return (-11.0/6.0*v[t]+3*v[t+1]-3.0/2.0*v[t+2]+1.0/3.0*v[t+3])/step;
	}
	template <class T> T df1o4(const std::vector<T>& v, T step, unsigned int t){
		return (-25.0/12.0*v[t]+4*v[t+1]-3*v[t+2]+4.0/3.0*v[t+3]-1.0/4.0*v[t+4])/step;
	}
	template <class T> T df1o5(const std::vector<T>& v, T step, unsigned int t){
		return (-137.0/60.0*v[t]+5*v[t+1]-5*v[t+2]+10.0/3.0*v[t+3]-5.0/4.0*v[t+4]+1.0/5.0*v[t+5])/step;
	}
	template <class T> T df1o6(const std::vector<T>& v, T step, unsigned int t){
		return (-49.0/20.0*v[t]+6*v[t+1]-15.0/2.0*v[t+2]+20.0/3.0*v[t+3]-15.0/4.0*v[t+4]+6.0/5.0*v[t+5]-1.0/6.0*v[t+6])/step;
	}
	
	/*
		Backward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T db1o1(const std::vector<T>& v, T step, unsigned int t){
		return (v[t-1]+v[t])/step;
	}
	template <class T> T db1o2(const std::vector<T>& v, T step, unsigned int t){
		return (1.0/2.0*v[t-2]-2*v[t-1]+3.0/2.0*v[t])/step;
	}
	template <class T> T db1o3(const std::vector<T>& v, T step, unsigned int t){
		return (-1.0/3.0*v[t-3]+3.0/2.0*v[t-2]-3*v[t-1]+11.0/6.0*v[t])/step;
	}
	template <class T> T db1o4(const std::vector<T>& v, T step, unsigned int t){
		return (1.0/4.0*v[t-4]-4.0/3.0*v[t-3]+3*v[t-2]-4*v[t-1]+25.0/12.0*v[t])/step;
	}
	template <class T> T db1o5(const std::vector<T>& v, T step, unsigned int t){
		return (-1.0/5.0*v[t-5]+5.0/4.0*v[t-4]-10.0/3.0*v[t-3]+5*v[t-2]-5*v[t-1]+137.0/60.0*v[t])/step;
	}
	template <class T> T db1o6(const std::vector<T>& v, T step, unsigned int t){
		return (1.0/6.0*v[t-6]-6.0/5.0*v[t-5]+15.0/4.0*v[t-4]-20.0/3.0*v[t-3]+15.0/2.0*v[t-2]-6*v[t-1]+49.0/20.0*v[t])/step;
	}
	
	//*********************************************************
	//Vectors - Eigen
	//*********************************************************
	
	/*
		Symmetric first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& dc1o2(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=0.5*(r[t+1]-r[t-1])/step;
	}
	inline Eigen::Vector3d& dc1o4(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/12.0*r[t-2]-2.0/3.0*r[t-1]+2.0/3.0*r[t+1]-1.0/12.0*r[t+2])/step;
	}
	inline Eigen::Vector3d& dc1o6(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/60.0*r[t-3]+3.0/20.0*r[t-2]-3.0/4.0*r[t-1]+3.0/4.0*r[t+1]-3.0/20.0*r[t+2]+1.0/60.0*r[t+3])/step;
	}
	inline Eigen::Vector3d& dc1o8(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/280.0*r[t-4]-4.0/105.0*r[t-3]+1.0/5.0*r[t-2]-4.0/5.0*r[t-1]+4.0/5.0*r[t+1]-1.0/5.0*r[t+2]+4.0/105.0*r[t+3]-1.0/280.0*r[t+4])/step;
	}
	
	/*
		Forward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& df1o1(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-r[t]+r[t+1])/step;
	}
	inline Eigen::Vector3d& df1o2(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-3.0/2.0*r[t]+2*r[t+1]-1.0/2.0*r[t+2])/step;
	}
	inline Eigen::Vector3d& df1o3(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-11.0/6.0*r[t]+3*r[t+1]-3.0/2.0*r[t+2]+1.0/3.0*r[t+3])/step;
	}
	inline Eigen::Vector3d& df1o4(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-25.0/12.0*r[t]+4*r[t+1]-3*r[t+2]+4.0/3.0*r[t+3]-1.0/4.0*r[t+4])/step;
	}
	inline Eigen::Vector3d& df1o5(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-137.0/60.0*r[t]+5*r[t+1]-5*r[t+2]+10.0/3.0*r[t+3]-5.0/4.0*r[t+4]+1.0/5.0*r[t+5])/step;
	}
	inline Eigen::Vector3d& df1o6(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-49.0/20.0*r[t]+6*r[t+1]-15.0/2.0*r[t+2]+20.0/3.0*r[t+3]-15.0/4.0*r[t+4]+6.0/5.0*r[t+5]-1.0/6.0*r[t+6])/step;
	}
	
	/*
		Backward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& db1o1(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(r[t-1]+r[t])/step;
	}
	inline Eigen::Vector3d& db1o2(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/2.0*r[t-2]-2*r[t-1]+3.0/2.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o3(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/3.0*r[t-3]+3.0/2.0*r[t-2]-3*r[t-1]+11.0/6.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o4(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/4.0*r[t-4]-4.0/3.0*r[t-3]+3*r[t-2]-4*r[t-1]+25.0/12.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o5(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/5.0*r[t-5]+5.0/4.0*r[t-4]-10.0/3.0*r[t-3]+5*r[t-2]-5*r[t-1]+137.0/60.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o6(const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/6.0*r[t-6]-6.0/5.0*r[t-5]+15.0/4.0*r[t-4]-20.0/3.0*r[t-3]+15.0/2.0*r[t-2]-6*r[t-1]+49.0/20.0*r[t])/step;
	}
	
	//*********************************************************
	//Arrays - Scalar
	//*********************************************************
	
	/*
		Symmetric first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T dc1o2(const T* v, T step, unsigned int t){
		return 0.5*(v[t+1]-v[t-1])/step;
	}
	template <class T> T dc1o4(const T* v, T step, unsigned int t){
		return (1.0/12.0*v[t-2]-2.0/3.0*v[t-1]+2.0/3.0*v[t+1]-1.0/12.0*v[t+2])/step;
	}
	template <class T> T dc1o6(const T* v, T step, unsigned int t){
		return (-1.0/60.0*v[t-3]+3.0/20.0*v[t-2]-3.0/4.0*v[t-1]+3.0/4.0*v[t+1]-3.0/20.0*v[t+2]+1.0/60.0*v[t+3])/step;
	}
	template <class T> T dc1o8(const T* v, T step, unsigned int t){
		return (1.0/280.0*v[t-4]-4.0/105.0*v[t-3]+1.0/5.0*v[t-2]-4.0/5.0*v[t-1]+4.0/5.0*v[t+1]-1.0/5.0*v[t+2]+4.0/105.0*v[t+3]-1.0/280.0*v[t+4])/step;
	}
	
	/*
		Forward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T df1o1(const T* v, T step, unsigned int t){
		return (-v[t]+v[t+1])/step;
	}
	template <class T> T df1o2(const T* v, T step, unsigned int t){
		return (-3.0/2.0*v[t]+2*v[t+1]-1.0/2.0*v[t+2])/step;
	}
	template <class T> T df1o3(const T* v, T step, unsigned int t){
		return (-11.0/6.0*v[t]+3*v[t+1]-3.0/2.0*v[t+2]+1.0/3.0*v[t+3])/step;
	}
	template <class T> T df1o4(const T* v, T step, unsigned int t){
		return (-25.0/12.0*v[t]+4*v[t+1]-3*v[t+2]+4.0/3.0*v[t+3]-1.0/4.0*v[t+4])/step;
	}
	template <class T> T df1o5(const T* v, T step, unsigned int t){
		return (-137.0/60.0*v[t]+5*v[t+1]-5*v[t+2]+10.0/3.0*v[t+3]-5.0/4.0*v[t+4]+1.0/5.0*v[t+5])/step;
	}
	template <class T> T df1o6(const T* v, T step, unsigned int t){
		return (-49.0/20.0*v[t]+6*v[t+1]-15.0/2.0*v[t+2]+20.0/3.0*v[t+3]-15.0/4.0*v[t+4]+6.0/5.0*v[t+5]-1.0/6.0*v[t+6])/step;
	}
	
	/*
		Backward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	template <class T> T db1o1(const T* v, T step, unsigned int t){
		return (v[t-1]+v[t])/step;
	}
	template <class T> T db1o2(const T* v, T step, unsigned int t){
		return (1.0/2.0*v[t-2]-2*v[t-1]+3.0/2.0*v[t])/step;
	}
	template <class T> T db1o3(const T* v, T step, unsigned int t){
		return (-1.0/3.0*v[t-3]+3.0/2.0*v[t-2]-3*v[t-1]+11.0/6.0*v[t])/step;
	}
	template <class T> T db1o4(const T* v, T step, unsigned int t){
		return (1.0/4.0*v[t-4]-4.0/3.0*v[t-3]+3*v[t-2]-4*v[t-1]+25.0/12.0*v[t])/step;
	}
	template <class T> T db1o5(const T* v, T step, unsigned int t){
		return (-1.0/5.0*v[t-5]+5.0/4.0*v[t-4]-10.0/3.0*v[t-3]+5*v[t-2]-5*v[t-1]+137.0/60.0*v[t])/step;
	}
	template <class T> T db1o6(const T* v, T step, unsigned int t){
		return (1.0/6.0*v[t-6]-6.0/5.0*v[t-5]+15.0/4.0*v[t-4]-20.0/3.0*v[t-3]+15.0/2.0*v[t-2]-6*v[t-1]+49.0/20.0*v[t])/step;
	}
	
	//*********************************************************
	//Arrays - Eigen
	//*********************************************************
	
	/*
		Symmetric first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& dc1o2(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=0.5*(r[t+1]-r[t-1])/step;
	}
	inline Eigen::Vector3d& dc1o4(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/12.0*r[t-2]-2.0/3.0*r[t-1]+2.0/3.0*r[t+1]-1.0/12.0*r[t+2])/step;
	}
	inline Eigen::Vector3d& dc1o6(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/60.0*r[t-3]+3.0/20.0*r[t-2]-3.0/4.0*r[t-1]+3.0/4.0*r[t+1]-3.0/20.0*r[t+2]+1.0/60.0*r[t+3])/step;
	}
	inline Eigen::Vector3d& dc1o8(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/280.0*r[t-4]-4.0/105.0*r[t-3]+1.0/5.0*r[t-2]-4.0/5.0*r[t-1]+4.0/5.0*r[t+1]-1.0/5.0*r[t+2]+4.0/105.0*r[t+3]-1.0/280.0*r[t+4])/step;
	}
	
	/*
		Forward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& df1o1(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-r[t]+r[t+1])/step;
	}
	inline Eigen::Vector3d& df1o2(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-3.0/2.0*r[t]+2*r[t+1]-1.0/2.0*r[t+2])/step;
	}
	inline Eigen::Vector3d& df1o3(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-11.0/6.0*r[t]+3*r[t+1]-3.0/2.0*r[t+2]+1.0/3.0*r[t+3])/step;
	}
	inline Eigen::Vector3d& df1o4(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-25.0/12.0*r[t]+4*r[t+1]-3*r[t+2]+4.0/3.0*r[t+3]-1.0/4.0*r[t+4])/step;
	}
	inline Eigen::Vector3d& df1o5(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-137.0/60.0*r[t]+5*r[t+1]-5*r[t+2]+10.0/3.0*r[t+3]-5.0/4.0*r[t+4]+1.0/5.0*r[t+5])/step;
	}
	inline Eigen::Vector3d& df1o6(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-49.0/20.0*r[t]+6*r[t+1]-15.0/2.0*r[t+2]+20.0/3.0*r[t+3]-15.0/4.0*r[t+4]+6.0/5.0*r[t+5]-1.0/6.0*r[t+6])/step;
	}
	
	/*
		Backward first-order derivatives
		o(n) corresponds to error in h (i.e. error~O(h^n))
	*/
	
	inline Eigen::Vector3d& db1o1(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(r[t-1]+r[t])/step;
	}
	inline Eigen::Vector3d& db1o2(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/2.0*r[t-2]-2*r[t-1]+3.0/2.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o3(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/3.0*r[t-3]+3.0/2.0*r[t-2]-3*r[t-1]+11.0/6.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o4(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/4.0*r[t-4]-4.0/3.0*r[t-3]+3*r[t-2]-4*r[t-1]+25.0/12.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o5(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(-1.0/5.0*r[t-5]+5.0/4.0*r[t-4]-10.0/3.0*r[t-3]+5*r[t-2]-5*r[t-1]+137.0/60.0*r[t])/step;
	}
	inline Eigen::Vector3d& db1o6(const Eigen::Vector3d*& r, Eigen::Vector3d& v, double step, unsigned int t){
		return v.noalias()=(1.0/6.0*r[t-6]-6.0/5.0*r[t-5]+15.0/4.0*r[t-4]-20.0/3.0*r[t-3]+15.0/2.0*r[t-2]-6*r[t-1]+49.0/20.0*r[t])/step;
	}
}

#endif