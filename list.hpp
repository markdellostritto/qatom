#ifndef LIST_HPP
#define LIST_HPP

#include <cstdlib>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace list{
	
	//***************************************************************
	//Sorting: Bubble
	//***************************************************************
	
	template <class T>
	std::vector<T>& bubbleSort(std::vector<T>& arr){
		bool swap;
		if(arr.size()>1){
			do{
				swap=false;
				for(unsigned int i=0; i<arr.size()-1; ++i){
					if(arr[i+1]<arr[i]){
						T temp=arr[i];
						arr[i]=arr[i+1];
						arr[i+1]=temp;
						swap=true;
					}
				}
			}while(swap==true);
		}
		
		return arr;
	}
	
	//***************************************************************
	//Sorting: Insertion
	//***************************************************************
	
	template <class T>
	std::vector<T>& insertionSort(std::vector<T>& arr){
		for(int i=1; i<arr.size(); ++i){
			T x=arr[i];
			int j=i-1;
			while(j>=0 && arr[j]>x){
				arr[j+1]=arr[j];
				j--;
			}
			arr[j+1]=x;
		}
		return arr;
	}
	
	template <class T>
	Eigen::Matrix<T,Eigen::Dynamic,1>& insertionSort(Eigen::Matrix<T,Eigen::Dynamic,1>& arr){
		for(int i=1; i<arr.size(); ++i){
			T x=arr[i];
			int j=i-1;
			while(j>=0 && arr[j]>x){
				arr[j+1]=arr[j];
				j--;
			}
			arr[j+1]=x;
		}
		return arr;
	}
	
	template <class T>
	void insertionSort(std::vector<T>& arr1, std::vector<T>& arr2){
		for(int i=1; i<arr1.size(); ++i){
			T x1=arr1[i];
			T x2=arr2[i];
			int j=i-1;
			while(j>=0 && arr1[j]>x1){
				arr1[j+1]=arr1[j];
				arr2[j+1]=arr2[j];
				j--;
			}
			arr1[j+1]=x1;
			arr2[j+1]=x2;
		}
	}
	
	template <class T, class U>
	void insertionSort(std::vector<T>& arr1, std::vector<U>& arr2){
		for(int i=1; i<arr1.size(); ++i){
			T x=arr1[i];
			int j=i-1;
			while(j>=0 && arr1[j]>x){
				arr1[j+1]=arr1[j];
				arr2[j+1]=arr2[j];
				j--;
			}
			arr1[j+1]=x;
		}
	}
	
	//***************************************************************
	//Sorting: Shell
	//***************************************************************
	
	template <class T>
	std::vector<T>& shellSort(std::vector<T>& arr){
		//find the maximum exponent of the power series of gaps by recursively generating the largest possible gap
		//the maximum gap should be size/3, rounded up [1]
		unsigned int max=arr.size()/3+((arr.size()%3)+1)%2;//this gets us arr/3, rounded up
		unsigned int pow=1,gap=1;
		while(gap<max){
			gap=gap*3+1;
			++pow;
		}
		
		//now that we have the maximum exponenet, perform shell sorting
		for(unsigned int i=pow; i>=1; --i){
			gap=(std::pow(3,i)-1)/2;
			for(int j=gap; j<arr.size(); j+=gap){
				T temp=arr[j];
				unsigned int k=j;
				while(k>0 && arr[k-gap]>temp){
					arr[k]=arr[k-gap];
					k-=gap;
				}
				arr[k]=temp;
			}
		}
		
		return arr;
	}
	
	//***************************************************************
	//Searching
	//***************************************************************
	
	//min/max
	
	template <class T>
	T& min(const std::vector<T>& arr){
		T min=0;
		for(unsigned int i=0; i<arr.size(); ++i) if(min>arr[i]) min=arr[i];
		return min;
	}
	
	template <class T>
	T& max(const std::vector<T>& arr){
		T max=0;
		for(unsigned int i=0; i<arr.size(); ++i) if(max<arr[i]) max=arr[i];
		return max;
	}
	
	//min/max index
	
	template <class T>
	unsigned int minIndex(const std::vector<T>& arr){
		T min=0;
		unsigned int n=0;
		for(unsigned int i=0; i<arr.size(); ++i) if(min>arr[i]){min=arr[i];n=i;}
		return n;
	}
	
	template <class T>
	unsigned int maxIndex(const std::vector<T>& arr){
		T max=0;
		unsigned int n=0;
		for(unsigned int i=0; i<arr.size(); ++i) if(max<arr[i]){max=arr[i];n=i;};
		return n;
	}
	
	//searching - linear 
	
	template <class T>
	unsigned int find_exact_lin(const T& x, std::vector<T>& arr){
		for(unsigned int i=0; i<arr.size(); i++) if(x==arr[i]) return i;
		throw std::runtime_error("No equivalent entry found in arr.");
	}
	
	template <class T>
	unsigned int find_approx_lin(const T& x, std::vector<T>& arr){
		for(unsigned int i=0; i<arr.size(); i++) if(std::fabs(x-arr[i])<1e-6) return i;
		throw std::runtime_error("No equivalent entry found in arr.");
	}
	
	/*References:
		[1] Knuth, Donald E. (1997). "Shell's method". The Art of Computer Programming. 
			Volume 3: Sorting and Searching (2nd ed.). Reading, Massachusetts: Addison-Wesley. 
			pp. 83–95. ISBN 0-201-89685-0.
	*/
	
}

#endif