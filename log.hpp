#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <memory>

namespace logging{
	
	inline void nullDeleter(void* p){}
	
	//**********************************************************************************
	//Sink class
	//**********************************************************************************
	
	class Sink{
	private:
		typedef std::shared_ptr<std::ostream> streamP;
		streamP stream_;
	public:
		//constructors/destructors
		Sink(){stream_=streamP(&std::cout,nullDeleter);}
		Sink(std::ostream* stream){stream_=(stream==&std::cout)?streamP(stream,nullDeleter):streamP(stream);}
		Sink(const Sink& sink):stream_(sink.stream()){};
		virtual ~Sink(){};
		
		//get/set functions
		const streamP& stream() const{return stream_;};
		
		//operators
		const Sink& operator=(const Sink& sink){stream_=sink.stream();}
	};
	
	//**********************************************************************************
	//Log class
	//**********************************************************************************
	
	class Log{
	public:
		static Log& get(){
			static Log log;
			return log;
		}
		//member functions
		void addSink(const Sink& sink){sinks.push_back(sink);};
		int numSinks() const{return sinks.size();};
		const Sink& sink(int i) const{return sinks.at(i);};
		template <class T> const Log& operator<<(const T& object) const;
	private:
		//constructors/destructors
		Log(){};
		~Log(){};
		
		//make sure we disable the copy constructor and assignment operator
		Log(const Log& log);
		void operator=(const Log& log);
		
		//member variables
		std::vector<Sink> sinks;
	};
	
	template <class T>
	const Log& Log::operator<<(const T& object) const{
		for(int i=0; i<sinks.size(); i++){
			*(sinks[i].stream())<<object;
		}
		return *this;
	}
	
	//**********************************************************************************
	//LoggerBase class
	//**********************************************************************************
	
	class LoggerBase{
	protected:
		Log& log_;
	public:
		//constructors/destructors
		LoggerBase():log_(Log::get()){};
		LoggerBase(const LoggerBase& l):log_(Log::get()){};
		virtual ~LoggerBase(){};
		
		//operators
		friend std::ostream& operator<<(std::ostream& stream, const LoggerBase& logger){return logger.print(stream);}
		template <class T> const Log& operator<<(const T& object) const;
		const LoggerBase& operator=(const LoggerBase& t){};
		
		//member functions
		Log& log(){return log_;};
		const Log& log()const{return log_;};
		void flush();
		virtual std::ostream& print(std::ostream& stream) const=0;
	};
	
	template <class T>
	const Log& LoggerBase::operator<<(const T& object) const{
		for(int i=0; i<log_.numSinks(); i++){
			*(log_.sink(i).stream())<<*this<<object;
		}
		return log_;
	}
	
	class DebugLogger: public LoggerBase{
	protected:
		std::string name_;
	public:
		//constructors/destructors
		DebugLogger():name_("DEBUG"){};
		DebugLogger(std::string name):name_(name){};
		DebugLogger(const char* name):name_(std::string(name)){};
		DebugLogger(const DebugLogger& l):name_(l.name()){};
		
		//get/set functions
		const std::string& name() const{return name_;};
		
		//operators
		const DebugLogger& operator=(const DebugLogger& l){name_=l.name();};
		
		//member functions
		std::ostream& print(std::ostream& stream) const{return stream<<name_<<": ";};
	};
	
}

/*namespace record{
	
	//**********************************************************************************
	//Stream class
	//**********************************************************************************
	
	class Stream{
	private:
		typedef std::shared_ptr<std::ostream> StreamP;
		StreamP stream_;
		std::string name_;
	public:
		//constructors/destructors
		Stream();
		Stream(std::ostream* stream);
		Stream(const Stream& stream);
		~Stream(){};
		
		//access
		StreamP& stream(){return stream_;};
		const StreamP& stream()const{return stream_;};
		std::string& name(){return name_;};
		const std::string& name()const{return name_;};
		
		//operators
		Stream& operator=(const Stream& stream);
	};
	
	//**********************************************************************************
	//Log class
	//**********************************************************************************
	
	class Log{
	private:
		//constructors/destructors
		Log(){};
		~Log(){};
		
		//make sure we disable the copy constructor and assignment operator
		Log(const Log& log);
		void operator=(const Log& log);
		
		//member variables
		std::vector<Stream> streams_;
	public:
		//static functions
		static Log& get(){static Log log; return log;}
		
		//operators
		template <class T> const Log& operator<<(const T& object) const;
		
		//access
		unsigned int nStreams()const{return streams_.size();};
		const Stream& stream(unsigned int i)const{return streams_.at(i);};
		
		//member functions
		void add(const Stream& stream){streams_.push_back(stream);};
		void flush();
		template <class T> const Log& append(const T& object);
	};
	
	template <class T>
	const Log& Log::operator<<(const T& object) const{
		for(unsigned int i=0; i<streams_.size(); i++) 
			*streams_[i].stream()<<"["<<streams_[i].name()<<"] "<<object;
		return *this;
	}
	
	template <class T> 
	const Log& Log::append(const T& object){
		for(unsigned int i=0; i<streams_.size(); i++) 
			*streams_[i].stream()<<object;
		return *this;
	}
	
}*/

#endif
