/*
 * utils.h
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <mutex>

#include "CLIArgs.h"
#include <gsl/gsl_matrix.h>

extern std::mutex log_handle_mutex;
extern CLIArgs prog_params;
extern std::ostream *logFileHandle;
void timed_log(std::string);
void timed_log(const gsl_vector *vec, std::string identifier);
void verbose_timed_log(std::string);
template<typename TYPE>
void print_stdvec(std::vector<TYPE> x, std::string identifier){
	std::cout << identifier <<":\t";
	for(size_t i=0;i<x.size();i++){
		std::cout<<x.at(i)<<"\t";
	}
	std::cout << std::endl;
};
template<typename TYPE>
void print_stdmat(std::vector<std::vector<TYPE> > x, std::string identifier){
	std::cout << identifier <<":\t"<<std::endl;
	for(size_t i=0;i<x.size();i++){
		for(size_t j=0;j<x.at(i).size();j++){
			std::cout<<x.at(i).at(j)<<"\t";
		}
		std::cout << std::endl;
	}
};

template<typename TYPE>
std::string getStrOfStdVector(std::vector<TYPE> arr){
	std::stringstream vecStrBuilder;
	for(size_t i=0;i<arr.size();i++){
		vecStrBuilder<<arr[i]<<",";
	}
	return vecStrBuilder.str();
};
void verbose_timed_log(const gsl_vector* vec, std::string identifier);
template<typename TYPE>
void timed_log(std::vector<TYPE> arr, std::string identifier){
	  time_t rawtime;
	  time (&rawtime);
	  std::string time_str(ctime(&rawtime));
	  std::stringstream vecStrBuilder;
	  for (size_t i=0;i<arr.size();i++){
		  vecStrBuilder << arr[i] << "\t";
	  }
    log_handle_mutex.lock();
	  (*logFileHandle) << "[" << time_str.substr(0,time_str.length()-1) << "]: " << identifier <<":\t"<<vecStrBuilder.str() << std::endl;
    log_handle_mutex.unlock();
}

template<typename TYPE>
void verbose_timed_log(std::vector<TYPE> arr, std::string identifier){
	if(prog_params.isVerbose()){
		timed_log(arr, identifier);
	}
}

#endif /* UTILS_H_ */
