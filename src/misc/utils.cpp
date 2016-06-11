/*
 * utils.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#include "utils.h"
#include <iostream>
#include <ctime>

void timed_log(std::string msg){
	  time_t rawtime;
	  time (&rawtime);
	  std::string time_str(ctime(&rawtime));
	log_handle_mutex.lock();
	  (*logFileHandle) << "[" << time_str.substr(0,time_str.length()-1) << "]: " << msg << std::endl;
    log_handle_mutex.unlock();
}

void verbose_timed_log(std::string msg) {
	if (prog_params.isVerbose()){
		timed_log(msg);
	}
}

void timed_log(const gsl_vector* vec, std::string identifier) {
	  time_t rawtime;
	  time (&rawtime);
	  std::string time_str(ctime(&rawtime));
    log_handle_mutex.lock();
	  (*logFileHandle) << "[" << time_str.substr(0,time_str.length()-1) << "]: " << identifier <<":\t";

	  for (size_t i=0;i<vec->size;i++){
		  (*logFileHandle) << gsl_vector_get(vec,i) << "\t";
	  }
	  (*logFileHandle) << std::endl;
    log_handle_mutex.unlock();
}

void verbose_timed_log(const gsl_vector* vec, std::string identifier) {
	if(prog_params.isVerbose()){
		timed_log(vec,identifier);
	}
}
