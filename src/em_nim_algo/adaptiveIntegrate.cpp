/*
 * adaptiveIntegrate.cpp
 *
 *  Created on: Nov 10, 2015
 *      Author: jiacheng
 */

#include "adaptiveIntegrate.h"
#include <cmath>
#include <string>
#include "../misc/utils.h"
//#include "../misc/cheng_quad.h"
#include <boost/math/special_functions/factorials.hpp>

double localLogit(double x){
	return 1.0-1.0/(1+std::exp(x));
}

double gX(double t,double sigma,double mu,double alpha, double beta,double kappa,double tau, double y){
	if (y<boost::math::max_factorial<double>::value){
		// if i can directly calculate y!
		double logPart = y*(alpha+beta*(t*sigma+mu))-std::exp(alpha+beta*(t*sigma+mu))-std::log(boost::math::factorial<double>(y));
		return localLogit(kappa+tau*t*sigma+tau*mu)*std::exp(logPart)*sigma;
	} else {
		//use stirling's approximation when y is large.
		double logPart = y*(alpha+beta*(t*sigma+mu))-std::exp(alpha+beta*(t*sigma+mu))-(y*std::log(y)-y);
		return localLogit(kappa+tau*t*sigma+tau*mu)*std::exp(logPart)*sigma;
	}
}

//			pars[1] = this->kappa.at(i);
//			pars[2] = this->tau.at(i);
//			pars[3] = this->alpha.at(i);
//			pars[4] = this->beta.at(i);
//			pars[5] = gsl_vector_get(theta_t,i);
//			pars[6] = sigma_t;
//			pars[7] = this->scaleFactors[i];
//			pars[8] = y.at(i);

int adaptiveIntegrate(gsl_function* F, double* pars, double epsabs,
		double epsrel, gsl_integration_cquad_workspace* workspace,
		double* result, double* error) {

	unsigned long max_eval= 0;
	int status = gsl_integration_cquad(F,*(pars+5)-40*(*(pars+6)),*(pars+5)+40*(*(pars+6)),epsabs,epsrel,workspace,result,error,&max_eval);
//	verbose_timed_log("GSLINTEGRATIONCQUAD number of evaluations:\t" + std::to_string((long long int)max_eval));

	return status;
}


