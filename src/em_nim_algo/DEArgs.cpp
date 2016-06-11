/*
 * NIMArgs.cpp
 *
 *  Created on: Jul 19, 2015
 *      Author: cheng
 */

#include "DEArgs.h"

DEArgs::DEArgs(std::vector<double> &kappa,
			   std::vector<double> &tau,
			   std::vector<double> &alpha,
			   std::vector<double> &beta,
			   std::vector<double> &scale_factors, std::vector<double> &y,
			   gsl_matrix *x, size_t sample_size, double epsabs, bool laplace) {
	// TODO Auto-generated constructor stub
	max_iter_cquad=2000;
	this->kappa = kappa;
	this->tau = tau;
	this->alpha = alpha;
	this->beta = beta;
	this->scale_factors = scale_factors;
	this->y = y;
	this->x = x;
	this->sample_size = sample_size;
	this->epsabs = epsabs;
	this->laplace = laplace;
}

DEArgs::~DEArgs() {
	// TODO Auto-generated destructor stub
}

