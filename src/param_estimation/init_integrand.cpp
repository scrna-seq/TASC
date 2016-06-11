/*
 * init_integrand.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: cheng
 */


#include "init_integrand.h"
#include "IntegrationArgs.h"
#include "InitMargLklhoodArgs.h"
#include "../cubature/cubature.h"
#include "../misc/utils.h"

#include <cmath>

#include <boost/math/distributions/negative_binomial.hpp>


int pdf_single_cell(unsigned ndim, const double *x, void *fdata,
		unsigned fdim, double *fval){
	double kappa = *x;
	double tau = *(x+1);
	double jacobian = (1+kappa*kappa)/(1-kappa*kappa)/(1-kappa*kappa)*(1+tau*tau)/(1-tau*tau)/(1-tau*tau);
	kappa = kappa/(1-kappa*kappa);
	tau = tau/(1-tau*tau);

	if((!std::isfinite(kappa)) || (!std::isfinite(tau))){
		*fval=0;
		return 0;
	}

	IntegrationArgs *int_args = (IntegrationArgs *)fdata;

	double prod = 1;
	double observed_prob;
	double expected_lambda;

	for (size_t i=0;i<int_args->y.size();++i){
		//		std::cout << int_args->y.at(i) << ":" << int_args->true_concentrations.at(i) << ":"<<prod <<" ";
		if (int_args->y.at(i)>0){
			prod = prod*expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
		} else {
			observed_prob = expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
			expected_lambda = std::exp(int_args->alpha+
									   int_args->beta*std::log(int_args->true_concentrations.at(i)));
			prod = prod*(1-observed_prob+observed_prob*std::exp(-expected_lambda));
			//			std::cout <<kappa+tau*log(int_args->true_concentrations.at(i)) <<":"<<observed_prob << ":" << expected_lambda <<":"<<prod<<"\t";
		}
	}


	gsl_vector *kappa_tau = gsl_vector_alloc(2);
	gsl_vector *mean_kappa_tau = gsl_vector_alloc(2);
	gsl_matrix *sigma = gsl_matrix_alloc(2,2);
	gsl_vector_set(kappa_tau,0,kappa);
	gsl_vector_set(kappa_tau,1,tau);
	gsl_vector_set(mean_kappa_tau,0,int_args->mean_kappa);
	gsl_vector_set(mean_kappa_tau,1,int_args->mean_tau);
	gsl_matrix_set(sigma,0,0,int_args->sd_kappa*int_args->sd_kappa);
	gsl_matrix_set(sigma,0,1,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,0,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,1,int_args->sd_tau*int_args->sd_tau);
	prod = prod * dmvnorm(kappa_tau, mean_kappa_tau, sigma);
	//std::cout << "prod2:"<< prod << "\t";
	*fval = prod*jacobian;
	//std::cout << std::endl<< prod << "\t";
	//std::cout << "kappa: " << kappa << "\ttau: "<<tau<<"\tfval: "<<*fval<<std::endl;
	gsl_vector_free(kappa_tau);
	gsl_vector_free(mean_kappa_tau);
	gsl_matrix_free(sigma);
	return 0;
}

int kappa_times_pdf_single_cell(unsigned ndim, const double* x, void* fdata,
		unsigned fdim, double* fval) {
	double kappa = *x;
	double tau = *(x+1);
	double jacobian = (1+kappa*kappa)/(1-kappa*kappa)/(1-kappa*kappa)*(1+tau*tau)/(1-tau*tau)/(1-tau*tau);
	kappa = kappa/(1-kappa*kappa);
	tau = tau/(1-tau*tau);

	if((!std::isfinite(kappa)) || (!std::isfinite(tau))){
		*fval=0;
		return 0;
	}


	IntegrationArgs *int_args = (IntegrationArgs *)fdata;

	double prod = 1;
	double observed_prob;
	double expected_lambda;
	for (size_t i=0;i<int_args->y.size();++i){
		//		std::cout << int_args->y.at(i) << ":" << int_args->true_concentrations.at(i) << ":"<<prod <<" ";
		if (int_args->y.at(i)>0){
			prod = prod*expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
		} else {
			observed_prob = expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
			expected_lambda = std::exp(int_args->alpha+
									   int_args->beta*std::log(int_args->true_concentrations.at(i)));
			prod = prod*(1-observed_prob+observed_prob*std::exp(-expected_lambda));
			//			std::cout <<kappa+tau*log(int_args->true_concentrations.at(i)) <<":"<<observed_prob << ":" << expected_lambda <<":"<<prod<<"\t";
		}
	}


	gsl_vector *kappa_tau = gsl_vector_alloc(2);
	gsl_vector *mean_kappa_tau = gsl_vector_alloc(2);
	gsl_matrix *sigma = gsl_matrix_alloc(2,2);
	gsl_vector_set(kappa_tau,0,kappa);
	gsl_vector_set(kappa_tau,1,tau);
	gsl_vector_set(mean_kappa_tau,0,int_args->mean_kappa);
	gsl_vector_set(mean_kappa_tau,1,int_args->mean_tau);
	gsl_matrix_set(sigma,0,0,int_args->sd_kappa*int_args->sd_kappa);
	gsl_matrix_set(sigma,0,1,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,0,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,1,int_args->sd_tau*int_args->sd_tau);
	prod = prod * dmvnorm(kappa_tau, mean_kappa_tau, sigma);
	//std::cout << std::endl<< prod << "\t";
	*fval = prod*jacobian*kappa;
	//std::cout << std::endl<< prod << "\t";
	//std::cout << "kappa: " << kappa << "\ttau: "<<tau<<"\tfval: "<<*fval<<std::endl;
	gsl_vector_free(kappa_tau);
	gsl_vector_free(mean_kappa_tau);
	gsl_matrix_free(sigma);
	return 0;
}

int tau_times_pdf_single_cell(unsigned ndim, const double* x, void* fdata,
		unsigned fdim, double* fval) {
	double kappa = *x;
	double tau = *(x+1);
	double jacobian = (1+kappa*kappa)/(1-kappa*kappa)/(1-kappa*kappa)*(1+tau*tau)/(1-tau*tau)/(1-tau*tau);
	kappa = kappa/(1-kappa*kappa);
	tau = tau/(1-tau*tau);

	if((!std::isfinite(kappa)) || (!std::isfinite(tau))){
		*fval=0;
		return 0;
	}


	IntegrationArgs *int_args = (IntegrationArgs *)fdata;

	double prod = 1;
	double observed_prob;
	double expected_lambda;
	for (size_t i=0;i<int_args->y.size();++i){
		//		std::cout << int_args->y.at(i) << ":" << int_args->true_concentrations.at(i) << ":"<<prod <<" ";
		if (int_args->y.at(i)>0){
			prod = prod*expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
		} else {
			observed_prob = expit(kappa+tau*std::log(int_args->true_concentrations.at(i)));
			expected_lambda = std::exp(int_args->alpha+
								  int_args->beta*std::log(int_args->true_concentrations.at(i)));
			prod = prod*(1-observed_prob+observed_prob*std::exp(-expected_lambda));
			//			std::cout <<kappa+tau*log(int_args->true_concentrations.at(i)) <<":"<<observed_prob << ":" << expected_lambda <<":"<<prod<<"\t";
		}
	}

	gsl_vector *kappa_tau = gsl_vector_alloc(2);
	gsl_vector *mean_kappa_tau = gsl_vector_alloc(2);
	gsl_matrix *sigma = gsl_matrix_alloc(2,2);
	gsl_vector_set(kappa_tau,0,kappa);
	gsl_vector_set(kappa_tau,1,tau);
	gsl_vector_set(mean_kappa_tau,0,int_args->mean_kappa);
	gsl_vector_set(mean_kappa_tau,1,int_args->mean_tau);
	gsl_matrix_set(sigma,0,0,int_args->sd_kappa*int_args->sd_kappa);
	gsl_matrix_set(sigma,0,1,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,0,int_args->sd_kappa*int_args->sd_tau*int_args->rho_kappa_tau);
	gsl_matrix_set(sigma,1,1,int_args->sd_tau*int_args->sd_tau);
	prod = prod * dmvnorm(kappa_tau, mean_kappa_tau, sigma);
	//std::cout << std::endl<< prod << "\t";
	*fval = prod*jacobian*tau;
	//std::cout << std::endl<< prod << "\t";
	//std::cout << "kappa: " << kappa << "\ttau: "<<tau<<"\tfval: "<<*fval<<std::endl;
	gsl_vector_free(kappa_tau);
	gsl_vector_free(mean_kappa_tau);
	gsl_matrix_free(sigma);
	return 0;
}

double init_marginal_lklhood (const gsl_vector *x, void * params){
	InitMargLklhoodArgs * init_params = (InitMargLklhoodArgs *) params;
	size_t sample_size = init_params->y.size();
	verbose_timed_log(x,"init_marginal_x:");

	double val[sample_size];
	double err[sample_size];
	double xmin[2] = {-1,-1}, xmax[2] = {1,1};



#pragma omp parallel for schedule(dynamic,1)
	for(size_t idx = 0; idx<sample_size;++idx){
		IntegrationArgs single_lklhood_args;
		single_lklhood_args.SetAllArgs(
				init_params->alphas[idx],
				init_params->betas[idx],
				init_params->y[idx],
				init_params->true_concentrations,
				init_params->ercc_lens,
				gsl_vector_get(x,0),
				gsl_vector_get(x,1),
				exp(gsl_vector_get(x,2)),
				exp(gsl_vector_get(x,3)),
				transform_rho_from_real_domain(gsl_vector_get(x,4)));


		//		single_lklhood_args.displayAll();
		double epsrel = prog_params.getMinRelTolCubature();
		hcubature(1, pdf_single_cell, (void *) (&single_lklhood_args),
				2, xmin, xmax, init_params->cubatureMaxIter, 0,
				epsrel, ERROR_INDIVIDUAL, val+idx, err+idx);
		//verbose_timed_log("Single Likelihood:\t"+std::to_string((long double)val[idx]));
	}

	double log_lklhood = 0;
	for(size_t i=0;i<sample_size;++i){
		log_lklhood+=std::log(val[i]);
	}
	log_lklhood = -log_lklhood;
	verbose_timed_log("Marginal Log Likelihood:\t"+std::to_string((long double)(-log_lklhood)));
	return log_lklhood;
}
