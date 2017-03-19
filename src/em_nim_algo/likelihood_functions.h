/*
 * likelihood_functions.h
 *
 *  Created on: Jul 14, 2015
 *      Author: cheng
 */

#ifndef LIKELIHOOD_FUNCTIONS_H_
#define LIKELIHOOD_FUNCTIONS_H_
#include <gsl/gsl_matrix.h>
#include <vector>

double single_lklhood (double x, void * params);
double ew_integrand (double x, void * params);
double ew_sq_integrand (double x, void * params);
double em_marginal_lklhood (std::vector<double> kappa,
		std::vector<double> tau,
		std::vector<double> alpha,
		std::vector<double> beta,
		gsl_vector *theta,
		double sigma,
		std::vector<double> scale_factors,
		std::vector<double> y,
		size_t sample_size,
		double epsabs);
double de_marginal_lklhood(const gsl_vector *x, void *params);
double quant_marginal_lklhood(const gsl_vector *x, void *params);
#endif /* LIKELIHOOD_FUNCTIONS_H_ */
