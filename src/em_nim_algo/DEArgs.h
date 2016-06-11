/*
 * NIMArgs.h
 *
 *  Created on: Jul 19, 2015
 *      Author: cheng
 */

#ifndef NIMARGS_H_
#define NIMARGS_H_

#include <gsl/gsl_matrix.h>
#include <vector>
#include <memory>
class DEArgs {
public:
	DEArgs(std::vector<double> &kappa,
		   std::vector<double> &tau,
		   std::vector<double> &alpha,
		   std::vector<double> &beta,
		   std::vector<double> &scale_factors, std::vector<double> &y,
		   gsl_matrix *x, size_t sample_size, double epsabs, bool laplace);
	std::vector<double> kappa;
	std::vector<double> tau;
	std::vector<double> alpha;
	std::vector<double> beta;
	std::vector<double> scale_factors;
	std::vector<double> y;
	bool laplace;
	gsl_matrix *x;
	size_t sample_size;
	double epsabs;
	size_t max_iter_cquad;
	virtual ~DEArgs();
};

#endif /* NIMARGS_H_ */
