/*
 * init_integrand.h
 *
 *  Created on: Aug 6, 2015
 *      Author: cheng
 */

#ifndef INIT_INTEGRAND_H_
#define INIT_INTEGRAND_H_
#include <gsl/gsl_matrix.h>
#include "../misc/basic_math.h"


double init_marginal_lklhood (const gsl_vector *x, void * params);
int pdf_single_cell(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int kappa_times_pdf_single_cell(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int tau_times_pdf_single_cell(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
#endif /* INIT_INTEGRAND_H_ */
