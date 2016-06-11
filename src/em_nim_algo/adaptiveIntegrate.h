/*
 * adaptiveIntegrate.h
 *
 *  Created on: Nov 10, 2015
 *      Author: jiacheng
 */

#ifndef EM_NIM_ALGO_ADAPTIVEINTEGRATE_H_
#define EM_NIM_ALGO_ADAPTIVEINTEGRATE_H_

#include <gsl/gsl_integration.h>
int adaptiveIntegrate(gsl_function *F, double *pars, double epsabs, double epsrel, gsl_integration_cquad_workspace *workspace, double *result, double *error);


#endif /* EM_NIM_ALGO_ADAPTIVEINTEGRATE_H_ */

