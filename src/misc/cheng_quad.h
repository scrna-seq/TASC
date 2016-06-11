//
// Created by jiacheng on 12/7/15.
//

#include <gsl/gsl_integration.h>
int gsl_integration_cquad_cheng(const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace *ws, double *result, double *abserr, unsigned long *nevals);
