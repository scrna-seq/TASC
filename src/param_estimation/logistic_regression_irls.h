//
// Created by jiacheng on 11/25/15.
//

#ifndef NEWINTEGRATE_LOGISTIC_REGRESSION_IRLS_H
#define NEWINTEGRATE_LOGISTIC_REGRESSION_IRLS_H
#include <gsl/gsl_matrix.h>
#include <vector>
int logistic_regression(gsl_matrix *orig_x, gsl_vector *y, gsl_vector **beta, gsl_vector **yhat, size_t max_iter, double eps, bool *converge);
int fitBivariateNormal(std::vector<double> kappa, std::vector<double> tau, double *e_kappa, double *e_tau, double *sd_kappa, double *sd_tau, double *rho_kappa_tau);
#endif //NEWINTEGRATE_LOGISTIC_REGRESSION_IRLS_H
