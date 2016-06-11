//
// Created by jiacheng on 11/12/15.
//

#ifndef TASC_NIMMARGINALLIKELIHOOD_H
#define TASC_NIMMARGINALLIKELIHOOD_H

#include <gsl/gsl_vector_double.h>

double nim_marginal_lklhood(const gsl_vector *x, void *params);

#endif //TASC_NIMMARGINALLIKELIHOOD_H
