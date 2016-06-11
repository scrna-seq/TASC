//
// Created by jiacheng on 11/12/15.
//

#include "NimMarginalLikelihood.h"
#include "NimArgs.h"
#include "../em_nim_algo/likelihood_functions.h"
#include "../em_nim_algo/adaptiveIntegrate.h"
#include "../misc/utils.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <limits>


double nim_marginal_lklhood(const gsl_vector *x, void *params){
    // x is a vector of size 2;
    // mu x[0]
    // log(sigma) x[1]

    NimArgs *nim_args = (NimArgs *) params;

    double sigma = std::exp(gsl_vector_get(x,1));

    if (std::isinf(sigma)){
        return 1.0/0.0;
    }


    double pars[9];
    void * ptr_pars = pars;

    pars[0] = 0;
    pars[5] = gsl_vector_get(x,0);
    pars[6] = sigma;
    pars[7] = 1;

//    timed_log("sigma:"+std::to_string((long double)sigma));
    gsl_function F_P;
    F_P.function = &single_lklhood;
    F_P.params = ptr_pars;
    size_t workspace_size = 1000;
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(workspace_size);
    double log_lklhood=0;

    for (size_t i=0; i<nim_args->sample_size; i++){
        pars[1] = nim_args->kappa.at(i);
        pars[2] = nim_args->tau.at(i);
        pars[3] = nim_args->alpha.at(i);
        pars[4] = nim_args->beta.at(i);
        pars[8] = nim_args->y.at(i);
        double integrated_result_p;
        double error;
        double epsrel = 1e-7;
        int status = adaptiveIntegrate(&F_P, pars, 0, epsrel, workspace, &integrated_result_p, &error);
        if (status!=0){
            timed_log("Exception: gsl_integration_cquad fails in nim marginal likelihood at iteration "+std::to_string((long long int)i)+".");
            return nan("");
        }
        log_lklhood += std::log(integrated_result_p);
    }
    gsl_integration_cquad_workspace_free(workspace);
    return -log_lklhood;
}
