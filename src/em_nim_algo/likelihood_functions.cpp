/*
 * likelihood_functions.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#include <gsl/gsl_integration.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>


#include "../misc/utils.h"
#include "../misc/basic_math.h"
#include "DEArgs.h"
#include "likelihood_functions.h"
#include "adaptiveIntegrate.h"
#include "../misc/cheng_quad.h"
#include "../bio_var_computer/NimArgs.h"
#include <cmath>

double single_lklhood (double x, void * params){
	double* pars = (double *) params;
	/*pars[0] = phi;
			pars[1] = kappa;
			pars[2] = tau;
			pars[3] = alpha;
			pars[4] = beta;
			pars[5] = theta_t[i];
			pars[6] = sigma_t;
			pars[7] = scaleFactors[i];
			pars[8] = y.at(i);*/
	double kappa = *(pars + 1);
	double tau = *(pars + 2);
	double alpha = *(pars + 3);
	double beta = *(pars + 4);
	double theta = *(pars + 5);
	double sigma = *(pars + 6);
	double y = *(pars + 8);

	//dpois(y,lambda) = true_dpois(y,lambda)*factorial(y)
	//this is to prevent underflow
	//this does not change the optimization results because it inflates the results by prod_i(factorial(y_i)
	//and y_i are constants in our optimization process.

	//the reason it's equivalent to multiplying a constant is because when y!=0 it's multiplying a constant
	//when y==0, dpois(0,lambda) = true_dpois(0,lambda)*factorial(0) == true_dpois(0,lambda)
	//so although there is addition when y==0, in this case dpois(0,lambda)==true_dpois(0,lambda) so it works out.

	double dpois_this = dpois(y,exp(alpha+beta*x));
	double base = (1.0-1.0/(1.0+exp(kappa+tau*x)))*dpois_this;
	if (y==0) {
		base = base + 1.0/(1.0+exp(kappa+tau*x));
	}
	if (std::isnan(x)){
		return 1.0/0.0;
	}
	double dnorm_this = dnorm(x,theta,sigma);
	base = base*dnorm_this;
//	std::cout << "x,y,alpha,beta,theta,sigma,base,dnorm,dpois,scale_factor:"<<x<<"\ty:"
//				<<y<<"\talpha:"<<alpha<<"\t"<<beta<<"\t"<<theta<<"\t"<<sigma<<"\t"<<base<<"\t"
//		<<dnorm_this<<"\t"<<dpois_this<<"\t"<<scale_factor<<std::endl;
	return base;
}

double ew_integrand (double x, void * params){
	double* pars = (double *) params;
	double kappa = *(pars + 1);
	double tau = *(pars + 2);
	double alpha = *(pars + 3);
	double beta = *(pars + 4);
	double theta = *(pars + 5);
	double sigma = *(pars + 6);
	double y = *(pars + 8);
	double base = (1.0-1.0/(1.0+exp(kappa+tau*x)))*dpois(y,exp(alpha+beta*x));
	if (y==0) {
		base = base + 1.0/(1.0+exp(kappa+tau*x));
	}
	base = x*base*dnorm(x,theta,sigma);
	return base;
}

double ew_sq_integrand (double x, void * params){
	double* pars = (double *) params;
	double kappa = *(pars + 1);
	double tau = *(pars + 2);
	double alpha = *(pars + 3);
	double beta = *(pars + 4);
	double theta = *(pars + 5);
	double sigma = *(pars + 6);

	double y = *(pars + 8);
	double base = (1.0-1.0/(1.0+exp(kappa+tau*x)))*dpois(y,exp(alpha+beta*x));
	if (y==0) {
		base = base + 1.0/(1.0+exp(kappa+tau*x));
	}
	base = (x-theta)*(x-theta)*base*dnorm(x,theta,sigma);
	return base;
}

double em_marginal_lklhood(std::vector<double> kappa,
						   std::vector<double> tau,
						   std::vector<double> alpha,
						   std::vector<double> beta,
						   gsl_vector *theta,
						   double sigma,
						   std::vector<double> scale_factors,
						   std::vector<double> y, size_t sample_size, double epsabs) {

	double pars[9];
	void * ptr_pars = pars;
	pars[0] = 0;

	pars[6] = sigma;
	gsl_function F_P;
	F_P.function = &single_lklhood;
	F_P.params = ptr_pars;
	size_t workspace_size = 1000;
	gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(workspace_size);
	double log_lklhood=0;
	for (size_t i=0; i<sample_size; i++){
		pars[1] = kappa[i];
		pars[2] = tau[i];
		pars[3] = alpha[i];
		pars[4] = beta[i];
		pars[5] = gsl_vector_get(theta,i);
		pars[7] = scale_factors[i];
		pars[8] = y[i];
		double integrated_result_p;
		double error;
		int status = adaptiveIntegrate(&F_P, pars, epsabs,epsabs, workspace, &integrated_result_p, &error);
		if (status!=0){
			timed_log("Exception: gsl_integration_qagi fails in em_marginal_lklhood at iteration "+std::to_string((long long int)i)+".");
			throw 3;
		}
		log_lklhood += log(integrated_result_p);
	}
	gsl_integration_cquad_workspace_free(workspace);
	return log_lklhood;
}

double secondOrderDeriv(double a, double b, double k, double t, double sigma, double x, bool yIsZero){

	double cg = 0;
	if (yIsZero){
		cg = -(0.2e1 + std::exp((double) (3 * t * x) + (double) (3 * k) -0.2e1 * std::exp((double) (b * x + a))) * sigma * sigma * (double) (t* t) - std::exp((double) (3 * t * x) + (double) (3 * k) -std::exp((double) (b * x + a))) * sigma * sigma * (double) (t * t) +std::exp((double) (2 * t * x) + (double) (2 * k) - 0.2e1 *std::exp((double) (b * x + a)) + (double) (b * x) + (double) a) *(double) (b * b) * sigma * sigma + std::exp((double) (3 * t * x) +(double) (3 * k) - std::exp((double) (b * x + a)) + (double) (b * x) +(double) a) * (double) (b * b) * sigma * sigma - std::exp((double) (2* b * x) + (double) (2 * a) + (double) (t * x) - std::exp((double) (b* x + a)) + (double) k) * (double) (b * b) * sigma * sigma + 0.4e1 *std::exp((double) (b * x) + (double) (2 * t * x) - std::exp((double)(b * x + a)) + (double) a + (double) (2 * k)) * (double) b * sigma *sigma * (double) t + 0.2e1 * std::exp((double) (b * x) + (double) (t *x) - std::exp((double) (b * x + a)) + (double) a + (double) k) *(double) b * sigma * sigma * (double) t + 0.8e1 * std::exp((double) (2* t * x) + (double) (2 * k) - std::exp((double) (b * x + a))) + 0.4e1* std::exp((double) (t * x + k)) + 0.2e1 * std::exp((double) (4 * t *x) + (double) (4 * k) - 0.2e1 * std::exp((double) (b * x + a))) -std::exp((double) (t * x) + (double) k - std::exp((double) (b * x +a))) * sigma * sigma * (double) (t * t) + 0.2e1 * std::exp((double) (b* x) + (double) (2 * t * x) - std::exp((double) (b * x + a)) +(double) a + (double) (2 * k)) * (double) (b * b) * sigma * sigma +std::exp((double) (b * x) + (double) (t * x) - std::exp((double) (b *x + a)) + (double) a + (double) k) * (double) (b * b) * sigma * sigma+ 0.4e1 * std::exp((double) (3 * t * x) + (double) (3 * k) -std::exp((double) (b * x + a))) + 0.4e1 * std::exp((double) (3 * t *x) + (double) (3 * k) - 0.2e1 * std::exp((double) (b * x + a))) +0.2e1 * std::exp((double) (2 * t * x) + (double) (2 * k) - 0.2e1 *std::exp((double) (b * x + a))) + std::exp((double) (4 * t * x) +(double) (4 * k) - 0.2e1 * std::exp((double) (b * x + a)) + (double)(b * x) + (double) a) * (double) (b * b) * sigma * sigma + 0.2e1 *std::exp((double) (3 * t * x) + (double) (3 * k) - 0.2e1 *std::exp((double) (b * x + a)) + (double) (b * x) + (double) a) *(double) (b * b) * sigma * sigma - std::exp((double) (2 * b * x) +(double) (3 * t * x) + (double) (2 * a) + (double) (3 * k) -std::exp((double) (b * x + a))) * (double) (b * b) * sigma * sigma -0.2e1 * std::exp((double) (2 * b * x) + (double) (2 * t * x) +(double) (2 * a) + (double) (2 * k) - std::exp((double) (b * x + a)))* (double) (b * b) * sigma * sigma + 0.2e1 * std::exp((double) (3 * t* x) + (double) (3 * k) - std::exp((double) (b * x + a)) + (double) (b* x) + (double) a) * (double) b * sigma * sigma * (double) t + 0.2e1 *std::exp((double) (2 * t * x + 2 * k)) + std::exp((double) (t * x +k)) * sigma * sigma * (double) (t * t) + 0.4e1 * std::exp((double) (t* x) + (double) k - std::exp((double) (b * x + a)))) *std::pow(std::exp((double) (t * x) + (double) k - std::exp((double) (b* x + a))) + 0.1e1, -0.2e1) * std::pow(sigma, -0.2e1) * std::pow(0.1e1+ std::exp((double) (t * x + k)), -0.2e1);
	} else {
		cg = -(std::exp((double) (b * x + 2 * t * x + a + 2 * k)) * (double)(b * b) * sigma * sigma + 0.2e1 * (double) (b * b) * std::exp((double)(b * x + t * x + a + k)) * sigma * sigma + std::exp((double) (t * x +k)) * sigma * sigma * (double) (t * t) + (double) (b * b) *std::exp((double) (b * x + a)) * sigma * sigma + 0.2e1 *std::exp((double) (2 * t * x + 2 * k)) + 0.4e1 * std::exp((double) (t* x + k)) + 0.2e1) * std::pow(0.1e1 + std::exp((double) (t * x + k)),-0.2e1) * std::pow(sigma, -0.2e1);
	}
	return cg;
}

double logDPois(double k, double log_lambda){
	return k*log_lambda-std::lgamma((long int)(k+1.5))-std::exp(log_lambda);
}

double logExpit(double x){
	if ((1.0-1.0/(1+std::exp(x)))==0){
		return x;
	} else {
		return std::log(1.0-1.0/(1+std::exp(x)));
	}
}

double logDNorm(double x, double mu, double sigma){
	if (sigma==0){
		if (x==mu){
			return 0;
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	} else {
		return -std::log(sigma*sigma*2*pi)/2-(x-mu)*(x-mu)/2/sigma/sigma;
	}
}

double logDpois0(double log_mean){
	return -exp(log_mean);
}

double max2(double a, double b){
    double max = a;
    if (b>a){
        max = b;
    }
    return max;
}

double log_sum_exp(double a, double b){
    double max_el = max2(a,b);
    return max_el + std::log(std::exp(a - max_el) + std::exp(b - max_el));
}

double negLogLikelihood(const gsl_vector *t, void *params){
//    std::cout<<x<<"\t";
	double x = gsl_vector_get(t,0);
	double* pars = (double *) params;
	double kappa = *(pars + 1);
	double tau = *(pars + 2);
	double alpha = *(pars + 3);
	double beta = *(pars + 4);
	double theta = *(pars + 5);
	double sigma = *(pars + 6);
	double y = *(pars + 8);

    double result;

	if (y==0) {
        result = log_sum_exp(logExpit(-kappa-tau*x), logExpit(kappa+tau*x) + logDpois0(alpha+beta*x));
	} else {
        result = logExpit(kappa+tau*x) + logDPois(y, alpha+beta*x);
    }

    return -(result+logDNorm(x,theta, sigma));
}


double absSingleLikelihood(double x, void *params){
	double* pars = (double *) params;
	/*pars[0] = phi;
            pars[1] = kappa;
            pars[2] = tau;
            pars[3] = alpha;
            pars[4] = beta;
            pars[5] = theta_t[i];
            pars[6] = sigma_t;
            pars[7] = scaleFactors[i];
            pars[8] = y.at(i);*/
	double scale_factor = *(pars + 7);
	//std::cout << "scale_factor:\t"<<scale_factor << std::endl;

	gsl_vector *ix = gsl_vector_alloc(1);
	gsl_vector_set(ix,0,x);

	double result = std::exp(-negLogLikelihood(ix, params) + scale_factor);

	gsl_vector_free(ix);
	return result;
}


double newIntegrate(double (*absLikelihood)(double, void *), void *params, double epsrel, double scaleFactor, double lower, double upper, size_t max_iter_cquad){
	double *pars = (double *) params;


	*(pars+7) = scaleFactor;
	gsl_function F_P;
	F_P.function = absLikelihood;
	F_P.params = params;
	size_t workspace_size = 1000;
	gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(workspace_size);
	double integrated_result_p;
	double error;

    size_t neval = max_iter_cquad;
	int status = gsl_integration_cquad_cheng(&F_P, lower, upper,0,epsrel,workspace,&integrated_result_p,&error,&neval);

	if (neval==max_iter_cquad || status!=0){
		timed_log("Exception: gsl_integration_cquad fails in de marignal likelihood.");
		return nan("");
	}

	gsl_integration_cquad_workspace_free(workspace);
	return integrated_result_p;
}


double margSingleLikelihood(double (*fn1)(const gsl_vector *, void *), void *params, double epsrel, bool laplace, size_t max_iter_cquad) {

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *ix;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	double *pars = (double *) params;


	/* Starting point */
	ix = gsl_vector_alloc (1);

	gsl_vector_set_all(ix, 1.0);

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (1);

	if((*(pars+8))==0){
		gsl_vector_set_all (ss, *(pars+5));
	} else {
		gsl_vector_set_all (ss, std::log(*(pars+8)));
	}

    /* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = fn1;
	minex_func.params = (void *)params;
    s = gsl_multimin_fminimizer_alloc (T, 1);
    gsl_multimin_fminimizer_set (s, &minex_func, ix, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-7);

	}while (status == GSL_CONTINUE && iter < 500);
//	bool converged = false;
//	if ( status == GSL_SUCCESS || status == 0 ) {
//		converged = true;
//	}
//    std::cout<<("convergence: "+std::to_string(converged))<<std::endl;
	double argmaxX = gsl_vector_get(s->x,0);
	double maxLogLklhood = s->fval;
	double secondDeriv = secondOrderDeriv(*(pars+3),*(pars+4),*(pars+1),*(pars+2),*(pars+6),argmaxX, (*(pars+8)==0));
//    std::cout << "argmaxX:\t" << argmaxX << "\t";
//    std::cout << "maxLogLklhood:\t" << maxLogLklhood << "\t";
//    std::cout << "secondDeriv:\t" << secondDeriv << "\t";

//    std::cout << "lowerBound, -20sigma:\t"<<lowerB<<"\t";
//    std::cout << "upperBound, +20sigma:\t"<<upperB<<"\t";

	double result = 0;

	if (laplace) {
		result = std::log(2 * pi) / 2 - std::log(std::abs(secondDeriv)) / 2 - maxLogLklhood;
	} else {
		double lowerB = argmaxX - 20 / std::sqrt(std::abs(secondDeriv));
		double upperB = argmaxX + 20 / std::sqrt(std::abs(secondDeriv));

		result = std::log(newIntegrate(&absSingleLikelihood, params, epsrel, maxLogLklhood, lowerB, upperB, max_iter_cquad)) -
				 maxLogLklhood;

	}
//	double laplacian = std::log(2*pi)/2-std::log(std::abs(secondDeriv))/2-maxLogLklhood;
//    std::cout << "Laplacian Approx:\t" << laplacian << "\t";

//    std::cout << "newIntegrate:\t"<<cquadmarginal<<std::endl;


	gsl_vector_free(ix);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return result;
}



double de_marginal_lklhood(const gsl_vector *x, void *params){
	DEArgs *nim_args = (DEArgs *) params;
	gsl_vector *theta = gsl_vector_alloc(nim_args->sample_size);
	gsl_vector *beta = gsl_vector_alloc(x->size-1);

	for (size_t i = 0; i<x->size-1; ++i){
		gsl_vector_set(beta,i,gsl_vector_get(x,i));
	}

	double sigma = std::exp(gsl_vector_get(x,x->size-1));
	int status = gsl_blas_dgemv (CblasNoTrans, 1.0, nim_args->x, beta, 0.0, theta);
	if (status) {
		return nan("");
	}
	double pars[9];
	pars[0] = 0;
	pars[6] = sigma;
	double log_lklhood=0;

	for (size_t i=0; i<nim_args->sample_size; i++){
		//verbose_timed_log("inside de_marginal_lklhood, sample: "+std::to_string(i));

		pars[1] = nim_args->kappa.at(i);
		pars[2] = nim_args->tau.at(i);
		pars[3] = nim_args->alpha.at(i);
		pars[4] = nim_args->beta.at(i);

		pars[5] = gsl_vector_get(theta,i);
		pars[7] = 1;
		pars[8] = nim_args->y.at(i);

		log_lklhood += margSingleLikelihood(&negLogLikelihood, pars, nim_args->epsabs, nim_args->laplace, nim_args->max_iter_cquad);
	}
	gsl_vector_free(theta);
	gsl_vector_free(beta);
	return -log_lklhood;
}

double quant_marginal_lklhood(const gsl_vector *x, void *params){
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

        log_lklhood += margSingleLikelihood(&negLogLikelihood, pars, nim_args->epsabs, false, 2000);
    }
    gsl_integration_cquad_workspace_free(workspace);
    return -log_lklhood;
}

