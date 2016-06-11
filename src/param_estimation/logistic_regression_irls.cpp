#include "logistic_regression_irls.h"
#include "../misc/basic_math.h"
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_statistics.h>

/*
 * This file contains the functions
 * to perform logistic regressions
 * using IRLS (iteratively reweighted
 * least squares)
 * */
inline double linkinv(double t){
    return 1-1/(1+exp(t));
}

inline double gprime(double t){
    return 1/(exp(t/2)+exp(-t/2))/(exp(t/2)+exp(-t/2));
}

inline double gprime_sq(double t){
    return gprime(t)*gprime(t);
}

inline double varg(double t){
    return t*(1-t);
}


/*
 * function to perform logistic regression
 * orig_x: pointer to the gsl_matrix of dimension nxp
 * y: pointer to the gsl_vector of dimension nx1
 * beta: pointer to a pointer of gsl_vector, storing covariates
 * yhat: pointer to a pointer of gsl_vector, storing hat(y)
 * eps: max_tol for l2-norm of s_old - s_new
 * converge: pointer to a boolean, storing whether the algorith is converged
 * return 1 if fail, return 0 if success
 * */
int logistic_regression(gsl_matrix *orig_x, gsl_vector *y, gsl_vector **beta, gsl_vector **yhat, size_t max_iter, double eps, bool *converge){

    gsl_set_error_handler_off();
    size_t p = orig_x->size2;
    size_t n = orig_x->size1;
    if (n!=y->size){
        return 1;
    }

    gsl_matrix *x = gsl_matrix_alloc(n,p);
    gsl_matrix_memcpy(x,orig_x);

    *beta = gsl_vector_alloc(p);
    *yhat = gsl_vector_alloc(n);

    gsl_vector *t = gsl_vector_alloc(n);
    gsl_vector_set_zero(t);

    gsl_vector *work_vec = gsl_vector_alloc(p);
    gsl_matrix *work_mat = gsl_matrix_alloc(p,p);
    gsl_matrix *v = gsl_matrix_alloc(p,p);
    gsl_vector *d = gsl_vector_alloc(p);
    gsl_linalg_SV_decomp_mod(x,work_mat,v,d,work_vec);

    gsl_vector *z = gsl_vector_alloc(n);
    gsl_matrix *W = gsl_matrix_alloc(n,n);
    gsl_matrix_set_zero(W);
    gsl_vector *s = gsl_vector_alloc(p);
    gsl_vector_set_zero(s);
    gsl_vector *s_old = gsl_vector_alloc(p);
    gsl_matrix *UWU = gsl_matrix_alloc(p,p);
    gsl_matrix *WU = gsl_matrix_alloc(n,p);
    gsl_vector *Wz = gsl_vector_alloc(n);
    gsl_vector *UWz = gsl_vector_alloc(p);
    gsl_matrix *SVD_V = gsl_matrix_alloc(p,p);
    gsl_vector *SVD_S = gsl_vector_alloc(p);
    for(size_t step=0; step < max_iter; ++step){
        double ti;
        for(size_t i=0;i<n;++i){
            ti = gsl_vector_get(t,i);
            gsl_vector_set(z,i,ti+(gsl_vector_get(y,i)-linkinv(ti))/gprime(ti));
            gsl_matrix_set(W,i,i,gprime_sq(ti)/varg(linkinv(ti)));
        }
        gsl_vector_memcpy(s_old,s);
        gsl_matrix_set_zero(WU);
        gsl_matrix_set_zero(UWU);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,W,x,0.0,WU);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,x,WU,0.0,UWU);
        int chol_s = gsl_linalg_cholesky_decomp(UWU);
        gsl_blas_dgemv(CblasNoTrans,1.0,W,z,0.0,Wz);
        gsl_blas_dgemv(CblasTrans,1.0,x,Wz,0.0,UWz);

        if(chol_s){
            gsl_linalg_SV_decomp_jacobi (UWU,SVD_V,SVD_S);
            gsl_linalg_SV_solve(UWU, SVD_V, SVD_S, UWz, s);
        }else{
            gsl_linalg_cholesky_solve(UWU,UWz,s);
        }
        gsl_blas_dgemv(CblasNoTrans,1.0,x,s,0.0,t);
        gsl_vector_sub (s_old, s);
        if (l2norm(s_old)<eps){
            break;
        }
    }

    gsl_matrix_free(work_mat);
    gsl_matrix_free(W);
    gsl_matrix_free(UWU);
    gsl_matrix_free(WU);
    gsl_matrix_free(SVD_V);
    gsl_vector_free(work_vec);
    gsl_vector_free(s);
    gsl_vector_free(z);
    gsl_vector_free(s_old);
    gsl_vector_free(Wz);
    gsl_vector_free(UWz);
    gsl_vector_free(SVD_S);


    int status = gsl_linalg_SV_solve(x,v,d,t,*beta);
    gsl_matrix_free(x);
    gsl_matrix_free(v);
    gsl_vector_free(d);
    gsl_vector_free(t);

    if(status){
        return 1;
    }
    gsl_blas_dgemv(CblasNoTrans,1.0,orig_x,*beta,0.0,*yhat);
    for (size_t i=0;i<(*yhat)->size;++i){
        gsl_vector_set(*yhat,i,linkinv(gsl_vector_get(*yhat,i)));
    }
    bool no_nan = true;
    for (size_t j = 0; j < (*beta)->size; ++j) {
        if (std::isnan(gsl_vector_get((*beta),j))){
            no_nan = false;
            break;
        }
    }

    for (size_t j = 0; j < (*yhat)->size; ++j) {
        if (std::isnan(gsl_vector_get((*yhat),j))){
            no_nan = false;
            break;
        }
    }

    if (no_nan) {
        return 0;
    } else {
        return 1;
    }
}

/*
 * function to fit a bivariate normal
 * distribution.
 * kappa: vector<double> for kappa
 * tau: vector<double> for tau
 * e_kappa: pointer to a double, storing EKappa
 * e_tau: pointer to a double, storing ETau
 * sd_kappa: pointer to a double, storing SDKappa
 * sd_tau: pointer to a double, storing SDTau
 * rho_kappa_tau: pointer to a double, storing
 * correlation b/w kappa and tau
 * */
int fitBivariateNormal(std::vector<double> kappa, std::vector<double> tau, double *e_kappa, double *e_tau, double *sd_kappa, double *sd_tau, double *rho_kappa_tau){
    // assume equal length of kappa and tau. no error checking if not.
    gsl_matrix *kt_sample = gsl_matrix_alloc(kappa.size(), 2);
    gsl_vector *k_vec = gsl_vector_alloc(kappa.size());
    gsl_vector *t_vec = gsl_vector_alloc(tau.size());

    for (size_t i = 0; i < kappa.size(); ++i) {
        gsl_matrix_set(kt_sample, i, 0, kappa[i]);
        gsl_matrix_set(kt_sample, i, 1, tau[i]);
        gsl_vector_set(k_vec,i,kappa[i]);
        gsl_vector_set(t_vec,i,tau[i]);
    }

    (*e_kappa) = gsl_stats_mean(k_vec->data, 1, k_vec->size);
    (*e_tau) = gsl_stats_mean(t_vec->data,1,t_vec->size);



    gsl_matrix *kt_mean_vec = gsl_matrix_alloc(1,2);
    gsl_matrix_set(kt_mean_vec,0,0,*e_kappa);
    gsl_matrix_set(kt_mean_vec,0,1,*e_tau);

    gsl_matrix *ones = gsl_matrix_alloc(kappa.size(),1);
    gsl_matrix_set_all(ones,1);

    gsl_matrix *diff = gsl_matrix_alloc(kappa.size(),2);
    gsl_matrix_memcpy(diff,kt_sample);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,ones,kt_mean_vec,1.0,diff);

    gsl_matrix *sigma = gsl_matrix_alloc(2,2);
    gsl_matrix_set_all(sigma,0);

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0/(kappa.size()-1),diff,diff,0.0,sigma);

    (*sd_kappa)=std::sqrt(gsl_matrix_get(sigma,0,0));
    (*sd_tau)=std::sqrt(gsl_matrix_get(sigma,1,1));
    (*rho_kappa_tau)=gsl_matrix_get(sigma,0,1)/(*sd_kappa)/(*sd_tau);


    gsl_matrix_free(diff);
    gsl_matrix_free(kt_mean_vec);
    gsl_matrix_free(kt_sample);
    gsl_matrix_free(ones);
    gsl_matrix_free(sigma);
    gsl_vector_free(k_vec);
    gsl_vector_free(t_vec);
    return 0;
}