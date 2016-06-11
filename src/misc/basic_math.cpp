/*
 * basic_math.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#include "basic_math.h"
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

bool isAllEqual(std::vector<size_t> sampleSizeArr){
	if (sampleSizeArr.empty()){
		return false;
	}
	size_t value = sampleSizeArr[0];
	for(size_t i=1; i<sampleSizeArr.size(); ++i){
		if(value != sampleSizeArr[i]){
			return false;
		}
	}
	return true;
}


double dpois(double x, double lambda){
	if(!std::isfinite(lambda)){
		return 0;
	}

	if (lambda==0){
		if(x==0){
			return 1;
		} else {
			return 0;
		}
	}
//	return std::pow(lambda, x)*std::exp(-lambda);
	boost::math::poisson_distribution<double> poiss_obj(lambda);
	return(boost::math::pdf(poiss_obj,x));
}

double l2norm(double arr[], size_t size){
	double sum = 0;
	for (size_t i=0; i<size; i++){
		sum += arr[i]*arr[i];
	}
	return sqrt(sum);
}

double l2norm(gsl_vector *x){
	double sum=0;
	for(size_t i=0;i<x->size;++i){
		sum += gsl_vector_get(x,i)*gsl_vector_get(x,i);
	}
	return sqrt(sum);
}

double gsl_vector_l2norm(gsl_vector *vec, size_t size){
	double arr[size];
	for (size_t i=0; i<size; i++){
		arr[i] = gsl_vector_get (vec, i);
	}
	return l2norm(arr, size);
}

std::vector<double> gsl_to_std_vector(gsl_vector *vec, size_t size){
	std::vector<double> std_vec;
	for (size_t i=0; i<size; i++){
		std_vec.push_back(gsl_vector_get(vec,i));
	}
	return std_vec;
}

std::vector<double> arr_to_std_vector(double vec[], size_t size){
	std::vector<double> std_vec;
	std_vec.assign(vec,vec+size);
	return std_vec;
}

gsl_matrix * remove_column(gsl_matrix *x, size_t i){
	gsl_matrix *y = gsl_matrix_alloc(x->size1,x->size2-1);
	for(size_t j=0; j<i-1; ++j){
		for(size_t k=0; k<x->size1; ++k){
			gsl_matrix_set(y,k,j,gsl_matrix_get(x,k,j));
		}
	}
	for(size_t j=i; j<x->size2; ++j){
		for(size_t k=0; k<x->size1; ++k){
			gsl_matrix_set(y,k,j-1,gsl_matrix_get(x,k,j));
		}
	}
	return y;
}

gsl_matrix * remove_row(gsl_matrix *x, size_t i){
	gsl_matrix *y = gsl_matrix_alloc(x->size1-1,x->size2);
	for(size_t j=0; j<i-1; ++j){
		for(size_t k=0; k<x->size2; ++k){
			gsl_matrix_set(y,j,k,gsl_matrix_get(x,j,k));
		}
	}
	for(size_t j=i; j<x->size1; ++j){
		for(size_t k=0; k<x->size2; ++k){
			gsl_matrix_set(y,j-1,k,gsl_matrix_get(x,j,k));
		}
	}
	return y;
}

bool std_vector_all_equal(std::vector<size_t> vec, size_t *val){

	if (vec.empty()){
		return false;
	}else{
		*val = *vec.begin();
		for (std::vector<size_t>::const_iterator i = vec.begin()+1; i!=vec.end(); ++i){
			if (*i!=*val){
				return false;
			}
			*val = *i;
		}
		return true;
	}
}

double dmvnorm(gsl_vector *x,  gsl_vector *mean, gsl_matrix *sigma){

	for (size_t i=0;i<x->size;i++){
		if (!std::isfinite(gsl_vector_get(x,i))){
			return 0;
		}
	}
	gsl_vector *x_minus_mean = gsl_vector_alloc(x->size);
	gsl_vector_memcpy(x_minus_mean,x);
	gsl_vector_sub(x_minus_mean,mean);

	gsl_matrix *s = gsl_matrix_alloc(sigma->size1,sigma->size2);
	gsl_matrix_memcpy(s,sigma);

	int status = gsl_linalg_cholesky_decomp(s);
	if (status) {
		return -1;
	}
	gsl_vector *InvSigmaXMinusMean = gsl_vector_alloc(x->size);
	gsl_linalg_cholesky_solve(s,x_minus_mean,InvSigmaXMinusMean);

	double  result = 0;
	gsl_blas_ddot(x_minus_mean,InvSigmaXMinusMean,&result);

	double detSigma = 1;
	for (size_t i = 0; i < s->size1; ++i){
		detSigma = detSigma*gsl_matrix_get(s,i,i);
	}
	double k = x->size;
	gsl_vector_free(x_minus_mean);
	gsl_matrix_free(s);
	gsl_vector_free(InvSigmaXMinusMean);
	return std::pow(2*pi,-k/2)/std::abs(detSigma)*std::exp(-result/2);
}

void print_vec(const gsl_vector *x, size_t n, std::string x_lab){
	std::cout << x_lab << ":\t";
	for(size_t i=0; i<std::min(n,x->size); ++i){
		std::cout << gsl_vector_get(x,i)<<",";
	}
	std::cout << std::endl;
}

void print_mat(gsl_matrix *x, size_t n1, size_t n2, std::string x_lab){
	std::cout << x_lab << ":\t";
	for(size_t i=0; i<std::min(n1,x->size1); ++i){
		for(size_t j=0; j<std::min(n2,x->size2); ++j){
		std::cout << gsl_matrix_get(x,i,j)<<",";
		}
		std::cout << std::endl;
	}
}

double expit(double t){
	return 1-1/(1+exp(t));
}

std::vector<size_t> array_minus(std::vector<size_t> all, std::vector<size_t> some){
	std::vector<size_t> result;
	for(std::vector<size_t>::iterator it=all.begin();it!=all.end();it++){
		if(std::find(some.begin(),some.end(),*it)==some.end()){
			result.push_back(*it);
		}
	}
	return result;
}

double transform_rho_from_real_domain(double real_rho) {
	return 1.0-2.0/(1.0+exp(real_rho));
}

double transform_rho_to_real_domain(double gunit_interval_rho) {
	if(gunit_interval_rho>1 || gunit_interval_rho <-1){
		return std::nan("");
	} else{
		return log((1+gunit_interval_rho)/(1-gunit_interval_rho));
	}
}

double dnorm(double x, double mean, double sd){
	boost::math::normal_distribution<double> norm_obj (mean, sd);
//	std::cout <<"dnorm:\t"<<mean<<"\t"<<sd<<std::endl;
	return(boost::math::pdf(norm_obj,x));
}

double pchisq(double x){
	boost::math::chi_squared_distribution<double> chi_sq (1);
	return boost::math::cdf(chi_sq,x);
}
