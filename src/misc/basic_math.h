/*
 * basic_math.h
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#ifndef BASIC_MATH_H_
#define BASIC_MATH_H_
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <vector>
#include <string>
struct mvnorm_params{
	double e_kappa;
	double e_tau;
	double sd_kappa;
	double sd_tau;
	double rho_kappa_tau;
};
double dnorm(double x, double mean, double sd);
double pchisq(double x);
double mean(double arr[], size_t size);
double l2norm(double arr[], size_t size);
double gsl_vector_l2norm(gsl_vector *vec, size_t size);
std::vector<double> gsl_to_std_vector(gsl_vector *vec, size_t size);
std::vector<double> arr_to_std_vector(double vec[], size_t size);
gsl_matrix * remove_column(gsl_matrix *x,size_t i);
bool std_vector_all_equal(std::vector<size_t> vec, size_t *val);
double dmvnorm(gsl_vector *x,  gsl_vector *mean, gsl_matrix *sigma);
double l2norm(gsl_vector *x);
void print_mat(gsl_matrix *x, size_t n1, size_t n2, std::string x_lab);
void print_vec(const gsl_vector *x, size_t n, std::string x_lab);
const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036;
double expit(double t);
std::vector<size_t> array_minus(std::vector<size_t> all, std::vector<size_t> some);
template<typename TYPE>
std::vector<std::vector<TYPE> > transpose(std::vector<std::vector<TYPE> > x){

	std::vector<std::vector<TYPE> > t_x(x[0].size(), std::vector<TYPE>(x.size()));
	for (size_t i = 0; i < x[0].size(); i++){
		for (size_t j = 0; j < x.size(); j++) {
			t_x[i][j] = x[j][i];
		}
	}
	return t_x;
};
double dpois(double x, double lambda);
bool isAllEqual(std::vector<size_t> sampleSizeArr);
double transform_rho_from_real_domain(double real_rho);
double transform_rho_to_real_domain(double gunit_interval_rho);

template<typename TYPE>
double mean(std::vector<TYPE> arr){
	double sum = 0;
	for (size_t i=0; i<arr.size(); i++){
		sum += arr[i];
	}
	return sum/(double)arr.size();
}

#endif /* BASIC_MATH_H_ */
