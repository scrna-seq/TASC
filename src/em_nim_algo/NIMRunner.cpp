/*
 * NIMRunner.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: jiacheng
 */

#include <gsl/gsl_multimin.h>
#include <unistd.h>
#include <algorithm>
#include <sys/syscall.h>


#include "NIMRunner.h"
#include "../misc/utils.h"
#include "DEArgs.h"
#include "likelihood_functions.h"
#include "../misc/basic_math.h"


NIMRunner::NIMRunner() {
	sampleSize=0;
    numCovariates = 0;
    designMatrix = NULL;
    maxIter = 1000;
    initialSigma = 0;
    minRelTol = 1e-5;
    epsabs = 1e-7;
    laplace = false;
    maxIterCQUAD = 2000;
}

void NIMRunner::runModel() {

	this->nimResult.covariate_names = this->covariateNames;
	this->nimResult.gene_name = this->geneName;
	this->nimResult.test_cols = this->covariatesToBeTested;
	gsl_matrix *x = this->designMatrix;
	const size_t num_covariates = this->numCovariates;
	std::vector<double> y = this->geneCounts;
	size_t sample_size = this->sampleSize;

	if (this->initialSigma <=0 ){
		timed_log("Exception: initial value of sigma_t <= 0!");
		exit(1);
	}

	DEArgs nim_args(kappa, tau, alpha, beta, this->scaleFactors,
					y, x, sample_size, this->minRelTol/1000, this->laplace);
    nim_args.max_iter_cquad = this->maxIterCQUAD;
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *ix;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	ix = gsl_vector_alloc (num_covariates+1);
	gsl_vector_set_all(ix, 1.0);
	gsl_vector_set(ix, num_covariates, log(this->initialSigma));

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (num_covariates+1);
	gsl_vector_set_all (ss, 1.0);

	/* Initialize method and iterate */
	minex_func.n = num_covariates+1;
	minex_func.f = &de_marginal_lklhood;
	minex_func.params = (void *)&nim_args;

	s = gsl_multimin_fminimizer_alloc (T, num_covariates+1);
	gsl_multimin_fminimizer_set (s, &minex_func, ix, ss);

	do
	{
		verbose_timed_log("Full Model For Gene " + this->geneName + " Step " + std::to_string((long long int)iter) + " PID: " + std::to_string((long long int)syscall(SYS_gettid)));
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, epsabs);

	} while (status == GSL_CONTINUE && iter < this->maxIter);
	if ( status == GSL_SUCCESS || status == 0 ) {
		this->nimResult.full_converge=true;
	} else {
		this->nimResult.full_converge=false;
	}
	std::vector<double> beta_hat = gsl_to_std_vector(s->x,num_covariates+1);
	beta_hat.pop_back();
	this->nimResult.beta_hat = beta_hat;
	this->nimResult.full_sigma = gsl_vector_get(s->x,num_covariates);
	this->nimResult.full_lklhood = -s->fval;

	gsl_vector_free(ix);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);


	//reduced model
	ix = gsl_vector_alloc (num_covariates);
	ss = gsl_vector_alloc (num_covariates);
	s = gsl_multimin_fminimizer_alloc (T, num_covariates);

	for (size_t test_col_idx=0; test_col_idx<this->covariatesToBeTested.size();++test_col_idx){
		size_t test_col = this->covariatesToBeTested.at(test_col_idx);
		gsl_matrix *reduced_x = remove_column(x,test_col);
		std::string covariate_name_being_tested = this->nimResult.covariate_names.at(test_col-1);
		iter = 0;
		status = 0;
		size = 0;

		/* Starting point */

		gsl_vector_set_all(ix, 1.0);
		gsl_vector_set(ix, num_covariates-1, log(this->initialSigma));

		/* Set initial step sizes to 1 */

		gsl_vector_set_all (ss, 1.0);

		/* Initialize method and iterate */
		nim_args.x=reduced_x;
		minex_func.n = num_covariates;
		minex_func.f = &de_marginal_lklhood;
		minex_func.params = (void *)&nim_args;
		gsl_multimin_fminimizer_set (s, &minex_func, ix, ss);

		do
		{
			verbose_timed_log("Reduced Model For Gene " + this->geneName + " Step " + std::to_string((long long int)iter) + " PID: " + std::to_string((long long int)syscall(SYS_gettid)));
			iter++;

			status = gsl_multimin_fminimizer_iterate(s);

			if (status)
				break;

			size = gsl_multimin_fminimizer_size (s);

			status = gsl_multimin_test_size(size, this->minRelTol);


		} while (status == GSL_CONTINUE && iter < this->maxIter);

		std::pair<std::string, bool> converge;

		converge.first = covariate_name_being_tested;

		if ( status == GSL_SUCCESS || status == 0 ) {
			converge.second = true;
			this->nimResult.reduced_converge.push_back(converge);
		} else {
			converge.second = false;
			this->nimResult.reduced_converge.push_back(converge);
		}

		this->nimResult.reduced_sigma.push_back(std::pair<std::string,double>(covariate_name_being_tested,gsl_vector_get(s->x,num_covariates-1)));

		this->nimResult.reduced_lklhood.push_back(std::pair<std::string,double>(covariate_name_being_tested,-s->fval));

		double lrt_stat_for_this_covariate = 2 * (this->nimResult.full_lklhood - (-s->fval));

		if (lrt_stat_for_this_covariate<0) lrt_stat_for_this_covariate=0;

		double lrt_pval_for_this_covariate = 1-pchisq(lrt_stat_for_this_covariate);

		this->nimResult.lrt_stat.push_back(std::pair<std::string,double>(covariate_name_being_tested,lrt_stat_for_this_covariate));

		this->nimResult.lrt_pval.push_back(std::pair<std::string,double>(covariate_name_being_tested,lrt_pval_for_this_covariate));

		gsl_matrix_free(reduced_x);

	}

	gsl_vector_free(ix);

	gsl_vector_free(ss);

	gsl_multimin_fminimizer_free (s);

	this->nimResult.aborted=false;

}

void NIMRunner::setData(
		std::pair<std::vector<std::string>, gsl_matrix*> designMatrix,
		std::pair<std::string, std::vector<double> > geneNameAndCounts,
		std::vector<size_t> covariatesToBeTested, double initSigma,
		std::vector<double> alpha,
		std::vector<double> beta,
		std::vector<double> kappa,
		std::vector<double> tau,
		double maxIter,
		double minRelTol,
		bool laplace,
        size_t max_iter_cquad) {
	this->laplace = laplace;
	this->alpha = alpha;
	this->beta = beta;
	this->kappa = kappa;
	this->tau = tau;
	this->maxIter = maxIter;
	this->minRelTol = minRelTol;
	this->covariatesToBeTested = covariatesToBeTested;
	this->geneName = geneNameAndCounts.first;
	this->geneCounts = geneNameAndCounts.second;
	this->covariateNames = designMatrix.first;
	this->numCovariates = this->covariateNames.size();
	this->designMatrix = designMatrix.second;
	this->initialSigma = initSigma;
	std::vector<size_t> sampleSizeArr(6);
    sampleSizeArr[0] = this->designMatrix->size1;
    sampleSizeArr[1] = this->alpha.size();
    sampleSizeArr[2] = this->beta.size();
    sampleSizeArr[3] = this->kappa.size();
    sampleSizeArr[4] = this->tau.size();
    sampleSizeArr[5] = this->geneCounts.size();

	if(!isAllEqual(sampleSizeArr)){

		timed_log("Uneven sample size from X and Y. Now exiting.");
		exit(1);
	} else {
		this->sampleSize = this->geneCounts.size();
	}
    std::vector<double> sf(1);
    this->scaleFactors = sf;
    this->maxIterCQUAD = max_iter_cquad;
}

const DEResult & NIMRunner::getNimResult() const {
	return nimResult;
}

NIMRunner::~NIMRunner() {
	// TODO Auto-generated destructor stub
}

