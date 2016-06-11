/*
 * EMRunner.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: jiacheng
 */

#include <cmath>
#include <memory>

#include "gsl/gsl_integration.h"
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_blas.h"

#include "EMRunner.h"
#include "../misc/basic_math.h"
#include "../misc/utils.h"
#include "likelihood_functions.h"
#include "adaptiveIntegrate.h"

EMRunner::EMRunner() {
	numCovariates = 0;
    designMatrix = NULL;
    maxIter = 1000;
    initialSigma = 0;
    minRelTol = 1e-5;
    epsabs = pow(std::numeric_limits<double>::epsilon(),0.25L);
    sampleSize=0;
}

void EMRunner::setData(std::pair<std::vector<std::string>, gsl_matrix*>  designMatrix,
		std::pair<std::string, std::vector<double> > geneNameAndCounts,
		std::vector<size_t> covariatesToBeTested,
		double initSigma,
		std::vector<double> alpha,
		std::vector<double> beta,
		std::vector<double> kappa,
		std::vector<double> tau,
		double maxIter,
		double minRelTol) {
	//default scaleFactorBase = 0.98
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

	for(size_t i=0; i<this->geneCounts.size(); ++i){
		this->scaleFactors.push_back(1);
	}
	//	for(size_t i=0; i<this->geneCounts.size(); ++i){
	//		double sf = dpois(std::round(this->geneCounts[i]*scaleFactorBase),this->geneCounts[i]);
	//		if (sf == 0){
	//			timed_log("Scale factor no. "+std::to_string((long long)i) + " is zero. Please adjust scale factor base and retry. Now exiting.");
	//			exit(1);
	//		}else{
	//			this->scaleFactors.push_back(sf);
	//		}
	//	}
}

EMResult EMRunner::runModel() {
	EMResult emResult;
	emResult.covariate_names = this->covariateNames;
	emResult.gene_name = this->geneName;
	emResult.test_cols = this->covariatesToBeTested;
	gsl_matrix *x = this->designMatrix;
	//gsl_matrix_fprintf (stderr, x, "%1.3e");

	std::vector<double> y = this->geneCounts;
	//	timed_log(y, "Y");

	emResult.max_iter = this->maxIter;
	emResult.min_tol = this->minRelTol;
	double epsabs = 1e-300;
	double epsrel = pow(std::numeric_limits<double>::epsilon(),0.25L);
	const size_t num_covariates = x->size2;

	if (this->initialSigma <=0 ){
		timed_log("Exception: initial value of sigma_t <= 0!");
		throw 1;
	}

	// full model;
	gsl_vector *theta_t = gsl_vector_alloc(this->sampleSize);

	for (unsigned int i = 0; i<this->sampleSize; i++){
		gsl_vector_set(theta_t,i,log1p(y.at(i)));
	}

	verbose_timed_log(theta_t, "initial theta_t:");
	double sigma_t = this->initialSigma;

	size_t workspace_size = 1000;
	gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc (workspace_size);

	gsl_vector *ewj_gsl_vec = gsl_vector_alloc(this->sampleSize);
	std::vector<double> var_wj(this->sampleSize);

	gsl_vector *prev_beta_hat = gsl_vector_alloc(num_covariates);

	double pars[9];
	void * ptr_pars = pars;
	//used to be phi, now useless. set it to 0. DO NOT REFERENCE in integrand.
	pars[0] = 0;
	gsl_function F_EWSQ;
	gsl_function F_EW;
	gsl_function F_P;
	F_EWSQ.function = &ew_sq_integrand;
	F_EWSQ.params = ptr_pars;
	F_EW.function = &ew_integrand;
	F_EW.params = ptr_pars;
	F_P.function = &single_lklhood;
	F_P.params = ptr_pars;
	double delta_beta_t =  std::numeric_limits<double>::max();

	double chisq;
	gsl_vector *beta_hat = gsl_vector_alloc(num_covariates);

	gsl_matrix *cov = gsl_matrix_alloc(num_covariates,num_covariates);
	gsl_multifit_linear_workspace * wspc = gsl_multifit_linear_alloc(this->sampleSize, num_covariates);

	for (unsigned int step = 0; step < this->maxIter; step++){
		//e-step
		verbose_timed_log("Full EM Gene "+this->geneName+" Step "+std::to_string((long long int)step));
		pars[6] = sigma_t;
				std::vector<double> ew(this->sampleSize);
				std::vector<double> p_vec(this->sampleSize);
				std::vector<double> ew_sq(this->sampleSize);

//		verbose_timed_log("sample size:"+std::to_string(sampleSize));
		for (unsigned int i = 0; i < this->sampleSize; i++){
			pars[1] = this->kappa.at(i);
			pars[2] = this->tau.at(i);
			pars[3] = this->alpha.at(i);
			pars[4] = this->beta.at(i);
			pars[7] = this->scaleFactors[i];
			pars[8] = y.at(i);
			pars[5] = gsl_vector_get(theta_t,i);
			double integrated_result_ew;
			double integrated_result_ew_sq;
			double integrated_result_p;
			double error;

			int status = adaptiveIntegrate(&F_P,pars, epsabs,epsrel, workspace, &integrated_result_p, &error);
			if (status!=0){
				return emResult;
			}
			status = adaptiveIntegrate (&F_EWSQ, pars,epsabs,epsrel, workspace, &integrated_result_ew_sq, &error);
//			verbose_timed_log("Step:"+std::to_string(i));
			if (status){
				std::cout<<gsl_strerror (status);
				return emResult;
			}
			status = adaptiveIntegrate(&F_EW,pars, epsabs,epsrel, workspace, &integrated_result_ew, &error);
			if (status!=0){
				return emResult;
			}
			//			timed_log("Running EM."+std::to_string((long long) step)+"\t"+std::to_string((long long) i));

			gsl_vector_set(ewj_gsl_vec,i,integrated_result_ew/integrated_result_p);
			var_wj[i] = integrated_result_ew_sq/integrated_result_p;
						ew[i] = integrated_result_ew;
						p_vec[i] = integrated_result_p;
						ew_sq[i] = integrated_result_ew_sq;

						std::vector<double> intermediateResults;
//						intermediateResults.push_back(integrated_result_p);
//						intermediateResults.push_back(integrated_result_ew);
//						intermediateResults.push_back(integrated_result_ew_sq);
			//
			//						timed_log("Running EM."+std::to_string((long long) step));
			//						std::vector<double>pars_vec(pars[1], pars[8]);
			//
			//						intermediateResults.insert(intermediateResults.end(), pars_vec.begin(), pars_vec.end());
//						verbose_timed_log(intermediateResults,"intermediate results\t"+std::to_string(i)+":");


		}

		/*lm_res <- lm(formula = y~x, data.frame(y=ew_js[,1], x=x))
		theta_t <- predict(lm_res)
		sigma_t_sq <- mean(ew_js[,2]-2*theta_t*ew_js[,1]+theta_t^2,na.rm = T)
		if (sigma_t_sq<0) {
		print(paste('error! m_step_full: sigma_t_sq = ',sigma_t_sq,sep=''))
		}

		sigma_t <- sqrt(sigma_t_sq)*/

				verbose_timed_log(ew,"EW:");
				verbose_timed_log(ew_sq,"EW_SQ:");
				verbose_timed_log(p_vec,"P:");
		gsl_multifit_linear(x, ewj_gsl_vec, beta_hat, cov, &chisq, wspc);
		emResult.full_pars.push_back(gsl_to_std_vector(beta_hat,num_covariates));

		if(step>1){
			//compare beta_hat to previous beta_hat;
			gsl_vector_sub(prev_beta_hat,beta_hat);
			delta_beta_t = gsl_vector_l2norm(prev_beta_hat,num_covariates);
			emResult.full_path.push_back(delta_beta_t);
		}

		//		timed_log(beta_hat, "beta_hat");
		int status = gsl_blas_dgemv (CblasNoTrans, 1.0, x, beta_hat, 0.0, theta_t);
		if (status != 0) {
			timed_log("Exception: gsl_blas_dgemv error.");
			return emResult;
		}
		//		timed_log(theta_t, "theta_t");
		sigma_t = mean(var_wj);
		sigma_t = sqrt(sigma_t);
		if ((sigma_t <= 0) || (!std::isfinite(sigma_t))) {
			timed_log("Exception: var_wj <= 0.");
			return emResult;
		}
		if (delta_beta_t < this->minRelTol){
			emResult.full_converge=true;
			break;
		}
		gsl_vector_swap(prev_beta_hat,beta_hat);

	}

	emResult.ewj = gsl_to_std_vector(ewj_gsl_vec,this->sampleSize);
	emResult.var_wj = var_wj;
	emResult.full_lklhood = em_marginal_lklhood(kappa,tau,alpha,beta,
			theta_t,sigma_t,this->scaleFactors,y,this->sampleSize,epsabs);
	emResult.full_sigma = sigma_t;
	emResult.beta_hat = gsl_to_std_vector(prev_beta_hat, num_covariates);

	gsl_vector_free(prev_beta_hat);
	gsl_vector_free(beta_hat);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(wspc);


	//reduced model
	prev_beta_hat = gsl_vector_alloc(num_covariates-1);
	beta_hat = gsl_vector_alloc(num_covariates-1);
	cov = gsl_matrix_alloc(num_covariates-1,num_covariates-1);
	wspc = gsl_multifit_linear_alloc(this->sampleSize, num_covariates-1);

	for (size_t test_col_idx=0; test_col_idx<this->covariatesToBeTested.size();++test_col_idx){
		size_t test_col = this->covariatesToBeTested.at(test_col_idx);
		gsl_matrix *reduced_x = remove_column(x,test_col);
		std::vector<std::vector<double> > reduced_pars_for_this_covariate;
		std::vector<double> reduced_path_for_this_covariate;
		sigma_t = this->initialSigma;
		std::string covariate_name_being_tested = emResult.covariate_names.at(test_col-1);

		for (unsigned int i = 0; i<this->sampleSize; i++){
			gsl_vector_set(theta_t,i,log1p(y.at(i)));
		}
		delta_beta_t =  std::numeric_limits<double>::max();

		for (unsigned int step = 0; step < this->maxIter; step++){
			//e-step
			verbose_timed_log("Reduced EM Gene "+this->geneName+" Step "+std::to_string((long long int)step));
			pars[6] = sigma_t;
			for (unsigned int i = 0; i < this->sampleSize; i++){
				pars[1] = this->kappa.at(i);
				pars[2] = this->tau.at(i);
				pars[3] = this->alpha.at(i);
				pars[4] = this->beta.at(i);

				pars[7] = this->scaleFactors[i];
				pars[8] = y.at(i);
				pars[5] = gsl_vector_get(theta_t,i);
				double integrated_result_ew;
				double integrated_result_ew_sq;
				double integrated_result_p;
				double error;
				int status = gsl_integration_cquad(&F_EWSQ, pars[5]-pars[6]*40,pars[5]+pars[6]*40,epsabs,epsrel, workspace, &integrated_result_ew_sq, &error,NULL);
				if (status!=0){
					return emResult;
				}
				status = gsl_integration_cquad(&F_EW,pars[5]-pars[6]*40,pars[5]+pars[6]*40, epsabs,epsrel, workspace, &integrated_result_ew, &error, NULL);
				if (status!=0){
					return emResult;
				}
				gsl_integration_cquad(&F_P, pars[5]-pars[6]*40,pars[5]+pars[6]*40, epsabs,epsrel, workspace, &integrated_result_p, &error, NULL);
				if (status!=0){
					return emResult;
				}
				gsl_vector_set(ewj_gsl_vec,i,integrated_result_ew/integrated_result_p);
				var_wj[i] = integrated_result_ew_sq/integrated_result_p;
			}


			/*lm_res <- lm(formula = y~x, data.frame(y=ew_js[,1], x=x))
	theta_t <- predict(lm_res)
	sigma_t_sq <- mean(ew_js[,2]-2*theta_t*ew_js[,1]+theta_t^2,na.rm = T)
	if (sigma_t_sq<0) {
	print(paste('error! m_step_full: sigma_t_sq = ',sigma_t_sq,sep=''))
	}
	sigma_t <- sqrt(sigma_t_sq)*/
			gsl_multifit_linear(reduced_x, ewj_gsl_vec, beta_hat, cov, &chisq, wspc);
			reduced_pars_for_this_covariate.push_back(gsl_to_std_vector(beta_hat,num_covariates-1));


			verbose_timed_log(theta_t,"reduced theta_t:");
			if(step>1){
				//compare beta_hat to previous beta_hat;
				gsl_vector_sub(prev_beta_hat,beta_hat);
				delta_beta_t = gsl_vector_l2norm(prev_beta_hat,num_covariates-1);
				reduced_path_for_this_covariate.push_back(delta_beta_t);
			}


			int status = gsl_blas_dgemv (CblasNoTrans, 1.0, reduced_x, beta_hat, 0.0, theta_t);
			if (status != 0) {
				timed_log("Exception: gsl_blas_dgemv error.");
				return emResult;
			}

			sigma_t = mean(var_wj);
			if ((sigma_t <= 0) || (!std::isfinite(sigma_t)) ) {
				timed_log("Exception: var_wj <= 0.");
				return emResult;
			}
			sigma_t = sqrt(sigma_t);
			if (delta_beta_t < this->minRelTol){
				emResult.reduced_converge.push_back(std::pair<std::string, bool> (covariate_name_being_tested, true));
				break;
			}
			gsl_vector_swap(prev_beta_hat,beta_hat);
		}
		emResult.reduced_converge.push_back(std::pair<std::string, bool> (covariate_name_being_tested, false));
		double log_lklhood_reduced_for_this_covariate = em_marginal_lklhood(kappa,tau,alpha,beta,theta_t,sigma_t,this->scaleFactors,y,this->sampleSize,epsabs);
		emResult.reduced_lklhood.push_back(std::pair<std::string, double> (covariate_name_being_tested,log_lklhood_reduced_for_this_covariate));
		emResult.reduced_sigma.push_back(std::pair<std::string, double> (covariate_name_being_tested,sigma_t));
		double lrt_stat_for_this_covariate = 2*(emResult.full_lklhood-log_lklhood_reduced_for_this_covariate);
		if (lrt_stat_for_this_covariate<0) lrt_stat_for_this_covariate=0;
		double lrt_pval_for_this_covariate = 1-pchisq(lrt_stat_for_this_covariate);
		emResult.lrt_pval.push_back(std::pair<std::string, double> (covariate_name_being_tested,lrt_pval_for_this_covariate));
		emResult.lrt_stat.push_back(std::pair<std::string, double> (covariate_name_being_tested,lrt_stat_for_this_covariate));
		emResult.reduced_pars.push_back(std::pair<std::string, std::vector<std::vector<double> > > (covariate_name_being_tested, reduced_pars_for_this_covariate));
		emResult.reduced_path.push_back(std::pair<std::string, std::vector<double> > (covariate_name_being_tested,reduced_path_for_this_covariate));
		gsl_matrix_free(reduced_x);
	}

	gsl_vector_free(theta_t);
	gsl_vector_free(ewj_gsl_vec);
	gsl_vector_free(prev_beta_hat);
	gsl_vector_free(beta_hat);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(wspc);
	gsl_integration_cquad_workspace_free(workspace);
	emResult.aborted=false;
	return emResult;
}

const std::string EMRunner::getGeneName() const {
	return geneName;
}

EMRunner::~EMRunner() {
	// TODO Auto-generated destructor stub
}


