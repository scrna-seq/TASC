/*
 * EMResult.cpp
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#include "EMResult.h"
#include <cmath>

EMResult::EMResult() {
	full_converge = false;
	full_sigma = nan("");
	max_iter = 0;
	min_tol = nan("");
	aborted = true;
	full_lklhood = nan("");
	// TODO Auto-generated constructor stub
}

EMResult::~EMResult() {
	// TODO Auto-generated destructor stub
}

void EMResult::print_all(std::ostream &outputFileHandle, bool verbose) {
		outputFileHandle << "Gene Name:\t" << gene_name << std::endl;
	outputFileHandle << "EM Aborted:\t"  << aborted << std::endl;
	if (!beta_hat.empty()){
		outputFileHandle << "Full Beta Hat:\t";
		for (std::vector<double>::const_iterator i = beta_hat.begin(); i != beta_hat.end(); ++i){
			outputFileHandle << *i << "\t";
		}
		outputFileHandle << std::endl;
	} else {
		outputFileHandle << "Full Beta Hat:\tEmpty" << std::endl;
	}
	outputFileHandle << "Tested Covariate\tLRT Stat\tP-Val\tReduced Likelihood\tReduced Sigma\tReduced Convergence"<<std::endl;

	for(size_t i=0; i<reduced_lklhood.size(); ++i){
		outputFileHandle << lrt_stat.at(i).first << "\t" << lrt_stat.at(i).second << "\t" << lrt_pval.at(i).second << "\t" << reduced_lklhood.at(i).second << "\t" << reduced_sigma.at(i).second << "\t" << reduced_converge.at(i).second << std::endl;
	}

	outputFileHandle << "Full_Sigma:\t"<<full_sigma<<std::endl;
	outputFileHandle << "Full_Converge:\t"<<full_converge<<std::endl;
	outputFileHandle << "Full log-likelihood:\t"<<full_lklhood<<std::endl;

	if (verbose) {
		if (!full_pars.empty()){
			outputFileHandle << "Full Parameter Path:"<<std::endl;
			for (std::vector<std::vector<double> >::const_iterator i = full_pars.begin(); i != full_pars.end(); ++i){
				for(std::vector<double>::const_iterator j = i->begin(); j != i->end(); ++j){
					outputFileHandle << *j << "\t";
				}
				outputFileHandle << std::endl;
			}
		} else {
			outputFileHandle <<"Full Parameter Path:\tEmpty"<<std::endl;
		}
		if (!ewj.empty()){
			outputFileHandle << "EW_j:\t";
			for (std::vector<double>::const_iterator i = ewj.begin(); i != ewj.end(); ++i){
				outputFileHandle << *i << "\t";
			}
			outputFileHandle << std::endl;
		} else {
			outputFileHandle << "EW_j:\tEmpty"<<std::endl;
		}

		if (!var_wj.empty()){
			outputFileHandle << "VarW_j:\t";
			for (std::vector<double>::const_iterator i = var_wj.begin(); i != var_wj.end(); ++i){
				outputFileHandle << *i << "\t";
			}
			outputFileHandle << std::endl;
		} else {
			outputFileHandle << "VarW_j:\tEmpty"<<std::endl;
		}

		if (!theta_t.empty()){
			outputFileHandle << "Theta_t:\t";
			for (std::vector<double>::const_iterator i = theta_t.begin(); i != theta_t.end(); ++i){
				outputFileHandle << *i << "\t";
			}
			outputFileHandle << std::endl;
		} else {
			outputFileHandle << "Theta_t:\tEmpty"<<std::endl;
		}
		if (!full_path.empty()){
			outputFileHandle << "Full Path:\t";
			for (std::vector<double>::const_iterator i = full_path.begin(); i != full_path.end(); ++i){
				outputFileHandle << *i << "\t";
			}
			outputFileHandle << std::endl;
		} else {
			outputFileHandle << "Full Path:\tEmpty"<<std::endl;
		}
	}

//	outputFileHandle << std::endl;
}
