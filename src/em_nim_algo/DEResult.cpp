/*
 * NIMResult.cpp
 *
 *  Created on: Jul 19, 2015
 *      Author: cheng
 */

#include "DEResult.h"
#include <cmath>

DEResult::DEResult() {
	// TODO Auto-generated constructor stub
	aborted=true;
	full_lklhood=nan("");
	full_sigma=nan("");
	full_converge=false;
}

/*std::string gene_name;
std::vector<std::string> covariate_names;
std::vector<size_t> test_cols;
bool aborted;
std::vector<double> beta_hat;
std::vector<std::pair<std::string,double> > lrt_stat;
std::vector<std::pair<std::string,double> > lrt_pval;
double full_lklhood;
std::vector<std::pair<std::string,double > > reduced_lklhood;
std::vector<std::pair<std::string,double > > reduced_sigma;
double full_sigma;
bool full_converge;
std::vector<std::pair<std::string,bool > > reduced_converge;*/
void DEResult::print_all(std::ostream &outputFileHandle, bool verbose){

	outputFileHandle << gene_name << "\t";
	outputFileHandle << aborted << "\t";
	if (beta_hat.empty()) {
		for (size_t i = 0; i < this->covariate_names.size(); ++i) {
			outputFileHandle << "NA\t";
		}
	} else {
		for (size_t i = 0; i < beta_hat.size(); ++i) {
			outputFileHandle << beta_hat[i] << "\t";
		}
	}
	for(size_t i=0; i<reduced_lklhood.size(); ++i){
		outputFileHandle << lrt_stat.at(i).second << "\t" << lrt_pval.at(i).second << "\t" <<
		reduced_converge.at(i).second << "\t" << std::exp(reduced_sigma.at(i).second) << "\t" <<
		reduced_lklhood.at(i).second << "\t";
	}

	outputFileHandle << full_converge << "\t";
	outputFileHandle << std::exp(full_sigma) << "\t";
	outputFileHandle << full_lklhood << std::endl;

//} else {
//	outputFileHandle << "Gene Name:\t" << gene_name << std::endl;
//	outputFileHandle << "NIM Aborted:\t"  << aborted << std::endl;
//	if (!beta_hat.empty()){
//		outputFileHandle << "Full Beta Hat:\t";
//		for (std::vector<double>::const_iterator i = beta_hat.begin(); i != beta_hat.end(); ++i){
//			outputFileHandle << *i << "\t";
//		}
//		outputFileHandle << std::endl;
//	} else {
//		outputFileHandle << "Full Beta Hat:\tEmpty" << std::endl;
//	}
//	outputFileHandle << "Tested Covariate\tLRT Stat\tP-Val\tReduced Likelihood\tReduced Sigma\tReduced Convergence"<<std::endl;
//
//	for(size_t i=0; i<reduced_lklhood.size(); ++i){
//		outputFileHandle << this->lrt_stat.at(i).first << "\t" << lrt_stat.at(i).second << "\t" << lrt_pval.at(i).second << "\t" << reduced_lklhood.at(i).second << "\t" << reduced_sigma.at(i).second << "\t" << reduced_converge.at(i).second << std::endl;
//	}
//
//	outputFileHandle << "Full_Sigma:\t"<<full_sigma<<std::endl;
//	outputFileHandle << "Full_Converge:\t"<<full_converge<<std::endl;
//	outputFileHandle << "Full log-likelihood:\t"<<full_lklhood<<std::endl;
//}

}


DEResult::~DEResult() {
	// TODO Auto-generated destructor stub
}
