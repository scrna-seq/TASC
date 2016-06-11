/*
 * IntegrationArgs.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: cheng
 */

#include "IntegrationArgs.h"
#include <iostream>

/*
 * Definition of a class passing the args
 * for the integration
 * */

IntegrationArgs::IntegrationArgs() {
	this->alpha=0;
	this->beta=0;
	this->mean_kappa = 0;
	this->mean_tau = 0;
	this->sd_kappa = 0;
	this->sd_tau = 0;
	this->rho_kappa_tau = 0;
	// TODO Auto-generated constructor stub

}

IntegrationArgs::~IntegrationArgs() {
	// TODO Auto-generated destructor stub
}

void IntegrationArgs::SetAllArgs(
		double alpha,
		double beta,
		std::vector<double> y,
		std::vector<double> true_concentrations,
		std::vector<double> ercc_lens,
		double mean_kappa,
		double mean_tau,
		double sd_kappa,
		double sd_tau,
		double rho_kappa_tau){
	this->alpha = alpha;
	this->beta = beta;
	this->y = y;
	this->true_concentrations = true_concentrations;
	this->ercc_lens = ercc_lens;
	this->mean_kappa = mean_kappa;
	this->mean_tau = mean_tau;
	this->sd_kappa = sd_kappa;
	this->sd_tau = sd_tau;
	this->rho_kappa_tau = rho_kappa_tau;
}

void IntegrationArgs::displayAll() {
	std::cerr<<this->alpha<<std::endl;
	std::cerr<<this->beta<<std::endl;
	std::cerr<<this->y.size()<<std::endl;
	std::cerr<<this->true_concentrations.size()<<std::endl;
	std::cerr<<this->ercc_lens.size()<<std::endl;
	std::cerr<<this->mean_kappa<<std::endl;
	std::cerr<<this->mean_tau<<std::endl;
	std::cerr<<this->sd_kappa<<std::endl;
	std::cerr<<this->sd_tau<<std::endl;
	std::cerr<<this->rho_kappa_tau<<std::endl;
}
