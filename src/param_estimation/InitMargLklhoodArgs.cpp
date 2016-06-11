/*
 * InitMargLklhoodArgs.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: cheng
 */

#include "InitMargLklhoodArgs.h"
#include <iostream>

InitMargLklhoodArgs::InitMargLklhoodArgs(){
	cubatureMaxIter = 100000;
	cubatureMinRelTol = 1e-7;
}

void InitMargLklhoodArgs::setData(std::vector<std::vector<double> > y,
		std::vector<double> true_concentrations,
		std::vector<double> alphas,
		std::vector<double> betas,
		std::vector<double> ercc_lens,
		size_t cubatureMaxIter,
		double cubatureMinRelTol) {
	this->y = y;
	this->true_concentrations  = true_concentrations;
	this->alphas = alphas;
	this->betas = betas;
	this->ercc_lens = ercc_lens;
	this->cubatureMaxIter = cubatureMaxIter;
	this->cubatureMinRelTol = cubatureMinRelTol;
}



InitMargLklhoodArgs::~InitMargLklhoodArgs() {
	// TODO Auto-generated destructor stub
}

void InitMargLklhoodArgs::displayAll() {
	std::cerr << this->y.size()<<std::endl;
	std::cerr << this->y[0].size()<<std::endl;
	std::cerr << this->true_concentrations.size()<<std::endl;
	std::cerr << this->alphas.size()<<std::endl;
	std::cerr << this->betas.size()<<std::endl;
	std::cerr << this->ercc_lens.size()<<std::endl;
	std::cerr << this->cubatureMaxIter<<std::endl;
	std::cerr << cubatureMinRelTol<<std::endl;
}
