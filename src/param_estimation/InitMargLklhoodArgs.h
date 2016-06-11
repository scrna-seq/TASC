/*
 * InitMargLklhoodArgs.h
 *
 *  Created on: Sep 2, 2015
 *      Author: cheng
 */

#ifndef INITMARGLKLHOODARGS_H_
#define INITMARGLKLHOODARGS_H_
#include <vector>
#include <cstddef>
class InitMargLklhoodArgs {
public:
	std::vector<std::vector<double> > y;
	std::vector<double> true_concentrations;
	std::vector<double> alphas;
	std::vector<double> betas;
	std::vector<double> ercc_lens;
	size_t cubatureMaxIter;
	double cubatureMinRelTol;

	InitMargLklhoodArgs();
	void setData(std::vector<std::vector<double> > y,
	std::vector<double> true_concentrations,
	std::vector<double> alphas,
	std::vector<double> betas,
	std::vector<double> ercc_lens,
	size_t cubatureMaxIter,
	double cubatureMinRelTol
	);
	void displayAll();

	virtual ~InitMargLklhoodArgs();
};

#endif /* INITMARGLKLHOODARGS_H_ */
