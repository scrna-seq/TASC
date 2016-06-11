/*
 * IntegrationArgs.h
 *
 *  Created on: Aug 18, 2015
 *      Author: cheng
 */

#ifndef INTEGRATIONARGS_H_
#define INTEGRATIONARGS_H_
#include <vector>
class IntegrationArgs {
public:
	IntegrationArgs();
	double alpha;
	double beta;
	std::vector<double> y;
	std::vector<double> true_concentrations;
	std::vector<double> ercc_lens;
	double mean_kappa;
	double mean_tau;
	double sd_kappa;
	double sd_tau;
	double rho_kappa_tau;
	void displayAll();
	void SetAllArgs(double alpha,
			double beta,
			std::vector<double> y,
			std::vector<double> true_concentrations,
			std::vector<double> ercc_lens,
			double mean_kappa,
			double mean_tau,
			double sd_kappa,
			double sd_tau,
			double rho_kappa_tau);

	virtual ~IntegrationArgs();
};

#endif /* INTEGRATIONARGS_H_ */
