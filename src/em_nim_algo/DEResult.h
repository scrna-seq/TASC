/*
 * NIMResult.h
 *
 *  Created on: Jul 19, 2015
 *      Author: cheng
 */

#ifndef NIMRESULT_H_
#define NIMRESULT_H_
#include <vector>
#include <string>
#include <iostream>

class DEResult {
public:
	DEResult();
	std::string gene_name;
	std::vector<std::string> covariate_names;
	std::vector<size_t> test_cols;
	std::vector<double> theta_t;
	bool aborted;
	std::vector<double> beta_hat;
	std::vector<std::pair<std::string,double> > lrt_stat;
	std::vector<std::pair<std::string,double> > lrt_pval;
	double full_lklhood;
	std::vector<std::pair<std::string,double > > reduced_lklhood;
	std::vector<std::pair<std::string,double > > reduced_sigma;
	double full_sigma;
	bool full_converge;
	std::vector<std::pair<std::string,bool > > reduced_converge;
	void print_all(std::ostream &outputFileHandle, bool verbose);
	virtual ~DEResult();
};

#endif /* NIMRESULT_H_ */
