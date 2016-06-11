/*
 * EMResult.h
 *
 *  Created on: Jul 15, 2015
 *      Author: jiacheng
 */

#ifndef EMRESULT_H_
#define EMRESULT_H_

#include <vector>
#include <string>
#include <iostream>

class EMResult {
public:
	EMResult();
	std::string gene_name;
	std::vector<std::string> covariate_names;
	std::vector<size_t> test_cols;
	bool aborted;
	std::vector<double> beta_hat;
	std::vector<std::pair<std::string,double> > lrt_stat;
	std::vector<std::pair<std::string,double> > lrt_pval;
	std::vector<double> ewj;
	std::vector<double> var_wj;
	std::vector<double> theta_t;
	double full_lklhood;
	std::vector<std::pair<std::string,double > > reduced_lklhood;
	std::vector<std::pair<std::string,double > > reduced_sigma;
	double full_sigma;
	std::vector<std::pair<std::string,std::vector<double> > > reduced_path;
	std::vector<double> full_path;
	bool full_converge;
	std::vector<std::pair<std::string,bool > > reduced_converge;
	unsigned int max_iter;
	double min_tol;
	std::vector<std::pair<std::string,std::vector<std::vector<double> > > > reduced_pars;
	std::vector<std::vector<double> > full_pars;
	void print_all(std::ostream &outputFileHandle, bool verbose);
	virtual ~EMResult();
};

#endif /* EMRESULT_H_ */
