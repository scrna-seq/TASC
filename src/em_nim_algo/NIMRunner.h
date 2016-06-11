/*
 * NIMRunner.h
 *
 *  Created on: Oct 30, 2015
 *      Author: jiacheng
 */

#ifndef EM_NIM_ALGO_NIMRUNNER_H_
#define EM_NIM_ALGO_NIMRUNNER_H_

#include <memory>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <gsl/gsl_matrix.h>

#include "DEResult.h"



class NIMRunner {
private:
	DEResult nimResult;
	size_t sampleSize;
	std::vector<double> alpha;
	std::vector<double> beta;
	std::vector<double> kappa;
	std::vector<double> tau;
	std::vector<std::string> covariateNames;
	size_t numCovariates;
	gsl_matrix *designMatrix;
	std::string geneName;
	std::vector<double> geneCounts;
	std::vector<size_t> covariatesToBeTested;
	std::vector<double> scaleFactors;
	size_t maxIter;
	double initialSigma;
	double minRelTol;
	double epsabs;
	bool laplace;
	size_t maxIterCQUAD;
public:
	NIMRunner();
	void runModel();
	void setData(std::pair<std::vector<std::string>, gsl_matrix*> designMatrix,
				 std::pair<std::string, std::vector<double> > geneNameAndCounts,
				 std::vector<size_t> covariatesToBeTested, double initSigma,
				 std::vector<double> alpha,
				 std::vector<double> beta,
				 std::vector<double> kappa,
				 std::vector<double> tau,
				 double maxIter,
				 double minRelTol,
				 bool laplace,
				 size_t max_iter_cquad);
	virtual ~NIMRunner();
	const DEResult & getNimResult() const;
};

#endif /* EM_NIM_ALGO_NIMRUNNER_H_ */
