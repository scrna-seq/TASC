/*
 * EMRunner.h
 *
 *  Created on: Oct 29, 2015
 *      Author: jiacheng
 */

#ifndef EM_NIM_ALGO_EMRUNNER_H_
#define EM_NIM_ALGO_EMRUNNER_H_

#include <vector>
#include <memory>
#include <string>
#include <limits>
#include <cmath>
#include <gsl/gsl_matrix.h>

#include "../misc/CLIArgs.h"
#include "EMResult.h"
class EMRunner {
private:
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

public:
	EMRunner();
	EMResult runModel();
	void setData(std::pair<std::vector<std::string>, gsl_matrix*>  designMatrix,
			std::pair<std::string, std::vector<double> > geneNameAndCounts,
			std::vector<size_t> covariatesToBeTested,
			double initSigma,
			std::vector<double> alpha,
			std::vector<double> beta,
			std::vector<double> kappa,
			std::vector<double> tau,
			double maxIter,
			double minRelTol
			);
	virtual ~EMRunner();
	const std::string getGeneName() const;
};

#endif /* EM_NIM_ALGO_EMRUNNER_H_ */
