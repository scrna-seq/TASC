/*
 * ParsedData.h
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */

#ifndef MISC_PARSEDDATA_H_
#define MISC_PARSEDDATA_H_
#include <string>
#include <vector>
#include <map>
#include <memory>

#include<gsl/gsl_matrix.h>

#include "ParsedCountsY.h"
#include "ParsedERCC.h"
#include "ParsedCovariatesX.h"
#include "CLIArgs.h"
#include "basic_math.h"
#include "../param_estimation/InitMargLklhoodArgs.h"

class ParsedData {
private:
	std::vector<double> kappa;
	std::vector<double> tau;
	std::vector<double> alpha;
	std::vector<double> beta;
	std::vector<double> cellTotalCounts;
	std::vector<double> cellTotalCountsERCC;
	std::vector<double> sizeFactors;
	ParsedERCC erccParser;
	ParsedCountsY yParser;
	ParsedCovariatesX xParser;
	mvnorm_params initialHyperParamValues;
	mvnorm_params estimatedHyperParamValues;
	InitMargLklhoodArgs marginalLikelihoodArguments;
	void synchronizeAllParsedData(std::vector<size_t> invalidSampleIdx);
	void parseY(std::string fileNameY);
	void parseX(std::string fileNameX);
	void parseERCC(std::string fileNameERCC);
	void computeKappaTau(CLIArgs runningParams);
	void computeEmpiricalBayesMean(CLIArgs runningParams);
	bool isEqualSampleSize();
	size_t sampleSize;
	void estimateKappaTaus();
	void parseABKT(std::string abktFileName);
	bool disableAdjustmentForSize;

public:
	ParsedData();
	void parseAndEstimateHyperParams(CLIArgs runningParams);
	std::vector<double> getAlpha();
	std::vector<double> getBeta();
	std::vector<double> getKappa();
	std::vector<double> getTau();
	std::vector<std::string> getGeneNames();
	std::vector<double> getCountsOfGene(std::string geneName);
	gsl_matrix* getDesignMatrix();
	std::vector<std::string> getCovNames();
	std::vector<size_t> getTestColIdx();
	size_t getNumGenes();
	std::pair<std::vector<std::string>, gsl_matrix*> getCovNamesAndDesignMat();
	std::pair<std::string, std::vector<double> > getGeneNameAndCounts(std::string geneName);
	std::vector<size_t> getCovariatesToBeTested();
	std::vector<double> getCellTotalCounts();
	std::vector<double> getCellTotalCountsERCC();
	std::vector<double> getSizeFactors();
	virtual ~ParsedData();
};

#endif /* MISC_PARSEDDATA_H_ */

