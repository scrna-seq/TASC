/*
 * ParsedCovariatesX.h
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */

#ifndef MISC_PARSEDCOVARIATESX_H_
#define MISC_PARSEDCOVARIATESX_H_
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <gsl/gsl_matrix.h>

class ParsedCovariatesX {
private:
	//need this vector to record the order of the covariates

	std::map<std::string, std::vector<double> > parsedX;
	gsl_matrix *designMatrix;
	size_t numCovariates;
	size_t numSamples;

public:
	std::vector<size_t> idxCovToBeTested;
	std::vector<std::string> covariateNames;
	ParsedCovariatesX();
	void initX(std::string fileNameX, std::string covIdxToBeTested);
	gsl_matrix* getDesignMatrix() const;
	std::pair<std::vector<std::string>, gsl_matrix *> getCovNamesAndDesignMatrix();
	size_t getNumCovariates() const;
	size_t getNumSamples() const;
	bool insertSample(std::string covariateName, std::vector<double> covariateValues);
	void removeCols(std::vector<size_t> idxCols);
	std::vector<double> getCovariateValues(std::string covariateName);
	virtual ~ParsedCovariatesX();
	void addSizeFactorCol(std::vector<double> sf);
};

#endif /* MISC_PARSEDCOVARIATESX_H_ */
