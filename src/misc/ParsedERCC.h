/*
 * ParsedERCC.h
 *
 *  Created on: Oct 28, 2015
 *      Author: jiacheng
 */

#ifndef MISC_PARSEDERCC_H_
#define MISC_PARSEDERCC_H_

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_matrix.h>

#include "basic_math.h"
#include "utils.h"
class ParsedERCC {
private:
	std::vector<std::string> ERCCNames;
	std::vector<double> alpha;
	std::vector<double> beta;
	std::vector<double> erccTrueConcentration;
	std::vector<double> erccLengths;
	std::vector<std::vector<double> > erccCounts;
	size_t numSamples;
	size_t numErccGenes;
	std::vector<size_t> invalidSampleIdx;
	void estimateAlphaBetasByOLS();
	bool addERCCGene(std::vector<double> countsOfThisERCCGene);
	std::string fileName;
	std::vector<std::vector<double> > erccCountsTranspose;
    std::vector<double> totalCountsByCellERCC;
public:
	ParsedERCC();
	void parseFromFile(std::string fileNameERCC);
	void initERCC(std::string fileNameERCC);
	void removeCols(std::vector<size_t> idxCols);
	std::vector<double> getCountsOfERCCGene(std::string erccGeneName);
	std::vector<size_t> getInvalidSampleIdx();
	std::vector<double> getAlpha();
	std::vector<double> getBeta();
	std::vector<std::vector<double> > getErccCounts();
	std::vector<double> getErccLengths();
	std::vector<std::string> getErccNames();
	std::vector<double> getErccTrueConcentration();
	bool computeTransposedErccCounts();
	size_t getNumErccGenes();
	size_t getNumSamples();
	void logSelf();
	virtual ~ParsedERCC();
	std::vector<double> computeTotalCellCountsERCC();
	std::vector<std::vector<double> > getErccCountsTranspose();
};

#endif /* MISC_PARSEDERCC_H_ */
