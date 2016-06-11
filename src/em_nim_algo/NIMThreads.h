/*
 * NIMThreads.h
 *
 *  Created on: Nov 3, 2015
 *      Author: jiacheng
 */

#ifndef EM_NIM_ALGO_NIMTHREADS_H_
#define EM_NIM_ALGO_NIMTHREADS_H_
#include "DEResult.h"
#include "NIMRunner.h"
#include "../misc/CLIArgs.h"
#include "../misc/ParsedData.h"

class NIMThreads {
private:
	std::vector<DEResult> nimRes;
	std::vector<NIMRunner> nimRunners;
	size_t numGenes;
	size_t numThreads;
	std::string outputFileName;
	std::vector<std::string> nameCovariates;
	std::vector<size_t> testedCovariatesIdx;

	void printHeader(std::ostream &outputFileHandle);
public:
	NIMThreads(ParsedData &data, CLIArgs &prog_params);
	void runThreads();
	void printAllResults();
	virtual ~NIMThreads();
};

#endif /* EM_NIM_ALGO_NIMTHREADS_H_ */
