/*
 * EMThreads.h
 *
 *  Created on: Nov 3, 2015
 *      Author: jiacheng
 */

#ifndef EM_NIM_ALGO_EMTHREADS_H_
#define EM_NIM_ALGO_EMTHREADS_H_
#include "EMResult.h"
#include "EMRunner.h"
#include <string>
#include <memory>
#include "../misc/CLIArgs.h"
#include "../misc/ParsedData.h"

class EMThreads {
private:
	std::vector<EMResult> emRes;
	std::vector<EMRunner> emRunners;
	size_t numGenes;
	size_t numThreads;
	std::string outputFileName;
public:
	EMThreads(ParsedData &data, CLIArgs &prog_params);
	void runThreads();
	void printAllResults();
	virtual ~EMThreads();
};



#endif /* EM_NIM_ALGO_EMTHREADS_H_ */
