/*
 * NIMThreads.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: jiacheng
 */

#include "NIMThreads.h"
#include <omp.h>
#include <string>

NIMThreads::NIMThreads(ParsedData &data, CLIArgs &prog_params) {
	numGenes = 0;
	numThreads = 1;
	outputFileName = "";
	this->numGenes = data.getNumGenes();
	this->nameCovariates = data.getCovNamesAndDesignMat().first;
	this->testedCovariatesIdx = data.getCovariatesToBeTested();
	this->outputFileName = prog_params.getOutputFileName() + ".de";
	this->numThreads = prog_params.getNumThreads();
	this->nimRes = std::vector<DEResult>(data.getNumGenes());
	this->nimRunners = std::vector<NIMRunner>(data.getNumGenes());
	std::vector<std::string> geneNames = data.getGeneNames();
	for(size_t i=0;i<data.getNumGenes();++i){
		this->nimRunners[i] = NIMRunner();

		this->nimRunners[i].setData(data.getCovNamesAndDesignMat(), data.getGeneNameAndCounts(geneNames[i]),
									data.getCovariatesToBeTested(),
									prog_params.getInitSigma(), data.getAlpha(), data.getBeta(), data.getKappa(),
									data.getTau(),
									prog_params.getMaxIterNim(), prog_params.getMinRelTolNim(),
									prog_params.isLaplace(),prog_params.getMaxIterCQuad());
	}
}

void NIMThreads::runThreads() {
	omp_set_num_threads(this->numThreads);
	timed_log("Starting NIM Algorithm with "+std::to_string((long long)numThreads)+" threads.");
#pragma omp parallel for schedule(dynamic,1)
	for (size_t i=0; i<this->numGenes;i++){
		this->nimRunners[i].runModel();
		this->nimRes[i] = this->nimRunners[i].getNimResult();
	}
}

void NIMThreads::printHeader(std::ostream &outputFileHandle) {


	outputFileHandle << "Gene Name" << "\t";
	outputFileHandle << "Algorithm Aborted" << "\t";
	for (size_t i = 0; i < this->nameCovariates.size(); ++i) {
		outputFileHandle << this->nameCovariates[i] << "\t";
	}

	for (size_t i = 0; i < this->testedCovariatesIdx.size(); ++i) {
		outputFileHandle << this->nameCovariates[this->testedCovariatesIdx[i] - 1] + " LRT Stat" << "\t";
		outputFileHandle << this->nameCovariates[this->testedCovariatesIdx[i] - 1] + " LRT PVal" << "\t";
		outputFileHandle << this->nameCovariates[this->testedCovariatesIdx[i] - 1] + " Reduced Model Converged" << "\t";
		outputFileHandle << this->nameCovariates[this->testedCovariatesIdx[i] - 1] + " Reduced Model Sigma" << "\t";
		outputFileHandle << this->nameCovariates[this->testedCovariatesIdx[i] - 1] + " Reduced Model Log Likelihood" <<
		"\t";
	}

	outputFileHandle << "Full Model Converged" << "\t";
	outputFileHandle << "Full Model Sigma" << "\t";
	outputFileHandle << "Full Model Log Likelihood" << std::endl;

}

void NIMThreads::printAllResults() {
	std::ofstream file_out;
	std::ostream* outputDevice;
	bool printFileName = false;

	try{
		file_out.open(this->outputFileName);
		outputDevice = &file_out;
	} catch (std::ofstream::failure &writeErr){
		timed_log("Cannot write to "+this->outputFileName+". Now redirecting to stdout.");
		outputDevice = &std::cout;
	}
	if(printFileName){
		std::cout <<"NIM results:"<< std::endl;
	}


	this->printHeader(*outputDevice);
	for(size_t i=0; i<this->numGenes;++i){
		this->nimRes[i].print_all(*outputDevice, prog_params.isVerbose());
	}

	try{
		file_out.close();
	} catch (std::ofstream::failure &err){

	}
}

NIMThreads::~NIMThreads() {
	// TODO Auto-generated destructor stub
}

