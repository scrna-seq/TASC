/*
 * EMThreads.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: jiacheng
 */

#include "EMThreads.h"
#include "omp.h"
#include <string>
#include <vector>
#include <memory>
#include <fstream>
EMThreads::EMThreads(ParsedData &data, CLIArgs &prog_params) {
	numGenes = 0;
	numThreads = 1;
	outputFileName = "";
	this->numGenes = data.getNumGenes();
	this->outputFileName = prog_params.getOutputFileName()+".em";
	this->numThreads = prog_params.getNumThreads();
	this->emRes = std::vector<EMResult>(data.getNumGenes());
//	std::shared_ptr<std::pair<std::string, std::vector<double> > >
//  ptr(new std::pair<std::string, std::vector<double> >);
	this->emRunners = std::vector<EMRunner>(data.getNumGenes());
	std::vector<std::string> geneNames = data.getGeneNames();
	for(size_t i=0;i<data.getNumGenes();++i){
		EMRunner ptr;
		this->emRunners[i] = ptr;
		this->emRunners[i].setData(data.getCovNamesAndDesignMat(),
				data.getGeneNameAndCounts(geneNames[i]),
				data.getCovariatesToBeTested(),
				prog_params.getInitSigma(),
				data.getAlpha(),data.getBeta(),data.getKappa(),data.getTau(),
				prog_params.getMaxIterEm(),prog_params.getMinRelTolEm());
	}
}

void EMThreads::runThreads() {
	omp_set_num_threads(this->numThreads);
	timed_log("Starting EM Algorithm with "+std::to_string((long long)numThreads)+" threads.");
#pragma omp parallel for
	for (size_t i=0; i<this->numGenes;i++){
		this->emRes[i] = this->emRunners[i].runModel();
	}
}

void EMThreads::printAllResults() {
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
		std::cout <<"EM results:"<< std::endl;
	}

	for(size_t i=0; i<this->numGenes;++i){
		this->emRes[i].print_all(*outputDevice, prog_params.isVerbose());
	}

	try{
		file_out.close();
	} catch (std::ofstream::failure &err){

	}
}

EMThreads::~EMThreads() {

	// TODO Auto-generated destructor stub
}

