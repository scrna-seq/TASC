/*
 * ParsedERCC.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: jiacheng
 */

#include <algorithm>
#include <sstream>
#include <fstream>

#include "ParsedERCC.h"

ParsedERCC::ParsedERCC() {
	this->numErccGenes = 0;
	this->numSamples = 0;
}




bool ParsedERCC::addERCCGene(std::vector<double> countsOfThisERCCGene){

	if (countsOfThisERCCGene.empty()){
		return false;
	}
	if (this->numSamples == 0){
		this->erccCounts.push_back(countsOfThisERCCGene);
		this->numSamples = countsOfThisERCCGene.size();
		this->numErccGenes++;
		return true;
	} else if (this->numSamples == countsOfThisERCCGene.size()){
		this->erccCounts.push_back(countsOfThisERCCGene);
		this->numErccGenes++;
		return true;
	} else {
		return false;
	}
}


void ParsedERCC::removeCols(std::vector<size_t> idxCols) {
	if(idxCols.empty()){return;}

	std::vector<double> erccCountsBeforeRemoval;
	std::vector<double> erccCountsAfterRemoval;
	std::vector<std::vector<double> > newERCCGeneCounts;

    std::vector<std::vector<double> >::iterator it;

    for (it = this->erccCounts.begin(); it != this->erccCounts.end(); ++it) {
		erccCountsAfterRemoval.clear();
        for (size_t i = 0; i < it->size(); ++i) {
			if (std::find(idxCols.begin(),idxCols.end(),i)==idxCols.end()){
                erccCountsAfterRemoval.push_back(it->at(i));
			}
		}
		newERCCGeneCounts.push_back(erccCountsAfterRemoval);
	}

//	BOOST_FOREACH(erccCountsBeforeRemoval, this->erccCounts){
//		erccCountsAfterRemoval.clear();
//		for(size_t i=0; i<erccCountsBeforeRemoval.size(); ++i){
//			if (std::find(idxCols.begin(),idxCols.end(),i)==idxCols.end()){
//				erccCountsAfterRemoval.push_back(erccCountsBeforeRemoval.at(i));
//			}
//		}
//		newERCCGeneCounts.push_back(erccCountsAfterRemoval);
//	}

	this->erccCounts = newERCCGeneCounts;
	this->numSamples = erccCountsAfterRemoval.size();

	std::vector<double> newAlpha;
	std::vector<double> newBeta;
	for(size_t i=0; i<std::min(this->alpha.size(), this->beta.size()); ++i){
		if(std::find(idxCols.begin(), idxCols.end(), i)==idxCols.end()){
			newAlpha.push_back(this->alpha.at(i));
			newBeta.push_back(this->beta.at(i));
		}
	}

	this->alpha = newAlpha;
	this->beta = newBeta;
}

void ParsedERCC::parseFromFile(std::string fileNameERCC){
	std::fstream file_in(fileNameERCC);
	if (!file_in.good()){
		timed_log("Cannot access "+fileNameERCC+".\nNow exiting.");
		exit(1);
	}

	std::string line;
	std::string tokenBuffer;
	double tokenValueBuffer;
	std::vector<double> countsValueArr;

	while(std::getline(file_in, line)){

		// retrieve ercc names;
		std::istringstream iss(line);
		iss >> tokenBuffer;
		this->ERCCNames.push_back(tokenBuffer);

		// retrieve true concentration
		iss >> tokenBuffer;

		try{
			tokenValueBuffer = std::stod(tokenBuffer,NULL);
		} catch (const std::invalid_argument& ia){
			timed_log("Values in ERCC file are not numeric! Now exiting.");
			exit(1);
		} catch (const std::out_of_range& oor) {
			timed_log("Values in ERCC file are out of range! Now exiting.");
			exit(1);
		}

		this->erccTrueConcentration.push_back(tokenValueBuffer);


		// retrieve ercc lengths
		iss >> tokenBuffer;

		try{
			tokenValueBuffer = std::stod(tokenBuffer,NULL);
		} catch (const std::invalid_argument& ia){
			timed_log("Values in ERCC file are not numeric! Now exiting.");
			exit(1);
		} catch (const std::out_of_range& oor) {
			timed_log("Values in ERCC file are out of range! Now exiting.");
			exit(1);
		}

		this->erccLengths.push_back(tokenValueBuffer);


		// retrieve counts
		iss >>tokenBuffer;
		std::istringstream iss2(tokenBuffer);
		std::string countTokenBuffer;
		countsValueArr.clear();

		double countValueBuffer;
		while(std::getline(iss2, countTokenBuffer, ',')){
			try{
				countValueBuffer = std::stod(countTokenBuffer,NULL);
			} catch (const std::invalid_argument& ia){
				timed_log("Values in ERCC file are not numeric! Now exiting.");
				exit(1);
			} catch (const std::out_of_range& oor) {
				timed_log("Values in ERCC file are out of range! Now exiting.");
				exit(1);
			}
			countsValueArr.push_back(countValueBuffer);
		}

		if (!this->addERCCGene(countsValueArr)){
			timed_log("Sample sizes are not consistent across genes! Now exiting.");
			exit(1);
		}
	}
}

void ParsedERCC::estimateAlphaBetasByOLS() {
	timed_log("Linear regressions estimating cell-specific alpha and beta.");
	double eps = std::numeric_limits<double>::epsilon();
	double tmp, c0, c1, cov00, cov01, cov11, chisq;


	for(size_t i=0;i<this->numSamples;++i){
		std::vector<double> non_zero_y;
		std::vector<double> non_zero_x;
		size_t num_non_zero = 0;

		for(size_t j=0;j<this->numErccGenes;++j){
			tmp = this->erccCounts.at(j).at(i);
			if(tmp>eps){
				non_zero_y.push_back(log(tmp));
				non_zero_x.push_back(log(this->erccTrueConcentration.at(j)));
				num_non_zero++;
			}
		}

		if(num_non_zero<2){
			timed_log("Sample #"+std::to_string((long long)i)+" contains less than 2 non-zero ERCC gene counts. Marked for removal.");
			this->invalidSampleIdx.push_back(i);
			alpha.push_back(NAN);
			beta.push_back(NAN);
			continue;
		}

		if(gsl_fit_linear(&non_zero_x[0],1,&non_zero_y[0],1,num_non_zero,&c0,&c1,&cov00,&cov01,&cov11,&chisq)){
			timed_log("Error in computing linear coefficients (alpha and beta) for Sample #"+std::to_string((long long)i)+". Marked for removal.");
			this->invalidSampleIdx.push_back(i);
			alpha.push_back(NAN);
			beta.push_back(NAN);
			continue;
		}

		alpha.push_back(c0);
		beta.push_back(c1);
	}

}

void ParsedERCC::initERCC(std::string fileNameERCC) {
	this->parseFromFile(fileNameERCC);
	this->estimateAlphaBetasByOLS();
}

std::vector<size_t> ParsedERCC::getInvalidSampleIdx() {
	return this->invalidSampleIdx;
}

std::vector<std::vector<double> > ParsedERCC::getErccCountsTranspose() {
	return this->erccCountsTranspose;
}

ParsedERCC::~ParsedERCC() {
	// TODO Auto-generated destructor stub
}

std::vector<double> ParsedERCC::getCountsOfERCCGene(std::string erccGeneName) {
	size_t idxErccGene = std::find(this->ERCCNames.begin(), this->ERCCNames.end(), erccGeneName) - this->ERCCNames.begin();
	if (this->erccCounts.size() > idxErccGene){
		return this->erccCounts.at(idxErccGene);
	} else {
		return std::vector<double>();
	}

}

std::vector<double> ParsedERCC::getAlpha() {
	return alpha;
}

std::vector<double> ParsedERCC::getBeta() {
	return beta;
}

std::vector<std::vector<double> > ParsedERCC::getErccCounts() {
	return erccCounts;
}

std::vector<double> ParsedERCC::getErccLengths() {
	return erccLengths;
}

std::vector<std::string> ParsedERCC::getErccNames() {
	return ERCCNames;
}

std::vector<double> ParsedERCC::getErccTrueConcentration() {
	return erccTrueConcentration;
}

size_t ParsedERCC::getNumErccGenes() {
	return numErccGenes;
}

size_t ParsedERCC::getNumSamples() {
	return numSamples;
}

bool ParsedERCC::computeTransposedErccCounts() {

	if(this->numErccGenes==0 || this->numSamples == 0){
		return false;
	} else {
		for(size_t i=0; i<this->numSamples; ++i){
			std::vector<double> y_t;
			y_t.push_back(this->erccCounts.at(0).at(i));
			this->erccCountsTranspose.push_back(y_t);
		}
		for(size_t i=0; i<this->numSamples;++i){
			for(size_t j=1;j<this->numErccGenes;++j){
				this->erccCountsTranspose.at(i).push_back(this->erccCounts.at(j).at(i));
			}
		}
		return true;
	}
}

std::vector<double> ParsedERCC::computeTotalCellCountsERCC() {
    std::vector<double> results(this->getNumSamples());
	for(std::vector<std::vector<double> >::iterator it = this->erccCounts.begin();it!=this->erccCounts.end();++it){
        for(size_t i=0; i<it->size(); ++i){
            results[i] = results[i] + it->at(i);
        }
    }

//    for (auto& gene : this->erccCounts){
//        for(size_t i=0; i<gene.size(); ++i){
//            results[i] = results[i] + gene[i];
//        }
//    }
    this->totalCountsByCellERCC = results;
    return this->totalCountsByCellERCC;
}
