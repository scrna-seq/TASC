/*
 * ParsedCountsY.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */


#include <fstream>
#include <algorithm>
#include <sstream>

#include "ParsedCountsY.h"
#include "utils.h"

ParsedCountsY::ParsedCountsY() {
	rowDim = 0;
	colDim = 0;
}

bool ParsedCountsY::addGene(std::string geneName,
		std::vector<double> countsY){
	if (countsY.empty()){
		return false;
	}
	if (colDim == 0){
		this->parsedY.insert(std::pair<std::string, std::vector<double> >(geneName, countsY));
		colDim = countsY.size();
		rowDim++;
		return true;
	} else if (colDim == countsY.size()){
		if (parsedY.find(geneName) == parsedY.end()) {
			this->parsedY.insert(std::pair<std::string, std::vector<double> >(geneName, countsY));
			rowDim++;
		} else {
			timed_log("Duplicate entries in Y file, please modify your file.\nNow exiting.\n");
			exit(1);
		}

		return true;
	} else {
		return false;
	}
}
size_t ParsedCountsY::getColDim() {
	return this->colDim;
}

size_t ParsedCountsY::getRowDim() {

	return this->rowDim;
}

std::vector<std::string> ParsedCountsY::getGeneNames() {
	std::pair<std::string, std::vector<double> > gene;
	std::vector<std::string> geneNames;
    std::map<std::string, std::vector<double> >::iterator it;
    for (it = this->parsedY.begin(); it != this->parsedY.end(); ++it) {
        geneNames.push_back(it->first);
    }
//	BOOST_FOREACH(gene, this->parsedY){
//		geneNames.push_back(gene.first);
//	}
	return geneNames;
}

std::vector<double> ParsedCountsY::getCountsFromAllSamplesOfGene(
		std::string geneName) {

	return this->parsedY.at(geneName);
}

void ParsedCountsY::removeCols(std::vector<size_t> idxCols) {
	if (idxCols.empty()){return;}
	std::pair<std::string, std::vector<double> > geneBeforeElementRemoval;
	std::vector<double> countsBuffer;
	std::map<std::string, std::vector<double> > newParsedY;

    std::map<std::string, std::vector<double> >::iterator it;

    for (it = this->parsedY.begin(); it != this->parsedY.end(); ++it) {
        countsBuffer.clear();
        for (size_t i = 0; i < it->second.size(); ++i) {
            if (std::find(idxCols.begin(), idxCols.end(), i) == idxCols.end()) {
                countsBuffer.push_back(it->second.at(i));
            }
        }
        newParsedY.insert(std::pair<std::string, std::vector<double> >(it->first, countsBuffer));
    }

//	BOOST_FOREACH(geneBeforeElementRemoval, this->parsedY){
//		countsBuffer.clear();
//		for(size_t i=0; i<geneBeforeElementRemoval.second.size(); ++i){
//			if (std::find(idxCols.begin(),idxCols.end(),i)==idxCols.end()){
//				countsBuffer.push_back(geneBeforeElementRemoval.second.at(i));
//			}
//		}
//		newParsedY.insert(std::pair<std::string, std::vector<double> >(geneBeforeElementRemoval.first, countsBuffer));
//	}

	this->parsedY = newParsedY;
	this->colDim = countsBuffer.size();
}

void ParsedCountsY::initY(std::string fileNameY) {
	timed_log("Parsing "+fileNameY+".");
	std::fstream file_in(fileNameY);
	if (!file_in.good()){
				timed_log("Cannot access "+fileNameY+".\nNow exiting.");
				exit(1);
			}
	std::string line;
	while(std::getline(file_in, line)){
		std::istringstream iss(line);
		std::string geneName;
		iss >> geneName;

		std::string geneCountsStr;
		iss >> geneCountsStr;

		std::istringstream iss2(geneCountsStr);
		std::string geneCountToken;

		std::vector<double> geneCountsArr;
		while(std::getline(iss2, geneCountToken, ',')){
			double geneCountDbl;
			try{
				geneCountDbl = std::stod(geneCountToken,NULL);
			} catch (const std::invalid_argument& ia){
				timed_log("Values in counts file (Y) are not numeric! Now exiting.");
				exit(1);
			} catch (const std::out_of_range& oor) {
				timed_log("Values in counts file (Y) are out of range! Now exiting.");
				exit(1);
			}
			geneCountsArr.push_back(geneCountDbl);
		}

		size_t non_zero_nums = 0;
		for (size_t i = 0; i < geneCountsArr.size(); ++i) {
			if(geneCountsArr.at(i)!=0){
				non_zero_nums++;
			}
		}

		if(non_zero_nums>2){
			if (!this->addGene(geneName, geneCountsArr)){
				timed_log("Sample sizes are not consistent across genes! Now exiting.");
				exit(1);
			}
		} else {
			timed_log(geneName + " has less than 3 non_zero counts. Ignoring.");
		}
	}
}

ParsedCountsY::~ParsedCountsY() {
	// TODO Auto-generated destructor stub
}

size_t ParsedCountsY::getNumGenes() {
	return this->rowDim;
}



std::pair<std::string, std::vector<double> > ParsedCountsY::getGeneNameAndCounts(std::string geneName) {
	std::pair<std::string, std::vector<double> > ptr;
	ptr.first = geneName;
	ptr.second = this->getCountsFromAllSamplesOfGene(geneName);
	return ptr;
}

std::vector<double> ParsedCountsY::computeCellTotalReadCounts() {
	std::vector<double> results(this->getColDim());
	for(std::map<std::string, std::vector<double> >::iterator it = this->parsedY.begin(); it !=this->parsedY.end(); ++it){
		for(size_t i=0; i<it->second.size(); ++i){
			results[i] = results[i] + it->second.at(i);
		}
	}
/*    for (auto& gene : this->parsedY){
        for(size_t i=0; i<gene.second.size(); ++i){
            results[i] = results[i] + gene.second.at(i);
        }
    }*/
    this->cellTotalCounts = results;
    return this->cellTotalCounts;
}
