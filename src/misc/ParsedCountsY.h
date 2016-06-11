/*
 * ParsedCountsY.h
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */

#ifndef MISC_PARSEDCOUNTSY_H_
#define MISC_PARSEDCOUNTSY_H_

#include <map>
#include <string>
#include <vector>
#include <memory>
class ParsedCountsY {
private:
	std::map<std::string, std::vector<double> > parsedY;
	size_t rowDim;
	size_t colDim;
	bool addGene(std::string geneName, std::vector<double>);
	std::vector<double> cellTotalCounts;
public:
	ParsedCountsY();
	void initY(std::string fileNameY);
	size_t getColDim();
	size_t getRowDim();
	std::vector<double> computeCellTotalReadCounts();
	size_t getNumGenes();
	std::vector<std::string> getGeneNames();
	std::vector<double> getCountsFromAllSamplesOfGene(std::string geneName);
	void removeCols(std::vector<size_t> idxCols);
	std::pair<std::string, std::vector<double> > getGeneNameAndCounts(std::string geneName);
	virtual ~ParsedCountsY();
};

#endif /* MISC_PARSEDCOUNTSY_H_ */
