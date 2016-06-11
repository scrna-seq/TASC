/*
 * ParsedCovariatesX.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */


#include "boost/foreach.hpp"
#include <algorithm>
#include <sstream>
#include <fstream>

#include "ParsedCovariatesX.h"
#include "utils.h"

ParsedCovariatesX::ParsedCovariatesX() {
    this->numCovariates = 0;
    this->numSamples = 0;
    this->designMatrix = NULL;
}

void ParsedCovariatesX::initX(std::string fileNameX, std::string covIdxToBeTested) {
    timed_log("Parsing " + fileNameX + '.');
    std::fstream file_in(fileNameX);
    if (!file_in.good()) {
        timed_log("Cannot access " + fileNameX + ".\nNow exiting.");
        exit(1);
    }
    std::string line;
    std::vector<double> covariateValuesBuffer;
    std::string covariateNameBuffer;
    std::string covariateValuesStrBuffer;


    while (std::getline(file_in, line)) {

        covariateValuesBuffer.clear();
        std::istringstream iss(line);
        iss >> covariateNameBuffer;
        this->covariateNames.push_back(covariateNameBuffer);

        iss >> covariateValuesStrBuffer;
        std::istringstream iss2(covariateValuesStrBuffer);

        std::string covariateValueToken;
        while (std::getline(iss2, covariateValueToken, ',')) {
            double covariateValueDbl;
            try {
                covariateValueDbl = std::stod(covariateValueToken, NULL);
            } catch (const std::invalid_argument &ia) {
                timed_log("Values in covariate file (X) are not numeric! Now exiting.");
                exit(1);
            } catch (const std::out_of_range &oor) {
                timed_log("Values in covariate file (X) are out of range! Now exiting.");
                exit(1);
            }
            covariateValuesBuffer.push_back(covariateValueDbl);
        }

        if (!(this->insertSample(covariateNameBuffer, covariateValuesBuffer))) {
            timed_log("Sample sizes are not consistent across covariates! Now exiting.");
            exit(1);
        };
    }

    this->designMatrix = gsl_matrix_alloc(numSamples, numCovariates);

    for (size_t i = 0; i < this->numCovariates; ++i) {
        for (size_t j = 0; j < this->numSamples; ++j) {
            gsl_matrix_set(this->designMatrix, j, i, this->parsedX.at(this->covariateNames.at(i)).at(j));
        }
    }

    // now parse covIdxtobetested;
    // this is fucked up. from reading the downstream code,
    // i need this idx to be 1-based not 0-based!
    // separated by comma

    std::istringstream testedCovSS(covIdxToBeTested);
    std::string covIdxToken;
    while (std::getline(testedCovSS, covIdxToken, ',')) {
        long covIdx;
        try {
            covIdx = std::stol(covIdxToken);
        } catch (const std::invalid_argument &ia) {
            timed_log("Values in tested covariates Idx are not numeric! Now exiting.");
            exit(1);
        } catch (const std::out_of_range &oor) {
            timed_log("Values in tested covariates Idx are out of range! Now exiting.");
            exit(1);
        }
        this->idxCovToBeTested.push_back((size_t) covIdx);
    }

    for (size_t i = 0; i < this->idxCovToBeTested.size(); i++) {
        if (this->idxCovToBeTested[i] > this->numCovariates) {
            timed_log("Index for covariates to be tested out of range. Now exiting.");
            exit(1);
        }
    }
}

gsl_matrix *ParsedCovariatesX::getDesignMatrix() const {
    return designMatrix;
}

size_t ParsedCovariatesX::getNumCovariates() const {
    return numCovariates;
}

size_t ParsedCovariatesX::getNumSamples() const {
    return numSamples;
}

bool ParsedCovariatesX::insertSample(std::string covariateName, std::vector<double> covariateValues) {
    if (covariateValues.empty()) {
        return false;
    }
    if (numSamples == 0) {
        this->parsedX.insert(std::pair<std::string, std::vector<double> >(covariateName, covariateValues));
        numSamples = covariateValues.size();
        numCovariates++;
        return true;
    } else if (numSamples == covariateValues.size()) {
        this->parsedX.insert(std::pair<std::string, std::vector<double> >(covariateName, covariateValues));
        numCovariates++;
        return true;
    } else {
        return false;
    }
}

void ParsedCovariatesX::removeCols(std::vector<size_t> idxCols) {
    if (idxCols.empty()) { return; }
    std::pair<std::string, std::vector<double> > covariateBeforeElementRemoval;
    std::vector<double> covariateValueBuffer;
    std::map<std::string, std::vector<double> > newParsedX;

    std::map<std::string, std::vector<double> >::iterator it;

    for (it = this->parsedX.begin(); it != this->parsedX.end(); ++it) {
        covariateValueBuffer.clear();
        for (size_t i = 0; i < it->second.size(); ++i) {
            if (std::find(idxCols.begin(), idxCols.end(), i) == idxCols.end()) {
                covariateValueBuffer.push_back(it->second.at(i));
            }
        }
        newParsedX.insert(std::pair<std::string, std::vector<double> >(it->first,
                                                                       covariateValueBuffer));
    }

//    BOOST_FOREACH(covariateBeforeElementRemoval, this->parsedX) {
//                    covariateValueBuffer.clear();
//                    for (size_t i = 0; i < covariateBeforeElementRemoval.second.size(); ++i) {
//                        if (std::find(idxCols.begin(), idxCols.end(), i) == idxCols.end()) {
//                            covariateValueBuffer.push_back(covariateBeforeElementRemoval.second.at(i));
//                        }
//                    }
//                    newParsedX.insert(std::pair<std::string, std::vector<double> >(covariateBeforeElementRemoval.first,
//                                                                                   covariateValueBuffer));
//                }

    this->parsedX = newParsedX;
    this->numSamples = covariateValueBuffer.size();

    //now refresh the pointer of gsl_matrix
    gsl_matrix_free(this->designMatrix);
    this->designMatrix = gsl_matrix_alloc(this->numSamples, this->numCovariates);
    for (size_t i = 0; i < this->numCovariates; ++i) {
        for (size_t j = 0; j < this->numSamples; ++j) {
            gsl_matrix_set(this->designMatrix, j, i, this->parsedX.at(this->covariateNames.at(i)).at(j));
        }
    }

}

std::vector<double> ParsedCovariatesX::getCovariateValues(
        std::string covariateName) {
    return this->parsedX.at(covariateName);
}

ParsedCovariatesX::~ParsedCovariatesX() {
//	gsl_matrix_free(this->designMatrix);
}

std::pair<std::vector<std::string>, gsl_matrix *> ParsedCovariatesX::getCovNamesAndDesignMatrix() {
    std::pair<std::vector<std::string>, gsl_matrix *> ptr;
    ptr.first = this->covariateNames;
    ptr.second = this->getDesignMatrix();
    return ptr;
}

void ParsedCovariatesX::addSizeFactorCol(std::vector<double> sf) {
    if (sf.size() != this->numSamples) {
        timed_log("Inconsistent size for the size factor and the X covariates. Now exiting.");
        exit(1);
    }

    this->parsedX.insert(std::pair<std::string, std::vector<double> >("AutoGeneratedSizeFactor", sf));

    this->numCovariates += 1;
    this->covariateNames.push_back("AutoGeneratedSizeFactor");

    gsl_matrix_free(this->designMatrix);
    this->designMatrix = gsl_matrix_alloc(numSamples, numCovariates);

    for (size_t i = 0; i < this->numCovariates; ++i) {
        std::vector<double> colValues = this->parsedX.at(this->covariateNames.at(i));
        for (size_t j = 0; j < this->numSamples; ++j) {
            gsl_matrix_set(this->designMatrix, j, i, colValues.at(j));
        }
    }
}
