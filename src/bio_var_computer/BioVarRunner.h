//
// Created by jiacheng on 11/12/15.
//

#ifndef TASC_BIOVARRUNNER_H
#define TASC_BIOVARRUNNER_H
#include <string>
#include <vector>
#include <stddef.h>
#include "NIMResult.h"
#include <cmath>
#include <limits>

class BioVarRunner {
public:
    BioVarRunner();
    void setData(std::pair<std::string, std::vector<double>> geneNameAndCounts, std::vector<double> alpha,
                 std::vector<double> beta, std::vector<double> kappa, std::vector<double> tau, double maxIter,
                 double minRelTol, double initSigma);
    void runModel();

    NIMResult getResult();
    virtual ~BioVarRunner();
private:
    NIMResult result;
    double epsabs;
    std::vector<double> scaleFactors;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> kappa;
    std::vector<double> tau;
    size_t maxIter;
    size_t sampleSize;
    double minRelTol;
    std::string geneName;
    std::vector<double> geneCounts;
    double initialSigma;



};


#endif //TASC_BIOVARRUNNER_H
