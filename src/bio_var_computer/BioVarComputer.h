//
// Created by jiacheng on 11/12/15.
//

#ifndef TASC_POSTPROCESSOR_H
#define TASC_POSTPROCESSOR_H

#include "../misc/ParsedData.h"
#include "../misc/CLIArgs.h"
#include "NIMResult.h"
#include "BioVarRunner.h"
#include <string>
#include <vector>

class BioVarComputer {
public:
    BioVarComputer(ParsedData parsedData, CLIArgs runningParam);

    void runThreads();

    void printAllResults();
    virtual ~BioVarComputer();

private:
    size_t numGenes;
    std::string outputFileName;
    size_t numThreads;
    std::vector<NIMResult> rawResults;
    std::vector<BioVarRunner> bioVarRunners;
    bool verbose;
};


#endif //TASC_POSTPROCESSOR_H
