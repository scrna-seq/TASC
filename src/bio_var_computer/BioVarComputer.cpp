//
// Created by jiacheng on 11/12/15.
//

#include <algorithm>
#include "BioVarComputer.h"
#include "omp.h"

BioVarComputer::BioVarComputer(ParsedData parsedData, CLIArgs runningParam) {
    this->numGenes = parsedData.getNumGenes();
    this->verbose = runningParam.isVerbose();
    this->outputFileName = runningParam.getOutputFileName() + ".vd";
    this->numThreads = prog_params.getNumThreads();
    this->rawResults = std::vector<NIMResult>(parsedData.getNumGenes());
    this->bioVarRunners = std::vector<BioVarRunner>(parsedData.getNumGenes());
    std::vector<std::string> geneNames = parsedData.getGeneNames();
    for (size_t i = 0; i < parsedData.getNumGenes(); ++i) {
        this->bioVarRunners[i] = BioVarRunner();
        this->bioVarRunners[i].setData(parsedData.getGeneNameAndCounts(geneNames[i]),
                                       parsedData.getAlpha(), parsedData.getBeta(),
                                       parsedData.getKappa(), parsedData.getTau(),
                                       prog_params.getMaxIterNim(), prog_params.getMinRelTolNim(),
                                       prog_params.getInitSigma());
    }
}


BioVarComputer::~BioVarComputer() { };


bool nimResultGreaterThan(const NIMResult &nimResult1, const NIMResult &nimResult2) {
    return (nimResult1.ratioOfBioVarToVarY > nimResult2.ratioOfBioVarToVarY);
};

void BioVarComputer::runThreads() {
    omp_set_num_threads(this->numThreads);
    timed_log("Starting NIM Algorithm with " + std::to_string((long long) numThreads) +
              " threads for variance decomposition.");
#pragma omp parallel for schedule(dynamic,1)
    for (size_t i = 0; i < this->numGenes; i++) {
        timed_log("Now running Gene No. "+std::to_string((long long)i));
        this->bioVarRunners[i].runModel();
        this->rawResults[i] = this->bioVarRunners[i].getResult();
    }

//
//    struct nimResultGreaterThan {
//        inline bool operator()(const NIMResult &nimResult1, const NIMResult &nimResult2) {
//            return (nimResult1.ratioOfBioVarToVarY > nimResult2.ratioOfBioVarToVarY);
//        }
//    };
//


    std::sort(this->rawResults.begin(), this->rawResults.end(), nimResultGreaterThan);
}

void BioVarComputer::printAllResults() {
    std::ofstream file_out;
    std::ostream *outputDevice;
    bool printFileName = false;

    try {
        file_out.open(this->outputFileName);
        outputDevice = &file_out;
    } catch (std::ofstream::failure &writeErr) {
        timed_log("Cannot write to " + this->outputFileName + ". Now redirecting to stdout.");
        outputDevice = &std::cout;
        printFileName = true;
    }
    if (printFileName) {
        std::cout << "NIM results:" << std::endl;
    }
    (*outputDevice) << "Gene Name" << "\t";
    (*outputDevice) << "NIM Aborted" << "\t";
    (*outputDevice) << "Mean Expression" << "\t";
    (*outputDevice) << "Biological Variation" << "\t";
    (*outputDevice) << "Total Variation" << "\t";
    (*outputDevice) << "Ratio Bio/Total" << "\t";
    (*outputDevice) << "Likelihood Converged" << "\t";
    (*outputDevice) << "Log Likelihood" << std::endl;
    for (size_t i = 0; i < this->numGenes; ++i) {
        this->rawResults[i].print_all(*outputDevice, this->verbose);
    }

    try {
        file_out.close();
    } catch (std::ofstream::failure &err) {

    }
}