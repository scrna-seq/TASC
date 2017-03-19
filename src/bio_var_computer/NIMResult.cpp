//
// Created by jiacheng on 11/12/15.
//

#include "NIMResult.h"
#include <cmath>

void NIMResult::print_all(std::ostream &outputFileHandle, bool verbose) {

    outputFileHandle << this->gene_name << "\t";
    outputFileHandle << this->aborted << "\t";
    outputFileHandle << this->theta_t << "\t";
    outputFileHandle << this->bioVar << "\t";
    outputFileHandle << this->varY << "\t";
    outputFileHandle << std::exp(this->sigma_t) << "\t";
    outputFileHandle << this->isConverged << "\t";
    outputFileHandle << this->logLikelihood << std::endl;
}

NIMResult::NIMResult(){};
NIMResult::~NIMResult(){};
