//
// Created by jiacheng on 11/12/15.
//

#ifndef TASC_NIMRESULT_H
#define TASC_NIMRESULT_H

#include <string>
#include <vector>
#include <iostream>

class NIMResult {
public:
    NIMResult();
    virtual ~NIMResult();
    std::string gene_name;
    double theta_t;
    double sigma_t;
    bool aborted;
    double logLikelihood;
    bool isConverged;
    double varY;
    double ratioOfBioVarToVarY;
    double bioVar;

    void print_all(std::ostream &outputDevice, bool verbose);
private:
};


#endif //TASC_NIMRESULT_H
