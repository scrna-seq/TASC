//
// Created by jiacheng on 11/12/15.
//

#include "NimArgs.h"

NimArgs::NimArgs(std::vector<double> alpha,
                 std::vector<double> beta,
                 std::vector<double> kappa,
                 std::vector<double> tau,
                 std::vector<double> y,
                 size_t sample_size,
                 double epsabs) {
    this->kappa = kappa;
    this->tau = tau;
    this->alpha = alpha;
    this->beta = beta;
    this->y = y;
    this->sample_size = sample_size;
    this->epsabs = epsabs;
}

NimArgs::~NimArgs() { };