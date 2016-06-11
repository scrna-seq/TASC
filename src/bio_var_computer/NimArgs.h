//
// Created by jiacheng on 11/12/15.
//

#ifndef TASC_NIMARGS_H
#define TASC_NIMARGS_H

#include <vector>
#include <stddef.h>


class NimArgs {
public:
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> kappa;
    std::vector<double> tau;
    std::vector<double> y;
    size_t sample_size;
    double epsabs;
    NimArgs(std::vector<double> alpha,
            std::vector<double> beta,
            std::vector<double> kappa,
            std::vector<double> tau,
            std::vector<double> y,
            size_t sample_size,
            double epsabs);
    virtual ~NimArgs();

};


#endif //TASC_NIMARGS_H
