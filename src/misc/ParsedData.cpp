/*
 * ParsedData.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: jiacheng
 */

#include <gsl/gsl_linalg.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <sstream>

#include "../misc/utils.h"
#include "../cubature/cubature.h"
#include "../param_estimation/InitMargLklhoodArgs.h"
#include "../param_estimation/init_integrand.h"
#include "../param_estimation/IntegrationArgs.h"
#include "../param_estimation/logistic_regression_irls.h"
#include "ParsedData.h"

ParsedData::ParsedData() {
    this->isGivenABKT = false;
    this->sampleSize = 0;
    this->initialHyperParamValues.e_kappa = -2;
    this->initialHyperParamValues.e_tau = 1;
    this->initialHyperParamValues.sd_kappa = 1;
    this->initialHyperParamValues.sd_tau = 1;
    this->initialHyperParamValues.rho_kappa_tau = transform_rho_from_real_domain(-1);
    // TODO Auto-generated constructor stub
}

void ParsedData::synchronizeAllParsedData(std::vector<size_t> invalidSampleIdx) {
    verbose_timed_log("Invalid samples: " + getStrOfStdVector(invalidSampleIdx));
    this->yParser.removeCols(invalidSampleIdx);
    this->xParser.removeCols(invalidSampleIdx);
    if (!this->isGivenABKT){
        this->erccParser.removeCols(invalidSampleIdx);
        this->alpha = this->erccParser.getAlpha();
        this->beta = this->erccParser.getBeta();
        this->sampleSize = this->alpha.size();
    } else {
        std::vector<double> alpha_after;
        std::vector<double> beta_after;
        std::vector<double> kappa_after;
        std::vector<double> tau_after;
        for (size_t i = 0; i < this->alpha.size(); ++i) {
            if(std::find(invalidSampleIdx.begin(),invalidSampleIdx.end(),i)==invalidSampleIdx.end()){
                alpha_after.push_back(this->alpha.at(i));
            }
        }

        for (size_t i = 0; i < this->beta.size(); ++i) {
            if(std::find(invalidSampleIdx.begin(),invalidSampleIdx.end(),i)==invalidSampleIdx.end()){
                beta_after.push_back(this->beta.at(i));
            }
        }

        for (size_t i = 0; i < this->kappa.size(); ++i) {
            if(std::find(invalidSampleIdx.begin(),invalidSampleIdx.end(),i)==invalidSampleIdx.end()){
                kappa_after.push_back(this->kappa.at(i));
            }
        }

        for (size_t i = 0; i < this->tau.size(); ++i) {
            if(std::find(invalidSampleIdx.begin(),invalidSampleIdx.end(),i)==invalidSampleIdx.end()){
                tau_after.push_back(this->tau.at(i));
            }
        }

        this->alpha = alpha_after;
        this->beta = beta_after;
        this->kappa = kappa_after;
        this->tau = tau_after;
    }

}

double safeConvertStrToDbl(std::string s, std::string valueMsg) {
    double r = 0;
    try {
        r = std::stod(s, NULL);
    } catch (const std::invalid_argument &ia) {
        timed_log("Values in " + valueMsg + " are not numeric! Now exiting.");
        exit(1);
    } catch (const std::out_of_range &oor) {
        timed_log("Values in " + valueMsg + " are out of range! Now exiting.");
        exit(1);
    }
    return r;
}

void ParsedData::parseABKT(std::string abktFileName) {
    std::string paramValueToken;
    std::fstream file_in(abktFileName);
    if (!file_in.good()) {
        timed_log("Cannot access " + abktFileName + ".\nNow exiting.");
        exit(1);
    }
    std::string line;
    while (std::getline(file_in, line)) {
        std::istringstream iss(line);
        iss >> paramValueToken;
        this->alpha.push_back(safeConvertStrToDbl(paramValueToken, "ABKT File"));
        iss >> paramValueToken;
        this->beta.push_back(safeConvertStrToDbl(paramValueToken, "ABKT File"));
        iss >> paramValueToken;
        this->kappa.push_back(safeConvertStrToDbl(paramValueToken, "ABKT File"));
        iss >> paramValueToken;
        this->tau.push_back(safeConvertStrToDbl(paramValueToken, "ABKT File"));
    }

}

void ParsedData::parseAndEstimateHyperParams(CLIArgs runningParams) {
    this->disableAdjustmentForSize = runningParams.isDisableSizeAdjustment();
    this->xParser.initX(runningParams.getXFileName(), runningParams.getCovaIdxStr());
    this->yParser.initY(runningParams.getYFileName());

    std::string abkt_file_name = runningParams.getAbktEstimatesFileName();

    if (!abkt_file_name.empty()) {
        this->isGivenABKT=true;
        this->disableAdjustmentForSize = true;
        timed_log("abkt_file_name not empty. ignoring ercc filename. ");
        struct stat buffer;
        if (stat(abkt_file_name.c_str(), &buffer) != 0) {
            timed_log("abkt_file_name not valid. Check your file name and run again. Now exiting.");
            exit(1);
        } else {
            timed_log("parsing " + abkt_file_name + ".");
            this->parseABKT(abkt_file_name);
            std::vector<size_t> sampleSizes(6);
            sampleSizes[0] = this->xParser.getNumSamples();
            sampleSizes[1] = this->yParser.getColDim();
            sampleSizes[2] = this->alpha.size();
            sampleSizes[3] = this->beta.size();
            sampleSizes[4] = this->kappa.size();
            sampleSizes[5] = this->tau.size();
            //this->erccParser.parseFromFile(runningParams.getErccFileName());

            if (!std_vector_all_equal(sampleSizes, &(this->sampleSize))) {
                timed_log("Inconsistent sample sizes across x,y, and abkt. Now exiting.");
                exit(1);
            }
        }
    } else {
        this->erccParser.initERCC(runningParams.getErccFileName());
        this->alpha = this->erccParser.getAlpha();
        this->beta = this->erccParser.getBeta();
        timed_log(this->alpha, "Raw Alpha:");
        timed_log(this->beta, "Raw Beta:");
        if (!this->isEqualSampleSize()) {
            timed_log("Sample sizes not consistent across all files! Now exiting.");
            exit(1);
        }

        this->synchronizeAllParsedData(this->erccParser.getInvalidSampleIdx());
        this->computeKappaTau(runningParams);
        if (!prog_params.isLogistic()) {
            this->computeEmpiricalBayesMean(runningParams);
        }
    }

    std::ostream *outputDevice;
    std::ofstream abktFileOut;
    bool printFileName = false;
    try {
        abktFileOut.open(runningParams.getOutputFileName() + ".abkt");
        outputDevice = &abktFileOut;
        timed_log("Writing estimated abkt to file: " + runningParams.getOutputFileName() + ".abkt");
    } catch (std::ofstream::failure &writeErr) {
        timed_log("Cannot write to " + runningParams.getOutputFileName() + ".abkt. Now redirecting to stdout.");
        outputDevice = &std::cout;
        (*outputDevice) << "ABKT Estimates:" << std::endl;
    }

    std::vector<size_t> invalid_idx;
    this->cellTotalCounts = this->yParser.computeCellTotalReadCounts();

    if (!this->disableAdjustmentForSize) {

        this->cellTotalCountsERCC = this->erccParser.computeTotalCellCountsERCC();

        if ((this->cellTotalCounts.size() != this->sampleSize) || (cellTotalCountsERCC.size() != this->sampleSize)) {
            timed_log("this->sampleSize," + std::to_string((long long) this->sampleSize));
            timed_log("this->cellTotalCounts.size()," + std::to_string((long long) this->cellTotalCounts.size()));
            timed_log("this->cellTotalCountsERCC.size()," + std::to_string((long long) cellTotalCountsERCC.size()));
            timed_log("Inconsistent sample size with cell total counts computation. Now exiting.");
            exit(1);
        }

        for (size_t j = 0; j < this->sampleSize; ++j) {
            if ((this->cellTotalCounts[j] == 0) || (cellTotalCountsERCC[j] == 0)) {
                invalid_idx.push_back(j);
            } else {
                this->sizeFactors.push_back(std::log(this->cellTotalCounts[j] / cellTotalCountsERCC[j]));
            }
        }

        double min_sf = *std::min_element(this->sizeFactors.begin(), this->sizeFactors.end());

        for (size_t l = 0; l < sizeFactors.size(); ++l) {
            this->sizeFactors[l] = this->sizeFactors[l] - min_sf;
        }

        timed_log(this->sizeFactors, "Size factors");

        this->xParser.addSizeFactorCol(this->sizeFactors);
        timed_log("Adding size factor covariate to the design matrix.");
    } else {

        for (size_t j = 0; j < this->sampleSize; ++j) {
            if (this->cellTotalCounts[j] == 0) {
                invalid_idx.push_back(j);
            }
        }

        timed_log("Size factor adjustment is disabled.");
    }

    timed_log("Removing invalid samples.");
    if (!invalid_idx.empty()) {

        this->synchronizeAllParsedData(invalid_idx);
    }


    (*outputDevice).precision(std::numeric_limits<double>::max_digits10);
    (*outputDevice) << "Alpha\tBeta\tKappa\tTau" << std::endl;

    for (size_t i = 0; i < this->sampleSize; i++) {
        (*outputDevice) << this->alpha[i] << "\t"
        << this->beta[i] << "\t"
        << this->kappa[i] << "\t"
        << this->tau[i] << std::endl;
    }

    timed_log("Total Number of Genes Parsed:\t" + std::to_string((long long int) (this->getNumGenes())));

    try {
        abktFileOut.close();
    } catch (std::ofstream::failure &e) { }

}

void ParsedData::computeKappaTau(CLIArgs runningParams) {
    timed_log("Computing Kappa Tau.");
    if (!erccParser.computeTransposedErccCounts()) {
        timed_log("Error transposing ercc counts matrix! Now exiting.");
        exit(1);
    }

    //compute Ekappa Etau Sigma_kappa Sigma_tau and Rho_kappa_tau
    gsl_vector *x = gsl_vector_alloc(5);
    gsl_vector_set(x, 0, this->initialHyperParamValues.e_kappa);
    gsl_vector_set(x, 1, this->initialHyperParamValues.e_tau);
    gsl_vector_set(x, 2, log(this->initialHyperParamValues.sd_kappa));
    gsl_vector_set(x, 3, log(this->initialHyperParamValues.sd_tau));
    gsl_vector_set(x, 4, transform_rho_to_real_domain(this->initialHyperParamValues.rho_kappa_tau));


    this->marginalLikelihoodArguments.setData(this->erccParser.getErccCountsTranspose(),
                                              this->erccParser.getErccTrueConcentration(),
                                              this->erccParser.getAlpha(),
                                              this->erccParser.getBeta(),
                                              this->erccParser.getErccLengths(),
                                              runningParams.getCubatureMaxIter(),
                                              runningParams.getMinRelTolCubature()
    );


    std::vector<double> kappa_is;
    std::vector<double> tau_is;
    bool converge = false;
    gsl_vector *beta_i = gsl_vector_alloc(2);
    gsl_vector *y_hat_i = gsl_vector_alloc(this->erccParser.getNumErccGenes());
    gsl_vector *y_counts = gsl_vector_alloc(this->erccParser.getNumErccGenes());
    gsl_matrix *x_mat = gsl_matrix_alloc(this->erccParser.getNumErccGenes(), 2);

    for (size_t j = 0; j < this->erccParser.getNumErccGenes(); ++j) {
        gsl_matrix_set(x_mat, j, 1, std::log(this->erccParser.getErccTrueConcentration().at(j)));
        gsl_matrix_set(x_mat, j, 0, 1);
    }

    if (prog_params.isLogistic()) {
        std::vector<size_t> failedSampleIdx;
        for (size_t i = 0; i < this->sampleSize; ++i) {
            for (size_t j = 0; j < this->erccParser.getNumErccGenes(); ++j) {
                gsl_vector_set(y_counts, j, this->marginalLikelihoodArguments.y.at(i).at(j) != 0);
            }
            int status = logistic_regression(x_mat, y_counts, &beta_i, &y_hat_i, 100, 1e-9, &converge);

            if (!status) {
                this->kappa.push_back(gsl_vector_get(beta_i, 0));
                this->tau.push_back(gsl_vector_get(beta_i, 1));
            } else {
                failedSampleIdx.push_back(i);
            }
        }

        std::stringstream failureSampleMsgBuilder;
        failureSampleMsgBuilder << failedSampleIdx.size() <<
        " samples failed in estimating Kappa and Tau. Their indexes are: ";
        for (size_t i = 1; i < failedSampleIdx.size(); ++i) {
            failureSampleMsgBuilder << failedSampleIdx.at(i) << ", ";
        }
        timed_log(failureSampleMsgBuilder.str());
        timed_log("Removing above samples from downstream analysis.");
        this->synchronizeAllParsedData(failedSampleIdx);
    } else {
        for (size_t i = 0; i < this->sampleSize; ++i) {
            for (size_t j = 0; j < this->erccParser.getNumErccGenes(); ++j) {
                gsl_vector_set(y_counts, j, this->marginalLikelihoodArguments.y.at(i).at(j) != 0);
            }
            int status = logistic_regression(x_mat, y_counts, &beta_i, &y_hat_i, 100, 1e-9, &converge);

            if (!status) {
                kappa_is.push_back(gsl_vector_get(beta_i, 0));
                tau_is.push_back(gsl_vector_get(beta_i, 1));
            }
        }


        fitBivariateNormal(kappa_is, tau_is, &(this->estimatedHyperParamValues.e_kappa),
                           &(this->estimatedHyperParamValues.e_tau),
                           &(this->estimatedHyperParamValues.sd_kappa),
                           &(this->estimatedHyperParamValues.sd_tau),
                           &(this->estimatedHyperParamValues.rho_kappa_tau));


        timed_log("EK = " + std::to_string((long double) this->estimatedHyperParamValues.e_kappa));
        timed_log("ET = " + std::to_string((long double) this->estimatedHyperParamValues.e_tau));
        timed_log("SD_K = " + std::to_string((long double) this->estimatedHyperParamValues.sd_kappa));
        timed_log("SD_T = " + std::to_string((long double) this->estimatedHyperParamValues.sd_tau));
        timed_log("Rho_KT = " + std::to_string((long double) this->estimatedHyperParamValues.rho_kappa_tau));
    }


//	timed_log("Minimum LogLikelihood = "+std::to_string((long double) -s->fval));
}

bool ParsedData::isEqualSampleSize() {
    std::vector<size_t> sampleSizes(3);
    sampleSizes[0] = this->xParser.getNumSamples();
    sampleSizes[1] = this->yParser.getColDim();
    sampleSizes[2] = this->erccParser.getNumSamples();

    return isAllEqual(sampleSizes);
}

void ParsedData::computeEmpiricalBayesMean(CLIArgs runningParams) {
    size_t sample_size = this->marginalLikelihoodArguments.y.size();
    int cubature_status_kappa_integr[sample_size];
    int cubature_status_tau_integr[sample_size];
    int cubature_status_marginal_lklhood[sample_size];
    double val_kappa[sample_size];
    double val_tau[sample_size];
    double val_marginal[sample_size];
    double err_kappa[sample_size];
    double err_tau[sample_size];
    double err_marginal[sample_size];
    double xmin[2] = {-1, -1}, xmax[2] = {1, 1};

    timed_log("Computing empirical Bayesian estimates for kappa and tau.");

    //initialize a container of indicators of whether integration fails at a specific sample.
    bool isSampleFailed[sample_size];
    for (size_t i = 0; i < sample_size; ++i) {
        isSampleFailed[i] = false;
    }

    //initialize containers for the idx of failed samples.
    std::vector<size_t> failedSampleIdx;

#pragma omp parallel for schedule(dynamic,1)
    for (size_t i = 0; i < sample_size; ++i) {
        IntegrationArgs single_lklhood_args;
        single_lklhood_args.SetAllArgs(
                this->marginalLikelihoodArguments.alphas[i],
                this->marginalLikelihoodArguments.betas[i],
                this->marginalLikelihoodArguments.y[i],
                this->marginalLikelihoodArguments.true_concentrations,
                this->marginalLikelihoodArguments.ercc_lens,
                this->estimatedHyperParamValues.e_kappa,
                this->estimatedHyperParamValues.e_tau,
                this->estimatedHyperParamValues.sd_kappa,
                this->estimatedHyperParamValues.sd_tau,
                this->estimatedHyperParamValues.rho_kappa_tau);
        unsigned fdim = 1;
        unsigned dim = 2;

        cubature_status_kappa_integr[i] = hcubature(fdim,
                                                    kappa_times_pdf_single_cell,
                                                    (void *) (&single_lklhood_args),
                                                    dim,
                                                    xmin,
                                                    xmax,
                                                    runningParams.getCubatureMaxIter(),
                                                    0,
                                                    runningParams.getMinRelTolCubature(),
                                                    ERROR_INDIVIDUAL,
                                                    val_kappa + i,
                                                    err_kappa + i);

        cubature_status_tau_integr[i] = hcubature(fdim,
                                                  tau_times_pdf_single_cell,
                                                  (void *) (&single_lklhood_args),
                                                  dim,
                                                  xmin,
                                                  xmax,
                                                  runningParams.getCubatureMaxIter(),
                                                  0,
                                                  runningParams.getMinRelTolCubature(),
                                                  ERROR_INDIVIDUAL,
                                                  val_tau + i,
                                                  err_tau + i);

        cubature_status_marginal_lklhood[i] = hcubature(fdim,
                                                        pdf_single_cell,
                                                        (void *) (&single_lklhood_args),
                                                        dim,
                                                        xmin,
                                                        xmax,
                                                        runningParams.getCubatureMaxIter(),
                                                        0,
                                                        runningParams.getMinRelTolCubature(),
                                                        ERROR_INDIVIDUAL,
                                                        val_marginal + i,
                                                        err_marginal + i);

        if (cubature_status_kappa_integr[i] || cubature_status_tau_integr[i] || cubature_status_marginal_lklhood[i]) {
            val_kappa[i] = NAN;
            err_kappa[i] = NAN;
            val_tau[i] = NAN;
            err_tau[i] = NAN;
            val_marginal[i] = NAN;
            err_marginal[i] = NAN;
            isSampleFailed[i] = true;
            continue;
        }
    }

    for (size_t i = 0; i < sample_size; ++i) {
        if (isSampleFailed[i]) {
            failedSampleIdx.push_back(i);
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " E_KAPPA:" +
                              std::to_string((long double) val_kappa[i]));
            verbose_timed_log(
                    "Sample #" + std::to_string((long long) i) + " E_TAU:" + std::to_string((long double) val_tau[i]));
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " MARGINAL:" +
                              std::to_string((long double) val_marginal[i]));
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " failed in estimating Kappa and Tau.");
        } else {
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " E_KAPPA:" +
                              std::to_string((long double) val_kappa[i]));
            verbose_timed_log(
                    "Sample #" + std::to_string((long long) i) + " E_TAU:" + std::to_string((long double) val_tau[i]));
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " MARGINAL:" +
                              std::to_string((long double) val_marginal[i]));
            verbose_timed_log("Sample #" + std::to_string((long long) i) + " succeeded in estimating Kappa and Tau.");
            this->kappa.push_back(val_kappa[i] / val_marginal[i]);
            this->tau.push_back(val_tau[i] / val_marginal[i]);
        }
    }

    std::stringstream failureSampleMsgBuilder;
    failureSampleMsgBuilder << failedSampleIdx.size() <<
    " samples failed in estimating Kappa and Tau. Their indexes are: ";
    for (size_t i = 1; i < failedSampleIdx.size(); ++i) {
        failureSampleMsgBuilder << failedSampleIdx.at(i) << ", ";
    }
    timed_log(failureSampleMsgBuilder.str());
    timed_log("Removing above samples from downstream analysis.");
    this->synchronizeAllParsedData(failedSampleIdx);
}

std::vector<std::string> ParsedData::getGeneNames() {
    return this->yParser.getGeneNames();
}

std::vector<double> ParsedData::getCountsOfGene(std::string geneName) {
    return this->yParser.getCountsFromAllSamplesOfGene(geneName);
}

gsl_matrix *ParsedData::getDesignMatrix() {
    return this->xParser.getDesignMatrix();
}

std::vector<std::string> ParsedData::getCovNames() {
    return this->xParser.covariateNames;
}

std::vector<size_t> ParsedData::getTestColIdx() {
    return this->xParser.idxCovToBeTested;
}

ParsedData::~ParsedData() {
    // TODO Auto-generated destructor stub
}


size_t ParsedData::getNumGenes() {
    return this->yParser.getNumGenes();
}

std::pair<std::vector<std::string>, gsl_matrix *> ParsedData::getCovNamesAndDesignMat() {
    return this->xParser.getCovNamesAndDesignMatrix();
}

std::pair<std::string, std::vector<double> > ParsedData::getGeneNameAndCounts(
        std::string geneName) {
    return this->yParser.getGeneNameAndCounts(geneName);
}

std::vector<size_t> ParsedData::getCovariatesToBeTested() {
    return this->xParser.idxCovToBeTested;
}

std::vector<double> ParsedData::getAlpha() {
    return this->alpha;
}

std::vector<double> ParsedData::getBeta() {
    return this->beta;
}

std::vector<double> ParsedData::getKappa() {
    return this->kappa;
}

std::vector<double> ParsedData::getTau() {
    return this->tau;
}

std::vector<double> ParsedData::getCellTotalCounts() {
    return this->cellTotalCounts;
}

std::vector<double> ParsedData::getCellTotalCountsERCC() {
    return this->cellTotalCountsERCC;
}

std::vector<double> ParsedData::getSizeFactors() {
    return this->sizeFactors;
}
