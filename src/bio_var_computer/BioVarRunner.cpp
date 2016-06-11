//
// Created by jiacheng on 11/12/15.
//

#include "BioVarRunner.h"
#include "../misc/basic_math.h"
#include "../misc/utils.h"
#include "NimArgs.h"
#include "NimMarginalLikelihood.h"
#include <gsl/gsl_multimin.h>

BioVarRunner::BioVarRunner() {
    epsabs = 1e-7;
    maxIter=10000;
    sampleSize = 0;
    minRelTol=1e-16;
    initialSigma = 0.1;

};
BioVarRunner::~BioVarRunner() { };

double computeSampleVariance(std::vector<double> x) {
    if (x.size() < 2) {
        return 0;
    }
    double sumX = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sumX += x[i];
    }
    double meanX = sumX / (double) x.size();
    double sumSqs = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sumSqs += (x[i] - meanX) * (x[i] - meanX);
    }
    return sumSqs / (double) (x.size() - 1);
}

void BioVarRunner::setData(
        std::pair<std::string, std::vector<double> > geneNameAndCounts,
        std::vector<double> alpha,
        std::vector<double> beta,
        std::vector<double> kappa,
        std::vector<double> tau,
        double maxIter,
        double minRelTol,
        double initSigma) {
    this->alpha = alpha;
    this->beta = beta;
    this->kappa = kappa;
    this->tau = tau;
    this->maxIter = maxIter;
    this->minRelTol = minRelTol;
    this->geneName = geneNameAndCounts.first;
    this->geneCounts = geneNameAndCounts.second;
    this->initialSigma = initSigma;
    std::vector<size_t> sampleSizeArr(5);
    sampleSizeArr[0] = this->alpha.size();
    sampleSizeArr[1] = this->beta.size();
    sampleSizeArr[2] = this->kappa.size();
    sampleSizeArr[3] = this->tau.size();
    sampleSizeArr[4] = this->geneCounts.size();


    if(!isAllEqual(sampleSizeArr)){
        timed_log("Uneven sample size from X and Y. Now exiting.");
        exit(1);
    } else {
        this->sampleSize = this->geneCounts.size();
    }

    for (size_t i = 0; i < this->sampleSize; ++i) {
        this->scaleFactors.push_back(1);
    }
}

void BioVarRunner::runModel() {

    this->result.gene_name = this->geneName;
    std::vector<double> y = this->geneCounts;


    if (this->initialSigma <=0 ){
        timed_log("Exception: initial value of sigma_t <= 0!");
        exit(1);
    }

    NimArgs nimArgs(this->alpha,this->beta,this->kappa,this->tau,this->geneCounts,this->sampleSize,this->epsabs);
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *ix;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */
    ix = gsl_vector_alloc (2);
    gsl_vector_set(ix, 0, 1.0);
    gsl_vector_set(ix, 1, log(this->initialSigma));

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);

    /* Initialize method and iterate */
    minex_func.n = 2;
    minex_func.f = &nim_marginal_lklhood;
    minex_func.params = (void *)&nimArgs;

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, ix, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, epsabs);

    } while (status == GSL_CONTINUE && iter < this->maxIter);
    if ( status == GSL_SUCCESS || status == 0 ) {
        this->result.isConverged = true;
    } else {
        this->result.isConverged = false;
    }
    this->result.theta_t = gsl_vector_get(s->x,0);
    this->result.sigma_t = gsl_vector_get(s->x,1);
    this->result.logLikelihood = s->fval;

    this->result.varY = computeSampleVariance(this->geneCounts);
    this->result.bioVar = std::exp(this->result.sigma_t * 2) * std::exp(2 * this->result.theta_t);
    this->result.ratioOfBioVarToVarY = this->result.bioVar / this->result.varY;
    gsl_vector_free(ix);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    if (std::isnan(this->result.theta_t)){
        this->result.aborted=true;
    } else {
        this->result.aborted=false;
    }
}

NIMResult BioVarRunner::getResult() {
    return this->result;
}
