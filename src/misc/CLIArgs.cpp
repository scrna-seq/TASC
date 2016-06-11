/*
 * ProgramParameters.cpp
 *
 *  Created on: Sep 14, 2015
 *      Author: cheng
 */

#include <iostream>
#include "CLIArgs.h"

CLIArgs::CLIArgs() {
    num_threads= 1;
    multi_threading= false;
    laplace= false;
    max_iter_em= 50;
    min_rel_tol_em= 1e-5;
    max_iter_nim= 500;
    min_rel_tol_nim= 1e-5;
    max_iter_ercc= 500;
    min_rel_tol_ercc= 1e-5;
    init_sigma= 0.2;
    run_de= false;
    run_vd= false;
    run_im= false;
    verbose= false;
    cubature_max_iter= 100000;
    min_rel_tol_cubature=1e-7;
    scale_factor_base= 0.5;
    size_t maxIterCQuad= 2000;
    disableSizeAdjustment= false;
}
CLIArgs::~CLIArgs() {
}

int CLIArgs::setAlmostAllParams(int num_threads,
								bool multi_threading,
								const std::string y_file_name,
								const std::string x_file_name,
								const std::string output_file_name,
								const std::string ercc_file_name,
								const std::string cova_idx_str,
								int max_iter_em,
								int max_iter_nim,
								int max_iter_ercc,
								double min_rel_tol_em,
								double min_rel_tol_nim,
								double min_rel_tol_ercc,
								double init_sigma,
								bool run_de,
								bool run_vd,
								bool run_im,
								bool verbose){
	this->num_threads=num_threads;
	this->multi_threading=multi_threading;
	this->y_file_name=y_file_name;
	this->x_file_name=x_file_name;
	this->output_file_name=output_file_name;
	this->ercc_file_name=ercc_file_name;
	this->cova_idx_str=cova_idx_str;
	this->max_iter_em=max_iter_em;
	this->max_iter_nim=max_iter_nim;
	this->max_iter_ercc=max_iter_ercc;
	this->min_rel_tol_em=min_rel_tol_em;
	this->min_rel_tol_nim=min_rel_tol_nim;
	this->min_rel_tol_ercc=min_rel_tol_ercc;
	this->init_sigma=init_sigma;
	this->run_de = run_de;
	this->run_vd = run_vd;
	this->run_im = run_im;
	this->verbose=verbose;
	return 0;
}

std::string CLIArgs::getCovaIdxStr()  {
	return cova_idx_str;
}

unsigned int CLIArgs::getCubatureMaxIter()  {
	return cubature_max_iter;
}

void CLIArgs::setCubatureMaxIter(unsigned int cubatureMaxIter) {
	cubature_max_iter = cubatureMaxIter;
}

std::string CLIArgs::getErccFileName()  {
	return ercc_file_name;
}

double CLIArgs::getInitSigma()  {
	return init_sigma;
}

unsigned int CLIArgs::getMaxIterEm()  {
	return max_iter_em;
}

unsigned int CLIArgs::getMaxIterErcc()  {
	return max_iter_ercc;
}

unsigned int CLIArgs::getMaxIterNim()  {
	return max_iter_nim;
}

void CLIArgs::setMaxIterNim(unsigned int maxIterNim) {
	max_iter_nim = maxIterNim;
}

double CLIArgs::getMinRelTolCubature()  {
	return min_rel_tol_cubature;
}

void CLIArgs::setMinRelTolCubature(double minRelTolCubature) {
	min_rel_tol_cubature = minRelTolCubature;
}

double CLIArgs::getMinRelTolEm()  {
	return min_rel_tol_em;
}

void CLIArgs::setMinRelTolEm(double minRelTolEm) {
	min_rel_tol_em = minRelTolEm;
}

double CLIArgs::getMinRelTolErcc()  {
	return min_rel_tol_ercc;
}

void CLIArgs::setMinRelTolErcc(double minRelTolErcc) {
	min_rel_tol_ercc = minRelTolErcc;
}

double CLIArgs::getMinRelTolNim()  {
	return min_rel_tol_nim;
}

void CLIArgs::setMinRelTolNim(double minRelTolNim ) {
	min_rel_tol_nim = minRelTolNim;
}

bool CLIArgs::isMultiThreading()  {
	return multi_threading;
}

void CLIArgs::setMultiThreading(bool multiThreading ) {
	multi_threading = multiThreading;
}

int CLIArgs::getNumThreads()  {
	return num_threads;
}

void CLIArgs::setNumThreads(int numThreads ) {
	num_threads = numThreads;
}

std::string CLIArgs::getOutputFileName()  {
	return output_file_name;
}

void CLIArgs::setOutputFileName( std::string outputFileName) {
	output_file_name = outputFileName;
}

bool CLIArgs::isRunVD()  {
	return run_vd;
}

void CLIArgs::setRunVD(bool runVD) {
	run_vd = runVD;
}

bool CLIArgs::isRunDE()  {
	return run_de;
}

void CLIArgs::setRunDE(bool runDE) {
	run_de = runDE;
}

bool CLIArgs::isVerbose()  {
	return verbose;
}

void CLIArgs::setVerbose(bool verbose ) {
	this->verbose = verbose;
}

std::string CLIArgs::getXFileName()  {
	return x_file_name;
}

void CLIArgs::setXFileName( std::string fileName) {
	x_file_name = fileName;
}

std::string CLIArgs::getYFileName()  {
	return y_file_name;
}

void CLIArgs::setYFileName( std::string fileName) {
	y_file_name = fileName;
}

std::string CLIArgs::getAbktEstimatesFileName()  {
	return abkt_estimates_file_name;
}

void CLIArgs::setAbktEstimatesFileName(
		std::string abktEstimatesFileName) {
	abkt_estimates_file_name = abktEstimatesFileName;
}

void CLIArgs::displayAllArgs(){
	std::cerr << "num_threads:\t" << num_threads << std::endl;
	std::cerr << "multi_threading:\t" << multi_threading << std::endl;
	std::cerr << "y_file_name:\t" << y_file_name << std::endl;
	std::cerr << "x_file_name:\t" << x_file_name << std::endl;
	std::cerr << "output_file_name:\t" << output_file_name << std::endl;
	std::cerr << "ercc_file_name:\t" << ercc_file_name << std::endl;
	std::cerr << "abkt_estimates_file_name:\t" << abkt_estimates_file_name << std::endl;
	std::cerr << "cova_idx_str:\t" << cova_idx_str << std::endl;

	std::cerr << "max_iter_em:\t" << max_iter_em << std::endl;
	std::cerr << "min_rel_tol_em:\t" << min_rel_tol_em << std::endl;
	std::cerr << "max_iter_nim:\t" << max_iter_nim << std::endl;
	std::cerr << "min_rel_tol_nim:\t" << x_file_name << std::endl;
	std::cerr << "max_iter_ercc:\t" << max_iter_ercc << std::endl;
	std::cerr << "min_rel_tol_ercc:\t" << min_rel_tol_ercc << std::endl;
	std::cerr << "init_sigma:\t" << init_sigma << std::endl;
	std::cerr << "run_de:\t" << run_de << std::endl;

	std::cerr << "run_vd:\t" << run_vd << std::endl;
	std::cerr << "verbose:\t" << verbose << std::endl;
	std::cerr << "cubature_max_iter:\t" << cubature_max_iter << std::endl;
	std::cerr << "min_rel_tol_cubature:\t" << min_rel_tol_cubature << std::endl;
}

double CLIArgs::getScaleFactorBase()  {
	return scale_factor_base;
}

void CLIArgs::setScaleFactorBase(double scaleFactorBase) {
	scale_factor_base = scaleFactorBase;
}

bool CLIArgs::isRunIM() {
	return this->run_im;
}

void CLIArgs::setRunIM(bool runIM) {
	this->run_im = runIM;
}

void CLIArgs::setLaplace(bool isLaplace) {
	this->laplace = isLaplace;
}

bool CLIArgs::isLaplace() {
	return this->laplace;
}

void CLIArgs::setMaxIterCQuad(size_t max_iter_cq) {
	this->maxIterCQuad = max_iter_cq;
}

size_t CLIArgs::getMaxIterCQuad() {
	return this->maxIterCQuad;
}

bool CLIArgs::isDisableSizeAdjustment() {
	return this->disableSizeAdjustment;
}

void CLIArgs::setDisableSizeAdjustment(bool dsa) {
	this->disableSizeAdjustment = dsa;
}
