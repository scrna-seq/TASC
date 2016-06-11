/*
 * ProgramParameters.h
 *
 *  Created on: Sep 14, 2015
 *      Author: cheng
 */

#ifndef CLIARGS_H_
#define CLIARGS_H_
#include <string>
class CLIArgs {
public:
	CLIArgs();
	virtual ~CLIArgs();
	int setAlmostAllParams(int num_threads,
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
						   bool verbose);
	std::string getCovaIdxStr() ;

	unsigned int getCubatureMaxIter() ;
	void setCubatureMaxIter(unsigned int cubatureMaxIter);
	std::string getErccFileName() ;

	double getInitSigma() ;
	unsigned int getMaxIterEm() ;
	unsigned int getMaxIterErcc() ;
	unsigned int getMaxIterNim() ;
	void setMaxIterNim(unsigned int maxIterNim);
	double getMinRelTolCubature() ;
	void setMinRelTolCubature(double minRelTolCubature);
	double getMinRelTolEm() ;
	void setMinRelTolEm(double minRelTolEm );
	double getMinRelTolErcc() ;
	void setMinRelTolErcc(double minRelTolErcc);
	double getMinRelTolNim() ;
	void setMinRelTolNim(double minRelTolNim);
	bool isMultiThreading() ;
	void setMultiThreading(bool multiThreading);
	int getNumThreads() ;
	void setNumThreads(int numThreads);
	std::string getOutputFileName() ;
	void setOutputFileName( std::string outputFileName);
	bool isRunVD() ;
	void setRunVD(bool runVD);
	bool isRunDE() ;
	void setRunDE(bool runDE);
	bool isRunIM() ;
	void setRunIM(bool runIM);
	bool isVerbose() ;
	void setVerbose(bool verbose);
	std::string getXFileName() ;
	void setXFileName( std::string fileName);
	std::string getYFileName() ;
	void setYFileName( std::string fileName);
	std::string getAbktEstimatesFileName() ;
	void setAbktEstimatesFileName( std::string abktEstimatesFileName);
	void displayAllArgs();
	double getScaleFactorBase() ;
	void setScaleFactorBase(double scaleFactorBase);

	void setLaplace(bool isLaplace);
	void setMaxIterCQuad(size_t max_iter_cq);
	size_t getMaxIterCQuad();
	bool isLaplace();
	bool isDisableSizeAdjustment();
	void setDisableSizeAdjustment(bool dsa);

private:
	int num_threads;
	bool multi_threading;
	bool laplace;
	std::string y_file_name;;
	std::string x_file_name;;
	std::string output_file_name;;
	std::string ercc_file_name;;
	std::string cova_idx_str;;
	std::string abkt_estimates_file_name;;
	unsigned int max_iter_em;
	double min_rel_tol_em;
	unsigned int max_iter_nim;
	double min_rel_tol_nim;
	unsigned int max_iter_ercc;
	double min_rel_tol_ercc;
	double init_sigma;
	bool run_de;
	bool run_vd;
	bool run_im;
	bool verbose;
	unsigned int cubature_max_iter;
	double min_rel_tol_cubature;
	double scale_factor_base;
	size_t maxIterCQuad;
	bool disableSizeAdjustment;
};

#endif /* CLIARGS_H_ */
