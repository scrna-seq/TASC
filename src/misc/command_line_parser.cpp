/*
 * command_line_parser.cpp
 *
 *  Created on: Sep 14, 2015
 *      Author: cheng
 */
#include "command_line_parser.h"
#include <boost/program_options.hpp>
#include "utils.h"

int parseCommandLine(CLIArgs& prog_params, int argc, char* argv[]){
	int num_threads = 1;
	bool multi_threading = false;
	std::string y_file_name;
	std::string x_file_name;
	std::string ercc_file_name;
	std::string output_file_name;
	std::string cova_idx_str;
	unsigned int max_iter_em = 50;
	double min_rel_tol_em = 1e-5;
	unsigned int max_iter_nim = 500;
	double min_rel_tol_nim = 1e-5;
	unsigned int max_iter_ercc = 500;
	double min_rel_tol_ercc = 1e-5;
	double init_sigma = 0.2;
	bool run_de = false;
	bool run_vd = false;
	bool run_im = false;
	bool verbose = false;
	std::string abkt_estimates_file_name;
	size_t max_iter_cubature_ercc = 100000;
	double min_rel_tol_cubature = 1e-16;
	double scale_factor_base = 0.98;
	bool laplace = false;
    size_t max_iter_cquad = 2000;
	bool disableSizeAdjustment = false;
	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
			("vd,a","run variance decomposition and variance ratio ranking")
			("de,b","run differential gene expression analysis")
			("laplace,f", "use laplace approximation instead of adaptive integration. faster but might be less accurate.")
			("disable_size_adjustment,e","disable adjustment for total sequencing counts.")
//			("em,a","run EM method")
			("num_threads,n",boost::program_options::value<int>(&num_threads)->default_value(1),"number of threads")
			("counts,y", boost::program_options::value<std::string>(&y_file_name)->default_value("../test/y.txt"), "file containing the read counts")
			("covariates,x", boost::program_options::value<std::string>(&x_file_name)->default_value("../test/x.txt"), "file containing the covariates")
			("ercc_spike_ins,k",boost::program_options::value<std::string>(&ercc_file_name)->default_value("../test/s.txt"),"file containing ercc spike-ins counts and true values")
			("abkt_estimates,w",boost::program_options::value<std::string>(&abkt_estimates_file_name)->default_value(""),"file containing the estimated alpha,beta,kappa,taus. Debug only.")
			("output_file,o", boost::program_options::value<std::string>(&output_file_name)->default_value("../test/result"), "output will be written into output.log, output.em, output.nim, output.abkt")
			("test_cols,t", boost::program_options::value<std::string>(&cova_idx_str)->default_value("2"), "indices of covariates to be tested")
//			("max_iter_em,m", boost::program_options::value<unsigned int>(&max_iter_em)->default_value(50),"maximum number of iterations for EM algorithm")
			("max_iter,i", boost::program_options::value<unsigned int>(&max_iter_nim)->default_value(500),"maximum number of iterations for NIM algorithm")
            ("max_iter_cquad,d", boost::program_options::value<size_t>(&max_iter_cquad)->default_value(2000),"maximum number of iterations for LibGSL CQUAD function.")
			("max_iter_ercc,r", boost::program_options::value<unsigned int>(&max_iter_ercc)->default_value(500),"maximum number of iterations for Estimation of Kappas and Taus")
			("max_iter_cubature_ercc,u", boost::program_options::value<size_t>(&max_iter_cubature_ercc)->default_value(100000),"maximum number of functional evaluations for Cubature Package")
//			("min_rel_tol_em,e", boost::program_options::value<double>(&min_rel_tol_em)->default_value(1e-5,"1e-5"),"maximum relative error for EM algorithm")
			("min_rel_tol,p", boost::program_options::value<double>(&min_rel_tol_nim)->default_value(1e-5, "1e-5"), "maximum relative error for NIM algorithm")
			("min_rel_tol_ercc,l", boost::program_options::value<double>(&min_rel_tol_ercc)->default_value(1e-8,"1e-8"),"maximum relative error for Estimation of Kappas and Taus")
			("min_rel_tol_cubature,q", boost::program_options::value<double>(&min_rel_tol_cubature)->default_value(1e-6,"1e-6"),"maximum relative error for Cubature Package")
			("init_sigma,s", boost::program_options::value<double>(&init_sigma)->default_value(0.2,"0.2"),"initial value of sigma")
//			("scale_factor_base,f",boost::program_options::value<double>(&scale_factor_base)->default_value(0.5,"0.5"),"base for the scale factors for computing marginal likelihood")
			("help,h", "display help message")
			("verbose,v","produce verbose debug info")
			("version","produce version info");

	boost::program_options::variables_map vm;
	if (argc<=1){
		std::cerr << desc << std::endl;
		exit(1);
	}
	try{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
		boost::program_options::notify(vm);
	} catch (boost::program_options::error &e){
		std::cerr << desc << std::endl;
		exit(1);
	}

	if (vm.count("disable_size_adjustment")) {
		disableSizeAdjustment = true;
	}
	if (vm.count("laplace")) {
		laplace = true;
	}
	if (vm.count("help")) {
		std::cerr << desc << std::endl;
		exit(1);
	}

	if (vm.count("de")) {
		run_de = true;
	}

	if (vm.count("vd")) {
		run_vd = true;
	}

	run_im = false;


	if (vm.count("version")){
		std::cerr << "Toolkit for Analysis of Single Cell data. TASC version 0.1."<< std::endl;
	}

	if (vm.count("verbose")) {
		verbose = true;
	}

	if (num_threads < 1) {
		timed_log("Invalid number of threads. Now exiting.");
		exit(1);
	} else if (num_threads == 1) {
		multi_threading = false;
	} else {
		multi_threading = true;
	}

	if (max_iter_em < 1 || max_iter_ercc < 1 || max_iter_nim < 1) {
		timed_log("Invalid max iterations. Now exiting.");
		exit(1);
	}

	if (min_rel_tol_em <= 0 || min_rel_tol_ercc <= 0 || min_rel_tol_nim <= 0) {
		timed_log("Invalid max relative error. Now exiting.");
		exit(1);
	}

	if(init_sigma<=0){
		timed_log("Invalid initial value of sigma. Now exiting.");
		exit(1);
	}


	prog_params.setAlmostAllParams(num_threads,
								   multi_threading,
								   y_file_name,
								   x_file_name,
								   output_file_name,
								   ercc_file_name,
								   cova_idx_str,
								   max_iter_em,
								   max_iter_nim,
								   max_iter_ercc,
								   min_rel_tol_em,
								   min_rel_tol_nim,
								   min_rel_tol_ercc,
								   init_sigma,
								   run_de,
								   run_vd,
								   run_im,
								   verbose);


	prog_params.setCubatureMaxIter(max_iter_cubature_ercc);
    prog_params.setMaxIterCQuad(max_iter_cquad);
	prog_params.setMinRelTolCubature(min_rel_tol_cubature);
	prog_params.setAbktEstimatesFileName(abkt_estimates_file_name);
	prog_params.setScaleFactorBase(scale_factor_base);
	prog_params.setLaplace(laplace);
	prog_params.setDisableSizeAdjustment(disableSizeAdjustment);
	return 0;
}



