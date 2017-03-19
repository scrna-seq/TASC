#include <iostream>
#include <vector>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <mutex>

#include "misc/utils.h"
#include "misc/command_line_parser.h"
#include "misc/ParsedCountsY.h"
#include "misc/ParsedERCC.h"
#include "misc/ParsedData.h"
#include "em_nim_algo/NIMThreads.h"
#include "bio_var_computer/BioVarComputer.h"

CLIArgs prog_params;
std::ostream *logFileHandle;
std::mutex log_handle_mutex;
int main(int argc, char *argv[]){

    //TODO: remove requirement for X when only quantification is requested

    //TODO: implement iterative method to perform the model fitting WITHOUT ERCC!

    gsl_set_error_handler_off();

    parseCommandLine(prog_params, argc, argv);
    std::ofstream logFileOut;
    try{
        logFileOut.open(prog_params.getOutputFileName()+".log");
        logFileHandle = &logFileOut;
        timed_log("Opened log file:"+prog_params.getOutputFileName()+".log for logging.");
    } catch (std::ofstream::failure &writeErr){
        logFileHandle = &std::cerr;
        timed_log("Cannot open log file:"+prog_params.getOutputFileName()+".log for writing. Now redirecting to stderr.");
    }
    omp_set_num_threads(prog_params.getNumThreads());
    timed_log("Using "+std::to_string((long long)prog_params.getNumThreads())+" threads for computation.");

    ParsedData dataParser;
    dataParser.parseAndEstimateHyperParams(prog_params);

    if (prog_params.isRunDE()){
        NIMThreads nimThreads(dataParser,prog_params);
        nimThreads.runThreads();
        nimThreads.printAllResults();
    }

    if (prog_params.isRunVD()) {
        BioVarComputer bioVarComputer(dataParser, prog_params);
        bioVarComputer.runThreads();
        bioVarComputer.printAllResults();
    }

    return 0;
}
