cmake_minimum_required(VERSION 2.6.4)
project(TASC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -Wall -std=c++0x")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

find_package(LibGSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.35.0 COMPONENTS program_options REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

set(SOURCE_FILES
        src/cubature/clencurt.h
        src/cubature/converged.h
        src/cubature/cubature.h
        src/cubature/hcubature.cpp
        src/cubature/pcubature.cpp
        src/cubature/vwrapper.h
        src/misc/cheng_quad.h
        src/misc/cquad_const.cpp
        src/misc/cheng_quad.cpp
        src/em_nim_algo/EMResult.cpp
        src/em_nim_algo/EMResult.h
        src/em_nim_algo/EMRunner.cpp
        src/em_nim_algo/EMRunner.h
        src/em_nim_algo/EMThreads.cpp
        src/em_nim_algo/EMThreads.h
        src/em_nim_algo/likelihood_functions.cpp
        src/em_nim_algo/likelihood_functions.h
        src/em_nim_algo/DEArgs.cpp
        src/em_nim_algo/DEArgs.h
        src/em_nim_algo/DEResult.cpp
        src/em_nim_algo/DEResult.h
        src/em_nim_algo/NIMRunner.cpp
        src/em_nim_algo/NIMRunner.h
        src/em_nim_algo/NIMThreads.cpp
        src/em_nim_algo/NIMThreads.h
        src/em_nim_algo/adaptiveIntegrate.h
        src/em_nim_algo/adaptiveIntegrate.cpp
        src/misc/basic_math.cpp
        src/misc/basic_math.h
        src/misc/CLIArgs.cpp
        src/misc/CLIArgs.h
        src/misc/command_line_parser.cpp
        src/misc/command_line_parser.h
        src/misc/ParsedCountsY.cpp
        src/misc/ParsedCountsY.h
        src/misc/ParsedCovariatesX.cpp
        src/misc/ParsedCovariatesX.h
        src/misc/ParsedData.cpp
        src/misc/ParsedData.h
        src/misc/ParsedERCC.cpp
        src/misc/ParsedERCC.h
        src/misc/utils.cpp
        src/misc/utils.h
        src/param_estimation/init_integrand.cpp
        src/param_estimation/init_integrand.h
        src/param_estimation/InitMargLklhoodArgs.cpp
        src/param_estimation/InitMargLklhoodArgs.h
        src/param_estimation/IntegrationArgs.cpp
        src/param_estimation/IntegrationArgs.h
        src/param_estimation/logistic_regression_irls.h
        src/param_estimation/logistic_regression_irls.cpp
        src/bio_var_computer/BioVarComputer.cpp
        src/bio_var_computer/BioVarComputer.h
        src/bio_var_computer/BioVarRunner.cpp
        src/bio_var_computer/BioVarRunner.h
        src/bio_var_computer/NIMResult.cpp
        src/bio_var_computer/NIMResult.h
        src/bio_var_computer/NimMarginalLikelihood.cpp
        src/bio_var_computer/NimMarginalLikelihood.h
        src/bio_var_computer/NimArgs.cpp
        src/bio_var_computer/NimArgs.h
        src/TASC.cpp)


add_executable(TASC ${SOURCE_FILES})
target_link_libraries(TASC ${LIBS})
