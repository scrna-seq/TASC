# TASC
Toolkit for Analysis of Single Cell data

### INSTALLATION
#### Before installation:
######Make sure you have the following libraries and packages installed on your system.
1. Boost C++ Library >= 1.35.0
2. GNU Scientific Library >= 1.15
3. cmake >= 2.6.4
4. gcc >= 4.4

The program has not been extensively tested on older versions of the above packages and libraries. Compilation and runtime error may occur if older versions are used. If you happen to work on a cluster running an outdated distro, please contact me for further support. But keep in mind I can only support reasonably old legacy distributions, for our program uses functions found only in more recent (as in post-2011) versions of certain libraries. 

#### Installation steps:
1. Download and unzip the project to a folder, for example, TASC-master/
2. Change the current directory to the root of the project: ```cd TASC-master```
3. Create a build folder: ```mkdir build```
4. Change the current directory to the build folder: ```cd build```
5. Prepare compilation configuration with cmake: ```cmake ../``` If the version of your boost library is old, use this instead: ```cmake -DBoost_NO_BOOST_CMAKE=ON ../``` If you wish to provide your own boost root, use this option ```-DBOOST_ROOT=/path/to/boost/root/``` Sometimes the compilation might fail if you don't provide a boost library, in this case, use this option ```-DBOOST_LIBRARYDIR=/path/to/boost/lib/``` 
6. Compile: ```make```
7. If the compilation is successful, the executable will appear in the build folder. If it fails, please check the versions of your Boost library and Gnu Scientific Library. Also make sure you only have ONE copy of the boost library. For example, if you also have rsem on your system, cmake might try to use the boost library embedded in rsem. The compilation will fail if that happens.

#### Binary:
If you have trouble compiling the program yourself, you can try the statically linked executable [here](https://raw.githubusercontent.com/scrna-seq/TASC/master/binary/TASC).
I have tested the executable on Red Hat Enterprise Linux Server release 6.6 (Santiago) and Ubuntu Linux 15.10 (Wily Werewolf).

### USAGE
In order to perform the DE gene analysis, TASC requires 3 headerless ASCII coded files, x.txt, y.txt and ercc.txt.

#### x.txt
This file codes the design matrix for your model. If you wish to test for differential expression between two groups of 4 cells, then the following x.txt can be used.

```
intercept   1,1,1,1,1,1,1,1
group       0,0,0,0,1,1,1,1
```

x.txt has no header. The first column contains the name of the covariates, and the second column contains the value of the covariates separated by comma. The two columns are separated by one tab.

If you have three groups of 3 cells, then the following x.txt can be used.

```
intercept   1,1,1,1,1,1,1,1,1
group_2     0,0,0,1,1,1,0,0,0
group_3     0,0,0,0,0,0,1,1,1
```

At this point, TASC can only test for the significance of a single covariate. We have yet to implement the testing procedure for contrasts.

#### y.txt
This file contains the read counts of each cell for all the genes. If you are testing for the differential expression of two genes G1 and G2, between two groups of 4 cells, y.txt should be formattted as, 

```
G1   79,5,91,67,42,88,0,63
G2   54,100,53,38,85,71,37,22
```

Again, y.txt has no header, so make sure the order of the columns correspond to the order in x.txt and ercc.txt. The read count for each cell is separated by commas, and the two columns are separated by a tab.

Note: y.txt does not allow entries with identical gene names. Please modify your file if different genes with duplicate names exist in your input.

#### ercc.txt
In order to provide the information on ERCC spike-ins, TASC depends on a separate file containing the requisite information. For example, your input has 4 spiked-in molecules (usually the number is around 60-90) with the name ERCC_1, ERCC_2, ERCC_3, ERCC_4. We need to provide the information for these spike-ins using the following format.

```
ERCC_1   1.5   1000   86,37,99,69,29,16,5,92
ERCC_2   2.9   1031   56,47,99,94,38,5,21,23
ERCC_3   3.1    905   67,98,10,78,5,91,27,54
ERCC_4   6.2    865   61,14,92,20,32,36,0,47
```

ercc.txt has four columns, each containing some important information w.r.t the spiked-in molecules. The first column is the name of the spike-in. The second column contains the true concentration of the ERCC molecules in the final sequencing mix. The third column contains the length of the specific molecules. The last column contains the read counts for all the cells in the order of x.txt and y.txt, separated by commas. The columns are separated by tabs. Note, it is preferred that the concentration of the ERCC molecules be measured in number of molecules per reaction. Other units could work in theory as well. For now TASC assumes the same amount of ERCC molecules are added to every reaction. Our model actually allows for more flexible definition of the concentrations, so in later updates we will relax this assumption in our software.

#### running TASC
After getting the required information ready, we can run TASC by invoking its command like this:

```
TASC -b -y y.txt -x x.txt -k ercc.txt -o output_dir/output_basename
```

You can try out TASC using the x.txt, y.txt and ercc.txt under the example/ folder in our github repo. 

Note: our method uses the Nelder-Mead simplex algorithm to numerically optimize the marginal likelihood. This requires a starting point on which the likelihood can be evaluated. Numerical underflow might happen if this starting point is not properly chosen. In this version, we set the starting point for all the coefficients in the linear model to be 1, which should be suitable for the majority of our users. However, if you are experiencing numerical instability, please [contact me](mailto:upenn-scrna-seq@outlook.com) for further assistance.

#### parallel acceleration
One of the important features of TASC is that it leverages thread-level parallelism to significantly speed up the computation. If you wish to use more than 1 core (which is highly recommended) to perform your analysis, please indicate that to TASC with the following option `-n`. For example, if you wish to use 8 threads to run the previous task, 

```
TASC -n 8 -b -y y.txt -x x.txt -k ercc.txt -o output_dir/output_basename
```

By default, TASC will use only 1 thread.

#### perform quantification
Under certain circumstances, one wishes to quantify the gene expression instead of performing DE analysis. You need to provide the ercc.txt and y.txt as mentioned above. In addition, you still have to provide a dummy x.txt in this case. The following x.txt is needed if you are quantifying the expression of the genes in 8 cells.

```
intercept   1,1,1,1,1,1,1,1
```

After these files are ready, use the `-a` option of TASC to perform quantification.
```
TASC -a -t 1 -y y.txt -x x.txt -k ercc.txt -o output_dir/output_basename
```
Notice that I have manually set the column for testing to be 1, this is because by default TASC tests for the second predictor (usually the experimental group if you are doing two sample comparison). This is not possible when the dummy design matrix x only contains 1 covariate. Switch to ``-t 1`` if you are only performing quantification with a dummy design matrix. If you are doing both, you should use a real design matrix x, and set ``-t`` to be whatever columns you wish to test for.

More on testing for columns, the ``-t`` option can take more than 1 covariates, for example, if you wish to test for the significance of both covariate No. 1 and No. 2, you can use the option like this: ``-t 1,2``

#### Laplacian approximation
The adaptive integration method we have implemented in TASC is extremely accurate yet slow, especially if you are not accessible to some 20-core machine in a cluster. In this case, please try the `-f` option, which forces TASC to use Laplacian approximation instead of adaptive qudrature to evaluate the integrals. In our experience, this greatly speeds up the program without sacrificing too much on accuracy. By default TASC uses the adaptive integration approach.

### Mac Users
TASC code CAN run on Mac OS X, but compiling TASC on Mac IS complicated. TASC is written with gcc/g++. However, since OS X Mavericks, Apple has abandoned gcc/g++ in their operating system, and will no longer provide libstdc++ with Xcode. If you want to compile TASC under Mac, you will be required to:

1. download and install gcc/g++ (homebrew). Remember to use the --without-multilib option, otherwise it won't support openMP. Don't be alamred if it runs for a while, homebrew is actually building gcc from source.
2. download and install boost library. Homebrew version does not work, as it is compiled with clang. You need to download the source of boost and compile it yourself using the gcc installed in the previous step following this [link](https://solarianprogrammer.com/2016/03/06/compiling-boost-gcc-5-clang-mac-os-x/).)
3. download and install gsl (homebrew)
4. download and install git and cmake (homebrew)
5. clone the TASC project and modify the CMakeLists.txt and change line 3 to 
  ```
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -Wall -std=c++0x -L/usr/local/boost-1.60.0/lib -lboost_program_options")
  ```
  The ``/usr/local/boost-1.60.0/lib`` should be replaced with the directory in which your boost library files are installed to.
6. build with the brewed g++. On my system, the cmake command is 
  ```
  cmake -DBOOST_ROOT=/usr/local/boost-1.60.0/ -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-5 -DCMAKE_C_COMPILER=/usr/local/bin/gcc-5 ../
  ```

clang can also compile TASC. However, the default clang on OS X does not support openMP, which means that you have to install clang-omp from homebrew, and change the source code of TASC to include the proper omp.h header and provide the location of the libiomp. Also you might need to change the c++ standard option in the CMakeFiles.txt. I have not tested TASC using clang, so unfortunately you are on your own. 

The above steps would require some understanding of the compiling and linking process for C/C++. Let me know if you are running into any trouble. 

I apologize for not being able to provide a version of statically linked binary for OS X. This is partially due to [Apple's strong disapproval of statically linked executables](https://developer.apple.com/library/mac/qa/qa1118/_index.html). Another reason is that Apple has abandoned gcc/g++ since Mavericks, and their default C++ compiler does not natively support openMP, the de-facto thread-level parallelization library required by TASC. The workarounds are tedious and not necessarily beginner-friendly. 

### Technical Support
Technical support is provided by [me](mailto:upenn-scrna-seq@outlook.com). Shoot me an E-mail if you have any question regarding our software.

### License
Copyright (c) <2016>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
