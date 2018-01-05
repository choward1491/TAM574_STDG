# Introduction
The purpose behind this GitHub project is to showcase materials related to the final project for a graduate course TAM 574, called **Advanced Finite Element Methods**, taught by Professor Robert Haber at UIUC.

# The Project
This GitHub project contains the source code, some post-processing scripts, some videos, and the final report for my submission regarding the final project.

The final project from the course essentially revolved around building a **Space-Time Discontinuous Galerkin** (STDG) Finite Element code for solving a small system of Hyperbolic Partial Differential Equations. 
When developing this code, one would aim to validate the code by finding the parabolic limit of the system of hyperbolic equations and comparing the numerical solution of this parabolic limit to the exact solution.
After the validation, the next step would be to gain an idea of the convergence of the STDG method by using some error estimators and plotting their values for various values of refinement and polynomial order within the elements.

### Comments
The report located in the **docs/** directory in the project helps to document and show the results of the analysis done for the above areas, while the source code helps show what I implemented at the time of the project to get the results within the report.

# Software and Visualization
To compile the software, just type `make` at the command line. After compilation takes place, you can then `cd` to the `bin/` directory and run the optimized executable `sim_exec` by typing `./sim_exec`. This executable needs to be run in a location where the files `GaussAbscissa.txt` and `GaussWeights.txt` exist, hence why you could start out by running it in the `bin/` directory. 

After running the simulation, some `x_*`,`u_*`,`q_*` data files will be generated in the same location as where the simulation was run. You can use/modify the matlab script in `matlab/scripts/plotSolution.m` to display the final solution for each field and the space-time solution for each field using these output data files.

Note that this software has not been touched since 2013 other than to add comments or rename a couple setup variables in `main.cpp`, so the software is not as flexible as it could be. If future work was to be done, it would be great to make it possible to pass in configuration files that can specify where quadrature data is, the file paths and names for data we want to output, the polynomial basis choice, and more. As of now, things are hard coded in the executable, so one would need to manually change things either in the code or post-execution of the software.
