# Introduction
The purpose behind this GitHub project is to showcase materials related to the final project for a graduate course TAM 574, called **Advanced Finite Element Methods**, taught by Professor Robert Haber at UIUC.

# The Project
This GitHub project contains the source code, some post-processing scripts, some videos, and the final report for my submission regarding the final project.

The final project from the course essentially revolved around building a **Space-Time Discontinuous Galerkin** (STDG) Finite Element code for solving a small system of Hyperbolic Partial Differential Equations. 
When developing this code, one would aim to validate the code by finding the parabolic limit of the system of hyperbolic equations and comparing the numerical solution of this parabolic limit to the exact solution.
After the validation, the next step would be to gain an idea of the convergence of the STDG method by using some error estimators and plotting their values for various values of refinement and polynomial order within the elements.

### Comments
The report located in the **docs/** directory in the project helps to document and show the results of the analysis done for the above areas, while the source code helps show what I implemented at the time of the project to get the results within the report.
