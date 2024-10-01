# Microarray patch vaccines for typhoid conjugate vaccines: a global cost-effectiveness analysis
=======================

Administered by Marina Antillon, PhD.  

Collaborating centers:

Department of Epidemiology and Public Health   
Swiss Tropical and Public Health (Swiss TPH) Institute   
(An affiliated institute of the University of Basel)   
4123 Allschwil, Switzerland  

PATH
Seattle, WA, USA

License: GPL-3.0 

For questions, please contact Marina Antillon via LinkedIn. 
---

# Project Objective 

In the current study, we leverage a previously developed model of typhoid transmission (LINKS) to evaluate the impact and cost-effectiveness of a potential micro-array patch typhoid conjugate vaccine not yet under development. The analysis assumes that traditional (needle-and-syringe) typhoid conjugate vaccines would be rolled out before 2032 in most or all low- and middle-income countries and that MAPs would be non-inferior in terms of the probability and duration of protection. The analysis is also novel in taking into account the 

The analysis encompasses two major geographical scopes:
A. a global analysis where the cost-effectiveness of TCV-MAPS is considered on a national basis.
B. a sub-national analysis where the cost-effectiveness of TCV-MAPS is considered on a district-by-district basis.
---

# Brief description

The analysis is broadly defined in four parts, each of four parts is executed in various R code files:  
I. Re-fitting incidence parameters to accommodate for 5 wealth quintiles in all countries.
II. Simulating vaccine impact for national-level analysis
    + Simulating the impact of N&S vaccines from deployment in the 2020's through the 2030's and 2040's as a comparator.
    + Simulating the impact of MAPS vaccines in the 2030's.
    + Cost-effectiveness analysis
III. Simulating vaccine impact for sub-national (district) level analysis in five example countries
    + Simulating the impact of N&S vaccines from deployment in the 2020's through the 2030's and 2040's as a comparator.
    + Simulating the impact of MAPS vaccines in the 2030's.
    + Cost-effectiveness analysis
V. Collating the results
    + National-level analyses
    + Subnational analyses

## File structure, output, and dependencies

** TO WRITE **
---
# Computational considerations

## Hardware

**Hardware needs:** For optimal performance, the model was run for different places and different scenarios using a high-performance computing cluster at the University of Basel (scicore, http://scicore.unibas.ch/). Scicore is run with a slurm scheduler, and the bash file is included in this repository. However, a user could use a parallel computing package in R or any other automated task management solution, but no such implementation is presented here.

**Memory needed:** Intermediate simulations and graphs need about 10 GB of storage for all outputs.

## Software
The tools for the analysis are coded in R. While we would highly recommend using the code within RStudio environment (in part because of its features to manage the project with .RProj and packrat) this is not strictly necessary and the benefits of git and renv are available from a classic R interface and a shell command line.

Some of the results tables for the project are produced automatically with the code in the project. This is done by employing the knitr and kableExtra packages. 

A single run of the code takes about 60 minutes in a MacBook Air (M1 chip, 2020 model) with 16 GB of RAM. 

# Installation to-do list (all free)
- [R](https://www.r-project.org) (required)
- [renv](https://rstudio.github.io/renv/index.html) (packages installed for R; required)
- [RStudio](https://www.rstudio.com/) (highly recommended though not required; it integrates well with the other management tools used in the project. The free version is enough.)
- Latex engine (not required, though some documentation created automatically won't run). See [here](https://support.rstudio.com/hc/en-us/articles/200532257-Customizing-LaTeX-Options) for the recommended engines.

---
# How to run this analysis 

** TO WRITE **

---

# Troubleshooting

**Issues with the renv repository:** If there is an issue with the repository, try typing renv:: rinit() and when prompted type the number 1, for Restored the project from lockfile.

**Adding a package to the renv library:** just use the command install.packages(). Note that doing this won't install the package for use with other projects. Then update the lockfile by typing renv::snapshot().
