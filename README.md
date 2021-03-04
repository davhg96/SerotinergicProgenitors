# SerotinergicProgenitors
# Reproducibility
To reproduce this analysis first create an R project in your local machine and pull the files from this GitHub repository. Reproducibility in this project is managed by the package renv. To start working on this firs install the package:
```
if (!requireNamespace("remotes"))
  install.packages("remotes")

remotes::install_github("rstudio/renv")
```
At this point you should have all necessary files downloaded in a directory in your local machine, check that you can see a file called `rent.lock`. Set the working directory to that folder and run `renv::restore()` to install all the necessary packages and versions to reproduce this analysis.

# Analisis files

All these scripts will access the data from a data directory located in your working directory and will create an output directory with all the output automatically.

- Functions.R Will load all the necessary packages and custom functions necessary to reproduce the analysis
- Analysis.R Contains the code for the differential expression analysis.
- Plotting.R Contains the code to reproduce all the images.