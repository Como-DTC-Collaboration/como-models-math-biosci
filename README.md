# Learning transmission dynamics modelling of COVID-19 using ComoModels

Interactive notebooks accompanying the article _Learning transmission dynamics modelling of COVID-19 using ComoModels_.

## Contents

This repository holds three R Markdown vignettes:

1. [S1: On the numerical solution of compartmental models](numerical_solution/numerical_solution_of_SEIRD.Rmd)
2. [S2: The importance of uncertainty in age-specific contact patterns for quantifying COVID-19 risk](contact_matrix/SEIRD_contact_matrix.Rmd)
3. [S3: The importance of sensitivity analysis and model inference](optimisation/SEIRD_optimisation.Rmd)

Within each folder, there are also rendered pdf files.

## Instructions

All vignettes require that the `comomodels` package be installed, which can be executed via the following R command:

```R
# install devtools if necessary through
# install.packages("devtools")

# install package
devtools::install_github("Como-DTC-Collaboration/como-models")
```

