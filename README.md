# Model checking for high dimensional generalized linear models based on random projections

## Contents of this repository
This repository contains:
1. R package PLStests
2. Code reproducing empirical results from Section 5.1 of paper

## Installation of R package PLStests

Installation in R from this github repository:

```
install.packages(devtools)
library(devtools)
install_github("CodingMyLife/model_check_glm_hd_rp/PLStests")
```
or download the project on your directory , and then using 
```
devtools::install_local("your_dir/model_check_glm_hd_rp/PLStests")
```
## Simulation of section 5.1 

### model_check_for_glm_study_01.R
- step 0: download this project, and unzip to .../yourhomedir/model_check_glm_hd_rp, where the model_check_for_glm_study_01.R and other files included in this directory.
- step 1: change your rstadio working directory to .../yourhomedir/model_check_glm_hd_rp
- step 2: source("./model_check_glm_hd_rp/model_check_for_glm_study_01.R") [or run in rstadio]
- step 3: the size or power  results be writen in fold "result", while runtime loggings in "tmp"

The outputs of this code will be saved in the the folds tmp and result. As the name indicate, logging files lay in tmp and the p value of our statistics lay in result fold. The size and power record in model_check_for_glm_study_01_xxx_agg_xxx.csv. In this table, you will find size or power of different combinations of n,p,rho,and a of different models .

It will take 10 minums to run this code. we comment out the GRP and RP test for time saving. If you want it, you should add it easily following our framwork. more details can be found in the comments.

### model_check_for_glm_study_02.R

same as before. more detail can be found in the comments in the file.
