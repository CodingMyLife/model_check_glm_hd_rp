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
or download the project, using 
```
devtools::install_local("your_dir/model_check_glm_hd_rp/PLStests")
```
## Simulation of section 5.1 

### model_check_for_glm_study_01.R
- step 0: download this project, and unzip to .../yourhomedir/model_check_glm_hd_rp-main
- step 1: change your r workdir to .../yourhomedir/model_check_glm_hd_rp-main
- step 2: create  two folds "tmp" and "result" in the fold model_check_glm_hd_rp
- step 3: source("./model_check_glm_hd_rp/model_check_for_glm_study_01.R") [run in rstadio]

the outputs of this code will save in the the fold tmp and result. as the name indicate, logging file lay in tmp and the p value of our statistic lay in fold result. The size and power record in model_check_for_glm_study_01_agg_xxx.csv. in the table, you will find size or power of different combinations of n,p,rho,and a of different models .

It will take 10 minums to run this code. we comment out the GRP and RP test for time saving. If you want it, you should add it easily following our framwork. more details can be found in the comments.

### model_check_for_glm_study_02.R

same as before. more detail can be found in the comments in the file.
