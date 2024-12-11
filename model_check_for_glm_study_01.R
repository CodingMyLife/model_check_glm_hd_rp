## Model checking for high dimensional generalized linear models based on random projections
# nohup Rscript ./main/model_check_for_glm_study_01.R > ./tmp/model_check_for_glm_study_01.nohup 2>&1 &
rm(list = ls())
library(MASS)
library(glmnet)
library(randomForest)
library(log4r)
library(tidyr)
library(dplyr)
library(psych)
library(pracma)
library(foreach)
library(doParallel)
library(harmonicmeanp)
library(PLStests)
# library(GRPtests)
# library(RPtests)

## 0 . this function for generating data expreesed in section 5.1 study 1
gen.data <- function(n, p, a,rho,func){
  
  beta_0 = rep(0,p)
  beta_0[1:5] = rep(1,5)/sqrt(5)
  beta_0 <-  matrix(beta_0,ncol = 1)   # parameters
  
  p1 = floor(2*(n^(1/3)))
  p1=ifelse(p1>p,p,p1)
  beta_1 <- rep(0,p)
  beta_1[1:p1]=1/sqrt(p1)
  beta_1 <-  matrix(beta_1,ncol = 1)
  
  x <- rho^c(0:(p-1))
  sigma0 <- toeplitz(x)
  
  X <- mvrnorm(n, rep(0,p), sigma0)
  
  if (func=="H11"){
    y = X%*%beta_0 + a*exp(-(X%*%beta_0)^2)+matrix(rnorm(n),ncol = 1)
  }
  if (func=="H12"){
    y = X%*%beta_0 + a*cos(0.6*pi*(X%*%beta_0))+matrix(rnorm(n),ncol = 1)
  }
  if (func=="H13"){
    y = X%*%beta_0 + a*(X%*%beta_1)^2+matrix(rnorm(n),ncol = 1)
  }
  if (func=="H14"){
    y = X%*%beta_0 + a*exp(X%*%beta_1)+matrix(rnorm(n),ncol = 1)
  }
  
  return(list(X = X, y = y))
}

hh = c("H11","H11","H13","H14")

##
## in our paper, there are four model to generate data, the looping variable hi can be set to 1:4.
## 

for (hi in 1:1) {
  
  pre_file_name = paste("model_check_for_glm_study_01",hh[hi],"_",
                        format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),sep="_")
  # logging file
  dir = getwd()
  log_file = file.path(dir,"tmp",paste(pre_file_name,".log",sep = ""))
  logger <- create.logger(logfile = log_file, level = "INFO")
  
  nexp = 1000  ## number of experiment per parameters.
  n_alpha = 10 ## the number of project direction
  num_to_count_time = nexp/2 ## per number of experiment for logging
  #set the number of cores need for your project. do not greater than totalCores
  totalCores = detectCores()
  number_of_cores = totalCores-2  # number of cores for parallel
  h = c(2)
  nexps = c(1:nexp)
  length.a = 5
  
  #### paramters for data generating n=200
  n=c(200)    ## number of samples in one test.
  func <- c("H11")     #
  rhos<-c(0,0.5)       # correlation of x
  a_1 <- seq(0,0.8,length.out=length.a)
  p <-  c(10,100,200,400)

  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H12")      
  a_2 <- seq(0,2,length.out=length.a)
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H13")      
  a_3 <- seq(0,0.4,length.out=length.a)
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H14")      
  a_4 <- seq(0,0.2,length.out=length.a)
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  exp_para_exp_02 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)
  
  ####  paramters for data generating n=800 p=500
  n=800    
  p=500
  
  func <- c("H11")
  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H12")      
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H13")      
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H14")      
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  exp_para_exp_03_1 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)
  
  
  #### paramters for data generating n=2000 p=3000
  n=2000    
  p=3000
  
  func <- c("H11")
  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H12")      
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H13")      
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  func <- c("H14")      
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)
  
  exp_para_exp_03_2 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)
  
  
  ## all paras..
  exp_para = rbind(exp_para_exp_02,exp_para_exp_03_1,exp_para_exp_03_2)
  colnames(exp_para)=c("n","func","p","a","rhos","nexps")
  exp_para = exp_para %>% arrange(desc(p))
  
  ####  each row of exp_para be the parameters of one simulation 
  ####  set n==200,  data generating function as "H11" if hi=1, a==0.
  ####  so p and rhos are all values defined before.
  ####  you can set n=800 or n=2000 to get other simulation result.
  exp_para = exp_para %>% filter(
    func==hh[hi],
    n==200
    ,a==0
    # ,rhos==0
    # ,p==400
  ) %>% arrange(rhos,desc(p),desc(n))
  
  start.time <- Sys.time()
  info(logger,paste("start time at ",start.time))

  result_file = file.path(dir,"tmp",paste(pre_file_name,".cluster",sep = ""))
  cluster <- makeCluster(number_of_cores,outfile=result_file)
  registerDoParallel(cluster)
  
  number_of_experiment_type = nrow(exp_para)
  number_e = 1:number_of_experiment_type
  
  result_all <- foreach(an = number_e,.combine=rbind,
                        .packages=c("MASS","glmnet","randomForest","log4r",
                                    "GRPtests","dplyr","pracma","tidyr","psych",
                                    "harmonicmeanp","PLStests")) %dopar% {
    i = an
    start.time <- Sys.time()
    params = exp_para[i,]

    fam  = "gaussian"
    func_ij = as.character(params$func)
    gd <- gen.data(n=params$n[1], p=params$p[1], a=params$a[1],rho=params$rhos[1],func=func_ij)
    x <- gd$X
    y <- gd$y
    
    # grp_test = tryCatch(GRPtest(x,
    #                             y,
    #                             fam = fam,
    #                             nsplits = 1,
    #                             output_all = TRUE),
    #                     error = function(e) return(list(pval=NA)))
    
    # rp_test = tryCatch(RPtest(x,
    #                           y,
    #                           test = "nonlin",
    #                           resid_type = "Lasso",
    #                           output_all = TRUE),
    #                    error = function(e) return(list("p-value"=NA))
    # )
    
    proj_test = PLStests(y, x,family = fam, b0=h, np =n_alpha)
    proj_test = as.data.frame(proj_test)
    num_h = length(h)
    result_temp <- params[rep(1,num_h),-ncol(params)]
    result_temp <- cbind(result_temp,proj_test)
    
    result_temp = result_temp %>% gather("test_stat" ,p_value,7:ncol(result_temp),factor_key =TRUE)
    result_all_a_param = result_temp
    
    
    ## logging size or power of this paras
    if (an%%num_to_count_time == 0) {
      info(logger,paste("round i=",an))
      
      infos = ""
      d_parms = dim(params)[2]
      for (i1 in 1:d_parms) {
        infos= paste(infos,colnames(params)[i1],":",params[1,i1]," ")
      }
      info(logger ,infos)
      
      end.time = Sys.time()
      time.taken = end.time - start.time
      time.taken.mins = as.numeric(time.taken, units = "mins")
      
      have_runed_turn = an
      info(logger,paste("time elspse:",time.taken.mins,"mins, and time left:",
                        time.taken.mins*(number_of_experiment_type-have_runed_turn)/number_of_cores,"minus"))
      # start.time <- Sys.time()
    }
    result_all_a_param
  }
  
  ## cluster over
  stopCluster(cluster)
  
  result_p_file = file.path(dir,"result",paste(pre_file_name,
                                               format(Sys.time(), "to_%m_%d_%H_%M_%S"),".csv",sep = "_"))
  result_all$greater005 <- result_all$p_value < 0.05
  write.csv(result_all,file =result_p_file,row.names = FALSE)
  
  result_all$h <- as.factor(result_all$h)
  
  result_all_agg<-aggregate(result_all$greater005, by=list(result_all$n,result_all$func, result_all$p,
                                                           result_all$a,result_all$rhos,result_all$h,result_all$test_stat), FUN=mean,na.rm=TRUE)
  colnames(result_all_agg)=c("n","func","p_dim","a","rhos","h","test_stat","size_power")
  result_all_agg = result_all_agg %>%
    mutate(succ_less007= ifelse(size_power<0.071,1,0),fail_bigger01 = ifelse(size_power>0.091,1,0))
  result_agg_file = file.path(dir,"result",paste(pre_file_name,"agg",format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),".csv",sep = "_"))
  write.csv(result_all_agg,file =result_agg_file,row.names = FALSE)
}
