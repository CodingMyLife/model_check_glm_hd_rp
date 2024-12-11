## Model checking for high dimensional generalized linear models based on random projections
# nohup Rscript ./main/main_god_linear_exp_n200_package.R > ./tmp/main_god_linear_exp_n200_package.nohup 2>&1 &
rm(list = ls())
library(MASS)
library(glmnet)
library(randomForest)
library(log4r)
library(GRPtests)
library(RPtests)
library(tidyr)
library(dplyr)
library(psych)
library(pracma)
library(foreach)
library(doParallel)
library(harmonicmeanp)
library(PLStests)


## 0 . this function for generating data expreesed in section 5.1 study 1
gen.data <- function(n, p, a,rho,func){

  beta_0 = rep(0,p)
  beta_0[1:5] = rep(1,5)/sqrt(5)
  beta_0 <-  matrix(beta_0,ncol = 1)   # parameters

  p1 = floor(2*(n^(1/3)))
  p1=ifelse(p1>p,p,p1)
  # p1 = floor(n/10)-5
  beta_1 <- rep(0,p)
  beta_1[1:p1]=1/sqrt(p1)
  beta_1 <-  matrix(beta_1,ncol = 1)

  # beta_2 <- rep(0,p)
  # beta_2[(p-p1+1):p]=1/sqrt(p1)
  # beta_2 <-  matrix(beta_2,ncol = 1)

  x <- rho^c(0:(p-1))
  sigma0 <- toeplitz(x)

  X <- mvrnorm(n, rep(0,p), sigma0)*1#matrix(runif(n*p),n,p)#

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
# hh = c("H13")
for (hi in 1:1) {

  set.seed(seeds[iv]+1000*hi)
  pre_file_name = paste("main_linear_all_n200_b02",hh[hi],"_",format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),sep="_")

  # logging file
  dir = getwd()
  log_file = file.path(dir,"tmp",paste(pre_file_name,".log",sep = ""))
  logger <- create.logger(logfile = log_file, level = "INFO")



  nexp = 1000   ## number of experiment per parameters.
  n_alpha = 10 ## the number of project direction
  num_to_count_time = nexp/2 ## number of experiment for count time.
  number_of_cores = 4  # number of cores for parallel
  h = c(2)
  nexps = c(1:nexp)
  length.a = 5

  #### paramters for data generating n=200
  n=c(200)    ## number of samples in one test.
  func <- c("H11")      #magitute of mis
  rhos<-c(0,0.5)       # correlation of x
  a_1 <- seq(0,0.8,length.out=length.a)
  p <-  c(10,100,200,400)
  # p <-  c(200,400)
  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H12")      #magitute of mis
  a_2 <- seq(0,2,length.out=length.a)
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)


  func <- c("H13")      #magitute of mis
  a_3 <- seq(0,0.4,length.out=length.a)
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H14")      #magitute of mis
  a_4 <- seq(0,0.2,length.out=length.a)
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  exp_para_exp_02 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)

  ####  paramters for data generating n=800 p=500
  n=800    ## number of samples in one test.
  p=500

  func <- c("H11")
  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H12")      #magitute of mis
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H13")      #magitute of mis
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H14")      #magitute of mis
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  exp_para_exp_03_1 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)


  #### paramters for data generating n=2000 p=3000
  n=2000    ## number of samples in one test.
  p=3000

  func <- c("H11")
  exp_para11 = expand.grid(n,func,p,a_1,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H12")      #magitute of mis
  exp_para12 = expand.grid(n,func,p,a_2,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H13")      #magitute of mis
  exp_para13 = expand.grid(n,func,p,a_3,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  func <- c("H14")      #magitute of mis
  exp_para14 = expand.grid(n,func,p,a_4,rhos,nexps,KEEP.OUT.ATTRS = FALSE)

  exp_para_exp_03_2 = rbind(exp_para11,exp_para12,exp_para13,exp_para14)


  ## all paras..
  exp_para = rbind(exp_para_exp_02,exp_para_exp_03_1,exp_para_exp_03_2)
  colnames(exp_para)=c("n","func","p","a","rhos","nexps")
  exp_para = exp_para %>% arrange(desc(p))

  #### set n=200, and data generating function as
  exp_para = exp_para %>% filter(
    func==hh[hi],
    n==200
    # ,a==0
    # ,rhos==0
    # ,p==400
  ) %>% arrange(rhos,desc(p),desc(n))

  test_stat <- c("T_cauchy","T_alpha","T_hmp","T_beta")

  start.time <- Sys.time()
  info(logger,paste("start time at ",start.time))

  # totalCores = detectCores()
  #set the number of cores need for your project. do not greater than totalCores
  result_file = file.path(dir,"tmp",paste(pre_file_name,".cluster",sep = ""))
  cluster <- makeCluster(number_of_cores,outfile=result_file)
  registerDoParallel(cluster)

  number_of_experiment_type = nrow(exp_para)
  an = 1:number_of_experiment_type

  result_all <- foreach(an = an,.combine=rbind,
                        .packages=c("MASS","glmnet","randomForest","log4r",
                                    "GRPtests","dplyr","pracma","tidyr","psych","harmonicmeanp")) %dopar% {
                                      i = an
                                      start.time <- Sys.time()
                                      result_all_a_param <- data.frame(n=numeric(0),
                                                                       func = factor(NULL,levels =func),
                                                                       p = logical(0),
                                                                       a = numeric(0),
                                                                       rhos = numeric(0),
                                                                       h = numeric(0),
                                                                       test_stat = factor(NULL,levels =test_stat ),
                                                                       p_value = numeric(0))

                                      params = exp_para[i,]

                                      len_test_stat = length(test_stat)
                                      for (j in 1:1) {
                                        # set.seed(j+2000+(as.numeric(params$func)-1)*1000)
                                        # set.seed(j+2000)
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
                                        grp_test_pval = 0

                                        # rp_test = tryCatch(RPtest(x,
                                        #                           y,
                                        #                           test = "nonlin",
                                        #                           resid_type = "Lasso",
                                        #                           output_all = TRUE),
                                        #                    error = function(e) return(list("p-value"=NA))
                                        # )
                                        rp_test_pval = 0

                                        # proj_test = god_proj_hd_multi_h_multiple(y,x,method.ls="SPARSE",family=fam,gived_proj = alpha,h=h)
                                        proj_test = PLStests(y, x,family = fam, b0=h, np =n_alpha)
                                        # LSTRP_grp_cauchy_double_tryCatch(y, x, b0=h, np =n_alpha)
                                        num_h = length(h)
                                        result_temp <- params[rep(1,num_h),-ncol(params)]
                                        result_temp <- cbind(result_temp,proj_test)
                                        result_temp = result_temp %>% gather("test_stat" ,p_value,7:10,factor_key =TRUE)
                                        result_all_a_param = rbind(result_all_a_param,result_temp)
                                      }

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
