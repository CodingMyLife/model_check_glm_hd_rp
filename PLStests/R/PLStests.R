#' @import glmnet
#' @import harmonicmeanp
#' @import MASS
#' @import psych
#'
PLStest_GLM <- function(y, x, b0=2, np=10) {

  n = nrow(x)
  p = ncol(x)
  b0_tmp = b0
  num_h = length(b0_tmp)

  LST_pval_cauchy <- rep(NA,num_h)
  LST_pval_single <- rep(NA,num_h)
  LST_pval_hmp <- rep(NA,num_h)
  LST_pval_beta_hat <- rep(NA,num_h)

  #lasso
  mod_cv0 <- glmnet::cv.glmnet(x, y, family = "binomial", maxit=100000,intercept = F)
  lamuda <- c("lambda.min", "lambda.1se")[1]
  beta_n <- coef(mod_cv0, s = lamuda)[-1]
  beta_none0 <- seq(1:ncol(x))[as.numeric(beta_n) != 0]

  for (id_h in 1:num_h) {
    b0 = b0_tmp[id_h]
    tryCatch({

      if (length(beta_none0) == 0) {

        pred.lasso = predict(mod_cv0, newx = x,type="response",s=mod_cv0$lambda.min)
        pred.lasso = matrix(unname(pred.lasso),ncol = 1)
        #residual
        U <- y
        #residual matrix ，ui*uj
        errormat <- U %*% t(U)
        #bandwidth
        h <- b0 * (n ^ (-1 / 4))
        #projections
        #new mu and sigma to distinguish the parameters that generate x
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)
        # bind beta_hat
        beta_hat =rnorm(p)
        beta_pro = rbind(beta_hat,beta_pro)

        for (k in 1:(np+1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }
        #p-value matrix
        pval_matrix <- matrix(nrow = np+1, ncol = 1)
        for (q in 1:(np+1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        pval_beta_hat = pval_matrix[1,1]

        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)
        #length(beta_none0)=1
      } else if (length(beta_none0) == 1) {

        X_S <- x[, beta_none0]
        fit.model = glm(y~X_S-1,family = binomial())
        pred = predict(fit.model, newx = X_S,type="response")
        pred = matrix(unname(pred),ncol = 1)

        beta_hat = unname(fit.model$coefficients)
        #residual
        # U = y - X_S * beta_hat
        U = y-pred
        #new x
        X_SC <- x[, -beta_none0]
        x_new <- cbind(X_S, X_SC)
        x <- x_new
        beta_hat <- c(beta_hat, rep(0, (p - length(beta_none0))))
        #residual
        #U=y-x%*%beta_n
        #residual matrix ，ui*uj
        errormat <- U %*% t(U)
        #bandwidth
        h <- b0 * (n ^ (-1 / 5))
        # projections
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)
        #beta1 <- beta_n
        # beta1 <- beta_hat
        beta_pro <- rbind(beta_hat, beta_pro)

        # bind beta_n
        # beta_pro = rbind(beta_pro,matrix(beta_n,nrow = 1))

        for (k in 1:(np + 1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }
        #p-value matrix
        pval_matrix <- matrix(nrow = np + 1, ncol = 1)
        for (q in 1:(np + 1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        pval_beta_hat = pval_matrix[1,1]
        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)
        #length(beta_none0)>1
      } else{
        X_S <- x[, beta_none0]
        fit.model = glm(y~X_S-1,family = binomial())
        pred = predict(fit.model, newx = X_S,type="response")
        pred = matrix(unname(pred),ncol = 1)
        beta_hat = unname(fit.model$coefficients)
        U = y-pred
        #new x
        X_SC <- x[, -beta_none0]
        x_new <- cbind(X_S, X_SC)
        x <- x_new
        beta_hat <- c(beta_hat, rep(0, (p - length(beta_none0))))
        #residual
        #U=y-x%*%beta_n
        #residual matrix ，ui*uj
        errormat <- U %*% t(U)
        #bandwidth
        len <- length(beta_none0)
        h <- b0 * (n ^ (-1 / (4 + len)))
        # 100 projections
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)

        beta_pro <- rbind(beta_hat, beta_pro)

        # # bind beta_n
        # beta_pro = rbind(beta_pro,matrix(beta_n,nrow = 1))

        for (k in 1:(np + 1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }
        #p-value matrix
        pval_matrix <- matrix(nrow = np + 1, ncol = 1)
        for (q in 1:(np + 1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        pval_beta_hat = pval_matrix[1,1]

        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)

      }
      #power
      result2 <- c(pval_cauchy, pval_single,pval_hmp,pval_beta_hat)
      LST_pval_cauchy[id_h] = result2[1]
      LST_pval_single[id_h] = result2[2]
      LST_pval_hmp[id_h] = result2[3]
      LST_pval_beta_hat[id_h] = result2[4]
      #error situation
    }, error = function(e) {
      #power
      print(e)
    })
  }

  #result
  return(
    data.frame(
      h = b0_tmp,
      T_cauchy = LST_pval_cauchy,
      T_alpha = LST_pval_single,
      T_hmp = LST_pval_hmp,
      T_beta = LST_pval_beta_hat
    )
  )
}

PLStest_LM <- function(y, x, b0=2, np=10) {

  n = nrow(x)
  p = ncol(x)
  b0_tmp = b0
  num_h = length(b0_tmp)
  LST_pval_cauchy <- rep(NA,num_h)
  LST_pval_single <- rep(NA,num_h)
  LST_pval_hmp <- rep(NA,num_h)
  LST_pval_beta_hat <- rep(NA,num_h)

  #lasso
  mod_cv0 <- glmnet::cv.glmnet(x, y, family = "gaussian",intercept = F)
  lamuda <- c("lambda.min", "lambda.1se")[1]
  beta_n <- coef(mod_cv0, s = lamuda)[-1]
  beta_none0 <- seq(1:ncol(x))[as.numeric(beta_n) != 0]
  for (id_h in 1:num_h) {
    b0 = b0_tmp[id_h]
    tryCatch({

      #length(beta_none0)=0
      if (length(beta_none0) == 0) {
        #residual
        U <- y
        errormat <- U %*% t(U)
        #bandwidth
        h <- b0 * (n ^ (-1 / 4))
        #projections
        #new mu and sigma to distinguish the parameters that generate x
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)
        # beta_hat is zeor, so choose randomly
        beta_hat = rnorm(p)
        beta_pro = rbind(beta_hat,beta_pro)
        for (k in 1:(np+1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }

        #p-value matrix
        pval_matrix <- matrix(nrow = np+1, ncol = 1) # 1 for beta_n
        for (q in 1:(np+1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        pval_matrix = matrix(pval_matrix,ncol = 1)
        # beta_hat
        pval_beta_hat = pval_matrix[1,1]
        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)
        #length(beta_none0)=1
      } else if (length(beta_none0) == 1) {
        X_S <- x[, beta_none0]
        beta_hat <- solve(t(X_S) %*% X_S) * (t(X_S) %*% y)
        #residual
        U = y - X_S %*% beta_hat
        #new x
        X_SC <- x[, -beta_none0]
        x_new <- cbind(X_S, X_SC)
        x <- x_new
        beta_hat <- c(beta_hat, rep(0, (p - length(beta_none0))))
        #residual
        #U=y-x%*%beta_n
        #residual matrix ，ui*uj
        errormat <- U %*% t(U)
        #bandwidth
        h <- b0 * (n ^ (-1 / 5))
        # 100 projections
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)
        #beta1 <- beta_n
        beta1 <- beta_hat
        beta_pro <- rbind(beta1, beta_pro)

        for (k in 1:(np + 1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }
        #p-value matrix
        pval_matrix <- matrix(nrow = np + 1, ncol = 1)
        for (q in 1:(np + 1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        pval_matrix = matrix(pval_matrix,ncol = 1) # remove the beta_n to recover the result like before.
        pval_beta_hat = pval_matrix[1,1]
        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)
        #length(beta_none0)>1
      } else{

        X_S <- x[, beta_none0]
        beta_hat <- solve(t(X_S) %*% X_S) %*% (t(X_S) %*% y)
        #residual
        U = y - X_S %*% beta_hat
        #new x
        X_SC <- x[, -beta_none0]
        x_new <- cbind(X_S, X_SC)
        x <- x_new
        beta_hat <- c(beta_hat, rep(0, (p - length(beta_none0))))
        #residual
        #U=y-x%*%beta_n
        #residual matrix ，ui*uj
        errormat <- U %*% t(U)
        #bandwidth
        len <- length(beta_none0)
        h <- b0 * (n ^ (-1 / (4 + len)))
        # 100 projections
        mu1 <- rep(0, p)
        sigma1 <- diag(1, p)
        beta_pro <- MASS::mvrnorm(n = np,
                            mu = mu1,
                            Sigma = sigma1)
        #beta1 <- beta_n
        beta_pro <- rbind(beta_hat, beta_pro)
        #p-value matrix
        pval_matrix <- matrix(nrow = np+1, ncol = 1) # 1 for beta_n

        for (k in 1:(np + 1)) {
          beta_pro[k, ] <- beta_pro[k, ] / sqrt(sum(beta_pro[k, ] ^ 2))
        }
        #p-value matrix
        pval_matrix <- matrix(nrow = np + 1, ncol = 1)
        for (q in 1:(np + 1)) {
          X_beta <- x %*% beta_pro[q, ]
          #kernel function matrix
          X_beta_mat <-
            ((X_beta) %*% matrix(1, 1, n) - matrix(1, n, 1) %*% (t(X_beta))) / h
          kermat <- (1 / sqrt(2 * pi)) * exp(-(X_beta_mat ^ 2) / 2)
          #test statistics
          Tn <-
            (sum(kermat * errormat) - psych::tr(kermat * errormat)) / sqrt(2 * (sum((
              kermat * errormat
            ) ^ 2) - psych::tr((
              kermat * errormat
            ) ^ 2)))
          pval <- 1 - pnorm(Tn)
          pval_matrix[q, ] <- pval
        }

        # beta_n
        pval_matrix = matrix(pval_matrix,ncol = 1) # remove the beta_n to recover the result like before.
        pval_beta_hat = pval_matrix[1,1]
        #cauchy combination
        pval_cauchy <- 1 - pcauchy(mean(tan((
          0.5 - pval_matrix
        ) * pi)))
        #single
        pval_single <- pval_matrix[2, 1]

        pval_matrix[pval_matrix==0,]=1e-15
        pval_hmp <- harmonicmeanp::pharmonicmeanp((np+1)/sum(1/pval_matrix),np+1)

      }
      #power
      result2 <- c( pval_cauchy, pval_single,pval_hmp,pval_beta_hat)
      LST_pval_cauchy[id_h] = result2[1]
      LST_pval_single[id_h] = result2[2]
      LST_pval_hmp[id_h] = result2[3]
      LST_pval_beta_hat[id_h] = result2[4]
      #error situation
    }, error = function(e) {
      print(e)
    })
  }

  #result
  return(
    data.frame(
      h = b0_tmp,
      T_cauchy = LST_pval_cauchy,
      T_alpha = LST_pval_single,
      T_hmp = LST_pval_hmp,
      T_beta = LST_pval_beta_hat
    )
  )
}


#' Model checking for high dimensional generalized linear models based on random projections
#'
#' The function can test goodness-of-fit of a low- or high-dimensional
#' generalized linear model (GLM) by detecting the presence of nonlinearity in
#' the conditional mean function of y given x using the statistics proposed by paper xx.
#' The outputs are p-value of  statisitics.
#'
#' @param y : y Input matrix with \code{n} rows, 1-dimensional response vector
#' @param x : x Input matrix with \code{n} rows, each a \code{p}-dimensional observation vector.
#' @param family : Must be "gaussian" or "binomial" for linear or logistic regression model.
#' @param b0 : a paramter to set bindwith, the default value may better for real data analysing.
#' @param np : the number of random projections.
#'
#' @return a list with five parameters returned. \code{h} stand for \code{b_0}.
#' T_alpha: the p value of our statistics by random projection. T_beta: the p value of our statistic by
#'  estimated projection. T_cauchy and T_hmp are p value of two combinational method proposed by
#'  Liu and Xie (2020) and Wilson (2019) respectively. each method combines p values of \code{np} random
#'  projections.
#' @export
#'
#' @examples
#'
#' set.seed(100)
#' data("sonar_mines")
#' x = sonar_mines[,-1]
#' y = sonar_mines$y
#'
#' ## make y as 0 or 1 for logistic regression
#' class1 = "R"
#' class2 ="M"
#' y = as.character(y)
#' y[y==class1]=1
#' y[y==class2]=0
#' y = as.numeric(y)
#' y = matrix(y,ncol = 1)
#'
#' ## scale x  and make data to be matrix
#' data_test_x = x
#' data_test_x = as.matrix(data_test_x)
#' data_test_y = as.matrix(y)
#' data_test_x = scale(data_test_x)
#' PLStests(data_test_y,data_test_x,family="binomial")
#'
#' ## add power2
#' x_x_2 = cbind(data_test_x,data_test_x^2)
#' x_x_2 = scale(x_x_2)
#' PLStests(data_test_y,x_x_2,family="binomial")
#'
PLStests <- function(y,x,family,b0=2,np=10){

  if (is.data.frame(x)) {
    message("Input is a data.frame. Converting to matrix...")
    x = as.matrix(x,ncol=1)
  } else if (is.matrix(x)) {
    message("Input x is already a matrix.")
  } else {
    stop("Input must be a data.frame or a matrix.")
  }

  if (is.data.frame(y)) {
    message("Input is a data.frame. Converting to matrix...")
    y = as.matrix(y,ncol=1)
  } else if (is.matrix(x)) {
    message("Input y is already a matrix.")
  } else {
    stop("Input must be a data.frame or a matrix.")
  }

  if (dim(x)[1]!=dim(y)[1]){
    stop("The rows of inputs x and y must be equal")
  }


  if(family=="gaussian"){
    result = PLStest_LM(y,x,b0,np)
  }

  if(family=="binomial"){
    result = PLStest_GLM(y,x,b0,np)
  }

  result <- as.list(result)


  return(result)
}
