library(doParallel)
#' Simulates 3 environments with 3-fold and 4-fold cross-validation (CV-O, CV-NO).
#'
#' @param omgaG_v0  A diagonal variance-covariance matrix across environments.
#' @param Y         A data.frame containing "trait", "environment", and "name".
#' @param mu        The mean value for each environment.
#' @param rho       The genetic correlation between environments.
#' @param K         An additive genetic relationship matrix.
#' @param cv        Different cross-validation problems, CV-O=1 and CV-NO=2.
#' @param fold      Different cross-validation methods, such as 3-fold and 4-fold.
#' @param random    A random seed number.
#' @export
#' @import sommer
#' @import doParallel
#' @import MASS
#' @import dplyr
#' @examples
#'
#' mu <- c(10,15,20)
#' omgaG_v0 <- diag(c(1,1.5,2))
#' sim_cvO_3_3E <- sim_3E(omgaG_v0 = omgaG_v0, Y = Yt2.s, mu = mu, rho = 0, K = K.t2, cv = 1,fold = 3,random = 0)





sim_3E = function(omgaG_v0, Y, mu, rho, K, cv, fold, random){
  
  ### set prior parameters
  k=nrow(K)
  mu <- rep(mu,each=k)
  
  ### Compute the covariance matrix and the diagonal matrix
  omgaG <- omgaG_v0
  omgaG[1,2] <- omgaG[2,1] <- rho * sqrt(omgaG[1,1] * omgaG[2,2])
  omgaG[1,3] <- omgaG[3,1] <- rho * sqrt(omgaG[1,1] * omgaG[3,3])
  omgaG[2,3] <- omgaG[3,2] <- rho * sqrt(omgaG[2,2] * omgaG[3,3])
  
  omgaE <- diag((1 - c(0.2, 0.5, 0.8)) / c(0.2, 0.5, 0.8) * diag(omgaG))
  
  SigmaG=kronecker(omgaG, K)
  SigmaE=kronecker(omgaE, diag(k))
  
  {count_res=function(ENV1,ENV2,ENV3){
    
    ### Evaluate the model based on given ENV1, ENV2, ENV3
    r1=c(ENV1,ENV2,ENV3)
    ### Add missing values
    Yna <- l
    Yna[r1, 1] <- NA
    
    ### AGBLUP model
    E <- diag(length(unique(Yna$Env)))
    rownames(E) <- colnames(E) <- unique(Yna$Env)
    EA <- kronecker(E, K, make.dimnames = TRUE)
    
    A=NULL 
    A <- mmer(Trait ~ Env,
              random = ~ vsr(Env:Name, Gu = EA),
              rcov = ~ vsr(dsr(Env), units),
              data = Yna, verbose = FALSE,
              method = 'AI',
              nIters = 10000)
    
    ### WGBLUP model
    W=NULL 
    W <- mmer(Trait ~ Env,
              random = ~ vsr(dsr(Env), Name, Gu = K),
              rcov = ~ vsr(dsr(Env), units),
              method = 'AI',
              nIters = 10000,
              data = Yna, verbose = FALSE)
    
    ### MGE model
    MGE=NULL 
    MGE <- tryCatch({
      mmer(Trait ~ Env,
           random = ~ vsr(Name, Gu = K) + vsr(dsr(Env), Name, Gu = K),
           rcov = ~ vsr(dsr(Env), units),
           data = Yna, verbose = FALSE,
           method = 'AI',
           nIters = 10000, tolParInv = 1e-06)
    }, error = function(e) {
      for (rep in 1:100) {
        tolParInv <- 1e-06 + 0.0005 * rep
        MGE <- tryCatch({
          mmer(Trait ~ Env,
               random = ~ vsr(Name, Gu = K) + vsr(dsr(Env), Name, Gu = K),
               rcov = ~ vsr(dsr(Env), units),
               data = Yna, verbose = FALSE,
               method = 'AI',
               nIters = 10000, tolParInv = tolParInv)
        }, error = function(e) NULL)
        if (!is.null(MGE)) break
      }
      MGE
    })
    
    ### MGBLUP model
    M=NULL 
    M <- tryCatch({
      mmer(Trait ~ Env,
           random = ~ vsr(usr(Env), Name, Gu = K),
           rcov = ~ vsr(dsr(Env), units),
           data = Yna,
           verbose = FALSE,
           nIters = 10000,
           method = 'AI',
           tolParInv = 1e-06)
    }, error = function(e) {
      for (rep in 1:100) {
        tolParInv <- 1e-06 + 0.05 * rep
        M <- tryCatch({
          mmer(Trait ~ Env,
               random = ~ vsr(usr(Env), Name, Gu = K),
               rcov = ~ vsr(dsr(Env), units),
               data = Yna,
               verbose = FALSE,
               nIters = 10000,
               method = 'AI',
               tolParInv = tolParInv)
        }, error = function(e) NULL)
        if (!is.null(M)) break
      }
      M
    })
    
    res <- matrix(NA, nrow = 1, ncol = 12)
    colnames(res) <- c("AGBLUP-1", "AGBLUP-2", "AGBLUP-3",
                       "WGBLUP-1", "WGBLUP-2", "WGBLUP-3",
                       "MGE-1", "MGE-2", "MGE-3",
                       "MGBLUP-1", "MGBLUP-2", "MGBLUP-3")
    
    ### Compute correlation coefficients
    Z0 <- rep(c(A[["Beta"]][["Estimate"]][1], A[["Beta"]][["Estimate"]][1] + A[["Beta"]][["Estimate"]][-1]), each = k)
    Z1 <- A[["U"]][["u:Env:Name"]][["Trait"]]
    res[1]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[2]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[3]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    
    Z0 <- rep(c(W[["Beta"]][["Estimate"]][1], W[["Beta"]][["Estimate"]][1] + W[["Beta"]][["Estimate"]][-1]), each = k)
    Z1 <- unlist(W[["U"]])
    res[4]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[5]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[6]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    
    Z0 <- rep(c(MGE[["Beta"]][["Estimate"]][1], MGE[["Beta"]][["Estimate"]][1] + MGE[["Beta"]][["Estimate"]][-1]), each = k)
    Z1 <- unlist(MGE[["U"]][-1]) + unlist(MGE[["U"]][1])
    res[7]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[8]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[9]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    
    Z0=rep(c(M[["Beta"]][["Estimate"]][1],sum(M[["Beta"]][["Estimate"]][1:2]),sum(M[["Beta"]][["Estimate"]][c(1,3)])),each=k)
    Z1=c(M[["U"]][["1:Name"]][["Trait"]]+(M[["U"]][["2:1:Name"]][["Trait"]][-(1:k)])+(M[["U"]][["3:1:Name"]][["Trait"]][-(1:k)]),
         M[["U"]][["2:Name"]][["Trait"]]+(M[["U"]][["2:1:Name"]][["Trait"]][1:k])+(M[["U"]][["3:2:Name"]][["Trait"]][-(1:k)]),
         M[["U"]][["3:Name"]][["Trait"]]+(M[["U"]][["3:1:Name"]][["Trait"]][1:k])+(M[["U"]][["3:2:Name"]][["Trait"]][1:k]))
    res[10]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[11]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[12]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    
    ### Ensure M model variance is not zero
    if(!is.null(M)){     
      if (M[["sigma"]][["1:Name"]]==0) {
        res[10]=NA
      }
      if (M[["sigma"]][["2:Name"]]==0) {
        res[11]=NA
      }
      if (M[["sigma"]][["3:Name"]]==0) {
        res[12]=NA
      }
    }
    return(res)
  }}
  
  for (i in 1:30) {
    cat(i,"...")
    set.seed(i*random)
    ### Simulate genetic effect
    ghat=mvrnorm(n=1,mu=rep(0,nrow(SigmaG)),Sigma = SigmaG)
    
    ### Start parallel computation
    cores <- detectCores()
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)
    
    answer <- foreach(j = 1:30, .combine = rbind, .packages = c('MASS', 'sommer', 'dplyr')) %dopar% {
      ### Simulate error effect
      ehat <- mvrnorm(n = 1, mu = rep(0, nrow(SigmaE)), Sigma = SigmaE)
      ### Create data frame
      l <- data.frame(Trait=(mu + ghat + ehat),TBV=(mu + ghat),Name=rep(Y[1:k,'Name'],3),Env=rep(1:3,each=k))
      l$Name <- as.factor(l$Name)
      l$Env <- as.factor(l$Env)
      
      l <- l %>%
        arrange(Env, Name)
      ### Set different cross-validation combinations
      if (cv == 1 & fold == 3) {
        sets <- c(rep(1:3, k / 3), seq(k - floor(k / 3) * 3))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 1) + k
        ENV3 <- which(sets != 1) + 2 * k
      } else if (cv == 2 & fold == 3) {
        sets <- c(rep(1:3, k / 3), seq(k - floor(k / 3) * 3))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 2) + k
        ENV3 <- which(sets != 3) + 2 * k
      } else if (cv == 1 & fold == 4) {
        sets <- c(rep(1:4, k / 4), seq(k - floor(k / 4) * 4))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 1) + k
        ENV3 <- which(sets != 1) + 2 * k
      } else if (cv == 2 & fold == 4) {
        sets <- c(rep(1:4, k / 4), seq(k - floor(k / 4) * 4))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 2) + k
        ENV3 <- which(sets != 3) + 2 * k
      } else {
        stop("Invalid combination of cv and fold. Allowed values: (cv=1, fold=3), (cv=2, fold=3), (cv=1, fold=4), (cv=2, fold=4).")
      }
      
      res=count_res(ENV1,ENV2,ENV3)
      res
    }
    
    if(i==1){result=answer}else{
      result=rbind(result,answer)
    }
    
    stopCluster(cl)
  }
  return(result)
}
