
#' Simulates 5 environments with 5-fold and 8-fold cross-validation (CV-O, CV-NO).
#'
#' @param omgaG_v0  A diagonal variance-covariance matrix across environments.
#' @param Y         A data.frame containing "trait", "environment", and "name".
#' @param mu        The mean value for each environment.
#' @param rho       The genetic correlation between environments.
#' @param K         An additive genetic relationship matrix.
#' @param cv        Different cross-validation problems, CV-O=1 and CV-NO=2.
#' @param fold      Different cross-validation methods, such as 5-fold and 8-fold.
#' @param random    A random seed number.
#' @export
#' @import sommer
#' @import doParallel
#' @import MASS
#' @import dplyr
#' @examples
#'
#' mu <- c(10,15,20,25,30)
#' omgaG_v0=diag(c(1,1.5,2,2.5,3))
#' sim_cvO_8_5E <- sim_5E(omgaG_v0 = omgaG_v0, Y = Yt2.s, mu = mu, rho = 0, K = K.t2, cv = 1,fold = 5,random = 0)
#'


sim_5E=function(omgaG_v0, Y, mu, rho, K, cv, fold, random){
  
  ### set prior parameters
  k=nrow(K)
  mu <- rep(mu,each=k)
  
  ### Compute the covariance matrix and the diagonal matrix
  omgaG <- omgaG_v0
  index_pairs <- combn(1:5, 2)
  
  for (i in 1:ncol(index_pairs)) {
    a <- index_pairs[1, i]
    b <- index_pairs[2, i]
    omgaG[a, b] <- omgaG[b, a] <- rho * sqrt(omgaG[a, a] * omgaG[b, b])
  }
  
  omgaE <- diag(c((1-0.2)/0.2*diag(omgaG)[1],(1-0.2)/0.2*diag(omgaG)[2],(1-0.5)/0.5*diag(omgaG)[3],
                  (1-0.5)/0.5*diag(omgaG)[4],(1-0.8)/0.8*diag(omgaG)[5]))
  
  SigmaG=kronecker(omgaG, K)
  SigmaE=kronecker(omgaE, diag(k))
  
  {count_res <- function(ENV1,ENV2,ENV3,ENV4,ENV5){
    
    if (k==301) {
      m_tolParInv <- 0.01 ; mge_tolParInv <- 1e-06
    }else{
      m_tolParInv <- 0.05 ; mge_tolParInv <- 5e-04
    }
    
    r1 <- c(ENV1,ENV2,ENV3,ENV4,ENV5)
    Yna <- l
    Yna[r1, 1] <- NA
    
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
    W=NULL 
    W <- mmer(Trait ~ Env,
              random = ~ vsr(dsr(Env), Name, Gu = K),
              rcov = ~ vsr(dsr(Env), units),
              method = 'AI',
              nIters = 10000,
              data = Yna, verbose = FALSE)
    
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
        tolParInv <- 1e-06 + mge_tolParInv * rep
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
        tolParInv <- 1e-06 + m_tolParInv * rep
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
    
    res <- matrix(NA, nrow = 1, ncol = 20)
    colnames(res) <- c("AGBLUP-1", "AGBLUP-2", "AGBLUP-3", "AGBLUP-4", "AGBLUP-5",
                       "WGBLUP-1", "WGBLUP-2", "WGBLUP-3", "WGBLUP-4", "WGBLUP-5",
                       "MGE-1", "MGE-2", "MGE-3", "MGE-4", "MGE-5",
                       "MGBLUP-1", "MGBLUP-2", "MGBLUP-3", "MGBLUP-4", "MGBLUP-5")
    
    #A
    Z0=rep(c(A[["Beta"]][["Estimate"]][1],A[["Beta"]][["Estimate"]][1]+A[["Beta"]][["Estimate"]][-1]),each=k)
    Z1=A[["U"]][["u:Env:Name"]][["Trait"]]
    
    res[1]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[2]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[3]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[4]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[5]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    
    #W
    Z0=rep(c(W[["Beta"]][["Estimate"]][1],W[["Beta"]][["Estimate"]][1]+W[["Beta"]][["Estimate"]][-1]),each=k)
    Z1=unlist(W[["U"]])
    res[6]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[7]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[8]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[9]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[10]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    #MGE
    Z0=rep(c(MGE[["Beta"]][["Estimate"]][1],MGE[["Beta"]][["Estimate"]][1]+MGE[["Beta"]][["Estimate"]][-1]),each=k)
    Z1=unlist(MGE[["U"]][-1])+unlist(MGE[["U"]][1])
    res[11]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[12]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[13]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[14]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[15]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    #M
    Z0=rep(c(M[["Beta"]][["Estimate"]][1],M[["Beta"]][["Estimate"]][1]+M[["Beta"]][["Estimate"]][-1]),each=k)
    a=M[["U"]][[1]][["Trait"]]+(M[["U"]][[2]][["Trait"]]+M[["U"]][[4]][["Trait"]]+M[["U"]][[7]][["Trait"]]+M[["U"]][[11]][["Trait"]])[-(1:k)]
    b=M[["U"]][[3]][["Trait"]]+M[["U"]][[2]][["Trait"]][1:k]+(M[["U"]][[5]][["Trait"]]+M[["U"]][[8]][["Trait"]]+M[["U"]][[12]][["Trait"]])[-(1:k)]
    c=M[["U"]][[6]][["Trait"]]+(M[["U"]][[4]][["Trait"]]+M[["U"]][[5]][["Trait"]])[1:k]+(M[["U"]][[9]][["Trait"]]+M[["U"]][[13]][["Trait"]])[-(1:k)]
    d=M[["U"]][[10]][["Trait"]]+(M[["U"]][[7]][["Trait"]]+M[["U"]][[8]][["Trait"]]+M[["U"]][[9]][["Trait"]])[1:k]+(M[["U"]][[14]][["Trait"]])[-(1:k)]
    e=M[["U"]][[15]][["Trait"]]+(M[["U"]][[11]][["Trait"]]+M[["U"]][[12]][["Trait"]]+M[["U"]][[13]][["Trait"]]+M[["U"]][[14]][["Trait"]])[1:k]
    
    Z1=c(a,b,c,d,e)
    
    res[16]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[17]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[18]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[19]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[20]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    if(!is.null(M)){  
      if (M[["sigma"]][["1:Name"]]==0) {
        res[16]=NA
      }
      if (M[["sigma"]][["2:Name"]]==0) {
        res[17]=NA
      }
      if (M[["sigma"]][["3:Name"]]==0) {
        res[18]=NA
      }
      if (M[["sigma"]][["4:Name"]]==0) {
        res[19]=NA
      }
      if (M[["sigma"]][["5:Name"]]==0) {
        res[20]=NA
      }
    }
    return(res)
  }}
  
  for (i in 1:30) {
    cat(i,"...")
    set.seed(i*random)
    ### Simulate genetic effect
    ghat <- mvrnorm(n=1, mu=rep(0, nrow(SigmaG)),Sigma = SigmaG)
    
    ### Start parallel computation
    cores <- detectCores()
    cl <- makeCluster(cores-1)
    registerDoParallel(cl)
    
    answer <- foreach(j = 1:30, .combine = rbind, .packages = c('MASS', 'sommer', 'dplyr')) %dopar% {
      
      ### Simulate error effect
      ehat <- mvrnorm(n = 1, mu = rep(0, nrow(SigmaE)), Sigma = SigmaE)
      ### Create data frame
      l <- data.frame(Trait=(mu + ghat + ehat),TBV=(mu + ghat),Name=rep(Y[1:k,'Name'],5),Env=rep(1:5,each=k))
      l$Name <- as.factor(l$Name)
      l$Env <- as.factor(l$Env)
      
      l <- l %>%
        arrange(Env, Name)
      ### Set different cross-validation combinations
      if (cv == 1 & fold == 5) {
        sets <- c(rep(1:5, k / 5), seq(k - floor(k / 5) * 5))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 1)+k
        ENV3 <- which(sets != 1)+2*k
        ENV4 <- which(sets != 1)+3*k
        ENV5 <- which(sets != 1)+4*k
      } else if (cv == 2 & fold == 5) {
        sets <- c(rep(1:5, k / 5), seq(k - floor(k / 5) * 5))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 2)+k
        ENV3 <- which(sets != 3)+2*k
        ENV4 <- which(sets != 4)+3*k
        ENV5 <- which(sets != 5)+4*k
      } else if (cv == 1 & fold == 8) {
        sets <- c(rep(1:8, k / 8), seq(k - floor(k / 8) * 8))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 1)+k
        ENV3 <- which(sets != 1)+2*k
        ENV4 <- which(sets != 1)+3*k
        ENV5 <- which(sets != 1)+4*k
      } else if (cv == 2 & fold == 8) {
        sets <- c(rep(1:8, k / 8), seq(k - floor(k / 8) * 8))
        sets <- sets[order(runif(k))]
        ENV1 <- which(sets != 1)
        ENV2 <- which(sets != 2)+k
        ENV3 <- which(sets != 3)+2*k
        ENV4 <- which(sets != 4)+3*k
        ENV5 <- which(sets != 5)+4*k
      } else {
        stop("Invalid combination of cv and fold. Allowed values: (cv=1, fold=5), (cv=2, fold=5), (cv=1, fold=8), (cv=2, fold=8).")
      }
      res=count_res(ENV1,ENV2,ENV3,ENV4,ENV5)
      res
    }
    if(i==1){result=answer}else{
      result=rbind(result,answer)
    }
    
    stopCluster(cl)
  }
  return(result)
}
