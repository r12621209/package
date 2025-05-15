#' Simulate Cross-Validation in Three Environments
#' 
#' Simulates 5 environments with 5-fold and 8-fold cross-validation (CV-O, CV-NO).
#'
#' @param omgaG_v0  A diagonal variance-covariance matrix across environments.
#' @param Y         A data.frame containing "Trait", "Name", and "ENV".
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
#' data(DST2_maize)
#' K.t2 <- kinship(snp_DST2_maize)
#' Yt2.s <- data.frame(
#'   Trait = as.vector(phe_DST2_maize),
#'   Name = factor(rep(rownames(K.t2), ncol(phe_DST2_maize))),
#'   Env  = factor(rep(seq_len(ncol(phe_DST2_maize)), each = nrow(phe_DST2_maize)))
#' )
#' mu <- c(10,15,20,25,30)
#' omgaG_v0=diag(c(1,1.5,2,2.5,3))
#' sim_cvO_8_5E <- sim_5E(omgaG_v0 = omgaG_v0, Y = Yt2.s, mu = mu, rho = 0.8, K = K.t2, cv = 2,fold = 8,random = 0)
#'


sim_5E=function(omgaG_v0, Y, mu, rho, K, cv, fold, random){
  
  library(doParallel)
  library(sommer)
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
  
  {count_res <- function(ENV1,ENV2,ENV3,ENV4,ENV5,l){
    
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
    
    
    ### AGBLUP modelA=NULL 
    A <- mmes(Trait ~ Env,
              random= ~vsm(ism(Env:Name), Gu=EA),              
              rcov=~ vsm(dsm(Env),ism(units)),
              data = Yna, verbose = FALSE)
    
    ### WGBLUP model
    W=NULL 
    W <- mmes(Trait ~ Env,
              random= ~ vsm(dsm(Env),ism(Name), Gu = K),
              rcov=~ vsm(dsm(Env),ism(units)),
              data = Yna, verbose = FALSE)
    
    MGE=NULL 
    MGE <- tryCatch({
      mmes(Trait ~ Env,
           random = ~ vsm(ism(Name), Gu = K) +  vsm(dsm(Env),ism(Name), Gu = K),
           rcov=~ vsm(dsm(Env),ism(units)),
           data = Yna, verbose = FALSE, tolParInv = 1e-06)
    }, error = function(e) {
      for (rep in 1:50) {
        tolParInv <- 1e-06 + mge_tolParInv * rep
        MGE <- tryCatch({
          mmes(Trait ~ Env,
               random = ~ vsm(ism(Name), Gu = K) +  vsm(dsm(Env),ism(Name), Gu = K),
               rcov=~ vsm(dsm(Env),ism(units)),
               data = Yna, verbose = FALSE, tolParInv = tolParInv)
        }, error = function(e) NULL)
        if (!is.null(MGE)) break
      }
      MGE
    })
    
    M=NULL 
    M <- tryCatch({
      mmes(Trait ~ Env,
           random = ~ vsm(usm(Env),ism(Name), Gu = K),
           rcov=~ vsm(dsm(Env),ism(units)),
           data = Yna, verbose = FALSE,
           tolParInv = 1e-06)
    }, error = function(e) {
      for (rep in 1:50) {
        tolParInv <- 1e-06 + m_tolParInv * rep
        M <- tryCatch({
          mmes(Trait ~ Env,
               random = ~ vsm(usm(Env),ism(Name), Gu = K),
               rcov=~ vsm(dsm(Env),ism(units)),
               data = Yna, verbose = FALSE,
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
    

    Z0 <- rep(c(A[["b"]][1], A[["b"]][1] + A[["b"]][-1]), each = k)
    Z1 <- A[["u"]]
    res[1]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[2]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[3]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[4]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[5]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    
    Z0 <- rep(c(W[["b"]][1], W[["b"]][1] + W[["b"]][-1]), each = k)
    Z1 <- W[["u"]]
    res[6]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[7]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[8]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[9]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[10]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    

    Z0 <- rep(c(MGE[["b"]][1], MGE[["b"]][1] + MGE[["b"]][-1]), each = k)
    Z1 <- unlist(MGE[["u"]][1:k]) + unlist(MGE[["u"]][-(1:k)])
    res[11]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[12]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[13]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[14]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[15]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    Z0 <- rep(c(M[["b"]][1], M[["b"]][1] + M[["b"]][-1]), each = k)
    Z1 <- M[["u"]]
    res[16]=cor((Z0+Z1)[ENV1],l$TBV[ENV1])
    res[17]=cor((Z0+Z1)[ENV2],l$TBV[ENV2])
    res[18]=cor((Z0+Z1)[ENV3],l$TBV[ENV3])
    res[19]=cor((Z0+Z1)[ENV4],l$TBV[ENV4])
    res[20]=cor((Z0+Z1)[ENV5],l$TBV[ENV5])
    
    ### Ensure M model variance is not zero
    var.diag <- diag(M[["theta"]][["vsm(usm(Env), ism(Name), Gu = K)"]])
    if(!is.null(M)){     
      if (var.diag[1]==0) {
        res[16]=NA
      }
      if (var.diag[2]==0) {
        res[17]=NA
      }
      if (var.diag[3]==0) {
        res[18]=NA
      }
      if (var.diag[4]==0) {
        res[19]=NA
      }
      if (var.diag[5]==0) {
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
      l <- NULL
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
      res=count_res(ENV1,ENV2,ENV3,ENV4,ENV5,l)
      res
    }
    if(i==1){result=answer}else{
      result=rbind(result,answer)
    }
    
    stopCluster(cl)
  }
  return(result)
}
