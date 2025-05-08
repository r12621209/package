#' Perform Cross-Validation on Real Data
#'
#' Applies 3-fold and 4-fold cross-validation on real phenotypic and genotypic data.
#' Supports both CV-O and CV-NO.
#'

#' @param Y         A data.frame containing "Trait", "Name", and "ENV".
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
#'data(DST2_maize)
#'K.t2 <- kinship(snp_DST2_maize)
#'Yt2.s <- data.frame(
#'  Trait = as.vector(phe_DST2_maize),
#'  Name = factor(rep(rownames(K.t2), ncol(phe_DST2_maize))),
#'  Env = factor(rep(seq_len(ncol(phe_DST2_maize)), each = nrow(phe_DST2_maize)))
#')
#' fold_3_cvO <- fold_3_4(Y=Yt2.s,K=K.t2,cv=1,fold=3,random=0)


fold_3_4 <- function(Y,K,cv,fold,random){

  library(doParallel)
  ### set prior parameters
  k=nrow(K)
  
  if (nlevels(Y$Env)==2) {
    count_res=function(ENV1,ENV2){
      
      ### Evaluate the model based on given ENV1, ENV2
      r1=c(ENV1,ENV2)
      ### Add missing values
      Yna <- l
      Yna[r1, 1] <- NA
      
      E <- diag(length(unique(Yna$Env)))
      rownames(E) <- colnames(E) <- unique(Yna$Env)
      EA <- kronecker(E, K, make.dimnames = TRUE)
      
      ### AGBLUP model
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
        for (rep in 1:200) {
          tolParInv <- 1e-06 + 0.2 * rep
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
      
      
      res <- matrix(NA, nrow = 1, ncol = 8)
      colnames(res) <- c("AGBLUP-1", "AGBLUP-2", "WGBLUP-1", "WGBLUP-2",
                         "MGE-1", "MGE-2", "MGBLUP-1", "MGBLUP-2")
      
      ### Compute correlation coefficients
      Z0 <- rep(c(A[["Beta"]][["Estimate"]][1], sum(A[["Beta"]][["Estimate"]][1:2])), each = k)
      Z1 <- A[["U"]][["u:Env:Name"]][["Trait"]]
      res[1]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[2]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])    
      
      Z0 <- rep(c(W[["Beta"]][["Estimate"]][1], W[["Beta"]][["Estimate"]][1] + W[["Beta"]][["Estimate"]][-1]), each = k)
      Z1 <- unlist(W[["U"]])
      res[3]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[4]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])    
      
      Z0 <- rep(c(MGE[["Beta"]][["Estimate"]][1], MGE[["Beta"]][["Estimate"]][1] + MGE[["Beta"]][["Estimate"]][-1]), each = k)
      Z1 <- unlist(MGE[["U"]][-1]) + unlist(MGE[["U"]][1])
      res[5]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[6]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])    
      
      Z0 <- rep(c(M[["Beta"]][["Estimate"]][1], sum(M[["Beta"]][["Estimate"]][1:2])), each = k)
      Z1 <- c(M[["U"]][["1:Name"]][["Trait"]] + M[["U"]][["2:1:Name"]][["Trait"]][-(1:k)],
              M[["U"]][["2:Name"]][["Trait"]] + M[["U"]][["2:1:Name"]][["Trait"]][(1:k)])
      res[7]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[8]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])    
      
      if(!is.null(M)){     
        if (M[["sigma"]][["1:Name"]]<=0) {
          res[7]=NA
        }
        if (M[["sigma"]][["2:Name"]]<=0) {
          res[8]=NA
        }
        
      }
      return(res)
    }
  }else if (nlevels(Y$Env)==3) {
    count_res=function(ENV1,ENV2,ENV3){
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
      res[1]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[2]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])
      res[3]=cor((Z0+Z1)[ENV3],l$Trait[ENV3])
      
      Z0 <- rep(c(W[["Beta"]][["Estimate"]][1], W[["Beta"]][["Estimate"]][1] + W[["Beta"]][["Estimate"]][-1]), each = k)
      Z1 <- unlist(W[["U"]])
      res[4]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[5]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])
      res[6]=cor((Z0+Z1)[ENV3],l$Trait[ENV3])
      
      Z0 <- rep(c(MGE[["Beta"]][["Estimate"]][1], MGE[["Beta"]][["Estimate"]][1] + MGE[["Beta"]][["Estimate"]][-1]), each = k)
      Z1 <- unlist(MGE[["U"]][-1]) + unlist(MGE[["U"]][1])
      res[7]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[8]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])
      res[9]=cor((Z0+Z1)[ENV3],l$Trait[ENV3])
      
      Z0=rep(c(M[["Beta"]][["Estimate"]][1],sum(M[["Beta"]][["Estimate"]][1:2]),sum(M[["Beta"]][["Estimate"]][c(1,3)])),each=k)
      Z1=c(M[["U"]][["1:Name"]][["Trait"]]+(M[["U"]][["2:1:Name"]][["Trait"]][-(1:k)])+(M[["U"]][["3:1:Name"]][["Trait"]][-(1:k)]),
           M[["U"]][["2:Name"]][["Trait"]]+(M[["U"]][["2:1:Name"]][["Trait"]][1:k])+(M[["U"]][["3:2:Name"]][["Trait"]][-(1:k)]),
           M[["U"]][["3:Name"]][["Trait"]]+(M[["U"]][["3:1:Name"]][["Trait"]][1:k])+(M[["U"]][["3:2:Name"]][["Trait"]][1:k]))
      res[10]=cor((Z0+Z1)[ENV1],l$Trait[ENV1])
      res[11]=cor((Z0+Z1)[ENV2],l$Trait[ENV2])
      res[12]=cor((Z0+Z1)[ENV3],l$Trait[ENV3])
      
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
    }
  }else {
    stop("Out of selectable environments, exiting the function.")
  }
  
  if (!(fold %in% c(3, 4))) {
    stop("fold must be 3 or 4.")
  }
  
  if (!(cv %in% c(1, 2))) {
    stop("cv must be 1 or 2.")
  }
  
  ### Start parallel computation
  cores <- detectCores()
  cl <- makeCluster(cores-1)
  registerDoParallel(cl)
  
  
  answer <- foreach(j = 1:500, .combine = rbind, .packages = c('MASS', 'sommer', 'dplyr')) %dopar% {
    cat(j,"...")
    
    ### Create data frame
    l <- Y
    l <- l %>%
      arrange(Env, Name)
    set.seed(j*random)
    
    ### For two environments
    if (nlevels(Y$Env) == 2) {
      
      ### Configure CV and K-fold combinations
      if (fold == 3) {
        sets <- c(rep(1:3, k %/% 3), seq(k - floor(k / 3) * 3))
        sets <- sets[order(runif(k))]
        
        ### Set different CV
        if (cv == 1) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 1) + k
          f_1 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 2) + k
          f_2 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 3) + k
          f_3 <- count_res(ENV1, ENV2)
          
          res <- rbind(f_1, f_2, f_3)
          
        } else if (cv == 2) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 2) + k
          f_1 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 3) + k
          f_2 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 1) + k
          f_3 <- count_res(ENV1, ENV2)
          
          res <- rbind(f_1, f_2, f_3)
          
        } else {
          stop("Out of cv.")
        }
        
      } else if (fold == 4) {
        sets <- c(rep(1:4, k %/% 4), seq(k - floor(k / 4) * 4))
        sets <- sets[order(runif(k))]
        
        if (cv == 1) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 1) + k
          f_1 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 2) + k
          f_2 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 3) + k
          f_3 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 4)
          ENV2 <- which(sets != 4) + k
          f_4 <- count_res(ENV1, ENV2)
          
          res <- rbind(f_1, f_2, f_3, f_4)
          
        } else if (cv == 2) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 2) + k
          f_1 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 3) + k
          f_2 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 4) + k
          f_3 <- count_res(ENV1, ENV2)
          
          ENV1 <- which(sets != 4)
          ENV2 <- which(sets != 1) + k
          f_4 <- count_res(ENV1, ENV2)
          
          res <- rbind(f_1, f_2, f_3, f_4)
          
        } 
      }
      
    }
    
    ### For three environments
    if (nlevels(Y$Env) == 3) {
      ### Configure CV and K-fold combinations
      if (fold == 3) {
        sets <- c(rep(1:3, k %/% 3), seq(k - floor(k / 3) * 3))
        sets <- sets[order(runif(k))]
        
        ### Set different CV
        if (cv == 1) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 1) + k
          ENV3 <- which(sets != 1) + 2 * k
          f_1 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 2) + k
          ENV3 <- which(sets != 2) + 2 * k
          f_2 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 3) + k
          ENV3 <- which(sets != 3) + 2 * k
          f_3 <- count_res(ENV1, ENV2, ENV3)
          
          res <- rbind(f_1, f_2, f_3)
          
        } else if (cv == 2) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 2) + k
          ENV3 <- which(sets != 3) + 2 * k
          f_1 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 3) + k
          ENV3 <- which(sets != 1) + 2 * k
          f_2 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 1) + k
          ENV3 <- which(sets != 2) + 2 * k
          f_3 <- count_res(ENV1, ENV2, ENV3)
          
          res <- rbind(f_1, f_2, f_3)
          
        } 
        
      } else if (fold == 4) {
        sets <- c(rep(1:4, k %/% 4), seq(k - floor(k / 4) * 4))
        sets <- sets[order(runif(k))]
        
        if (cv == 1) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 1) + k
          ENV3 <- which(sets != 1) + 2 * k
          f_1 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 2) + k
          ENV3 <- which(sets != 2) + 2 * k
          f_2 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 3) + k
          ENV3 <- which(sets != 3) + 2 * k
          f_3 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 4)
          ENV2 <- which(sets != 4) + k
          ENV3 <- which(sets != 4) + 2 * k
          f_4 <- count_res(ENV1, ENV2, ENV3)
          
          res <- rbind(f_1, f_2, f_3, f_4)
          
        } else if (cv == 2) {
          ENV1 <- which(sets != 1)
          ENV2 <- which(sets != 2) + k
          ENV3 <- which(sets != 3) + 2 * k
          f_1 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 2)
          ENV2 <- which(sets != 3) + k
          ENV3 <- which(sets != 4) + 2 * k
          f_2 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 3)
          ENV2 <- which(sets != 4) + k
          ENV3 <- which(sets != 1) + 2 * k
          f_3 <- count_res(ENV1, ENV2, ENV3)
          
          ENV1 <- which(sets != 4)
          ENV2 <- which(sets != 1) + k
          ENV3 <- which(sets != 2) + 2 * k
          f_4 <- count_res(ENV1, ENV2, ENV3)
          
          res <- rbind(f_1, f_2, f_3, f_4)
          
        } 
        
      }
      
    }
    res
  }
  
  return(answer)
  stopCluster(cl)
}
