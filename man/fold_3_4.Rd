% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fold_3_4.R
\name{fold_3_4}
\alias{fold_3_4}
\title{Perform Cross-Validation on Real Data}
\usage{
fold_3_4(Y, K, cv, fold, random)
}
\arguments{
\item{Y}{A data.frame containing "Trait", "Name", and "ENV".}

\item{K}{An additive genetic relationship matrix.}

\item{cv}{Different cross-validation problems, CV-O=1 and CV-NO=2.}

\item{fold}{Different cross-validation methods, such as 5-fold and 8-fold.}

\item{random}{A random seed number.}
}
\description{
Applies 3-fold and 4-fold cross-validation on real phenotypic and genotypic data.
Supports both CV-O and CV-NO.
}
\examples{

data(DST2_maize)
K.t2 <- kinship(snp_DST2_maize)
Yt2.s <- data.frame(
 Trait = as.vector(phe_DST2_maize),
 Name = factor(rep(rownames(K.t2), ncol(phe_DST2_maize))),
 Env = factor(rep(seq_len(ncol(phe_DST2_maize)), each = nrow(phe_DST2_maize)))
)
fold_3_cvO <- fold_3_4(Y=Yt2.s,K=K.t2,cv=1,fold=3,random=0)
}
