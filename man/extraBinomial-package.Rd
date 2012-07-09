\name{extraBinomial-package}
\alias{extraBinomial-package}
\alias{extraBinomial}
\docType{package}
\title{
Extra-binomial approach for pooled sequencing data
}
\description{
This package tests for differences in minor allele frequency between
groups and is based on extra-binomial variation model for pooled sequencing data. 
}
\details{
\tabular{ll}{
Package: \tab extraBinomial\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2012-07-09\cr
License: \tab GPL-3\cr
}
To use the function exbio, simply define two matrices R, R.alt with the
same dimensions (rows index SNPs and columns index pools), a vector cc
indicating the case and control status, number of chromosomes (n) and
then do: exbio(R, R.alt, cc, n) to yield the estimated allele frequencies and p-value based on extra-binomial model.
}

\author{
Xin Yang, Chris Wallace

Maintainer: Xin Yang <xin.yang@cimr.cam.ac.uk>
}
\references{
Yang et al. "Extra-binomial variation approach for analysis of pooled DNA sequencing data", under review.
}
\keyword{ package }

