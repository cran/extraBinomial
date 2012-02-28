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
A set of simulated pooled sequencing data can be generated using this package.
}
\details{
\tabular{ll}{
Package: \tab extraBinomial\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2012-02-28\cr
License: \tab GPL-3\cr
}
To use the function exbio, simply define two matrices R, R.alt with the
same dimensions (rows index SNPs and columns index pools), a vector cc
indicating the case and control status, number of chromosomes (n) and
then do: exbio(R, R.alt, cc, n) to yield the estimated allele frequencies and p-value based on extra-binomial model.
To use the function mysim, simply define the number of persons per pool (n), a vector cc
indicating the case and control status, number of SNPs needed (n.snp)
and their minor allele frequencies in controls (maf.ctl), assumed
sequencing error rate (e), relative risk ratio (rr), pool count scale parameter (z) and then do:
mysim(n,cc, n.snp,maf.ctl, e, rr,z) to yield two matrices R.err and
R.alt.err which give you the information about pool reads for both
alleles at each SNP position.
}
\author{
Xin Yang

Maintainer: Xin Yang <xin.yang@cimr.cam.ac.uk>
}
\references{
Yang et al. "Extra-binomial variation approach for reanalysis of pooled
DNA sequencing data", under review.
}
\keyword{ package }
