\name{exbio}
\alias{exbio}

\title{
Extra-binomial approach for pooled sequencing data
}
\description{
This funtion tests for differences in minor allele frequency between
groups and is based on extra-binomial variation model for pooled sequencing data. 
}
\usage{
exbio(R, R.alt, cc, n, tol = 0.001, a.start = 1, b.start = 1, max.it =
1000, digits = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
 A matrix with rows indexed by SNPs and columns by pools. The entries are counts of allele 1.
}
  \item{R.alt}{
A similarly formatted matrix containing the counts of allele 2.
}
  \item{cc}{
A case/control indicator vector with length = number of pools containing 0s (control pool) and 1s (case pool).
}
  \item{n}{
Number of chromosomes (twice the number of subjects) in each pooled sample.
}
  \item{tol}{
Maximum difference between coefficient values in successive glm before
we can stop, the default=0.001.
}
  \item{a.start}{
An intial value for the parameter a in linear regression, the default=1.
}
  \item{b.start}{
An intial value for the parameter b in linear regression, the default=1.
}
  \item{max.it}{
Maximum iterations, the default=1000.
}
  \item{digits}{
How many significant digits are to be used for allele frequency and p-value. The default, 'NULL', uses 'getOption(digits)'.
}
}
\details{
R and R.alt contain the read counts for the major allele
  and the alternative allele respectively and are required to have the
  same dimension. 

  The extra-binomial model defined:
  E(R/N)=p,
  Var(R/N)=p(1-p)(a/n+b/N) when N=R+R.alt

  We denote: W=1/(a/n+b/N), which may be interpreted as the adjusted
  depth of pool j for SNP i. Given the expected quantities:
  E(r2)=1/W=a/n+b/N,  the parameters a and b can be estimated by linear
  regression of r2 on 1/N, giving a/n as the intercept and b as the
  slope. This regression is carried out using generalized linear model
  (GLM) by first adopting Gaussian errors to estimate a relatively good
  start value of a and b, and then using these start values to do GLM
  with gamma errors and identity link because both a and b are
  positive.
  
  Since the estimated allele frequency p depends on a and b, the calculations are carried out iteratively.

  A chi-square test is performed on a 2*2 table using the weighted
  allele counts to calculate the p-value. 

}
\value{
A list containing the following components:
 \item{result}{a data.frame with three columns: the first shows the
   minor allele frequency of controls; the second shows the minor
   allele freqeuncy of cases; the third shows the p-value. Each row stands for a SNP.}  
\item{parameters}{a character vector indicating the values of the
  parameters a and b in the linear regression and and the times of iteration.}
}
\references{
Yang et al. "Extra-binomial variation approach for reanalysis of pooled
DNA sequencing data", under review.
}
\author{
Xin Yang
}

\examples{
R<-matrix(c(1409,1530,1490,1630,924,998,1000,1012),nrow=2,ncol=4,byrow=TRUE)
R.alt<-matrix(c(170,210,192,209,13,14,30,38),nrow=2,ncol=4,byrow=TRUE)
cc<-c(0,0,1,1)
n=96
exbio(R, R.alt, cc, n, max.it = 100, digits=3)
##=> p.value = 9.91e-01 for SNP1 and 4.01e-11 for SNP2,
##so association for SNP2 is established, but not for SNP1.
}
