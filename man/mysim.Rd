\name{mysim}
\alias{mysim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating pooled sequencing data
}
\description{
This function simulates a set of pooled sequencing data and yields the
pool reads for both major and minor allele.
}
\usage{
mysim(n, cc, n.snp, maf.ctl, e, rr, z)
}
\arguments{
  \item{n}{
Number of persons per pool.
}
  \item{cc}{
A  case/control indicator vector with length = number of pools containing 0s (control pool) and 1s (case pool).
}
  \item{n.snp}{
Number of SNPs.
}
  \item{maf.ctl}{
Minor allele frequency in controls.
}
  \item{e}{
Sequencing error rate.
}
  \item{rr}{
Relative risk ratio.
}
\item{z}{
A scale parameter to keep the expected total reads per pool stable. z=1
when there are 20 pools with 50 persons per pool. z=2 if there are 40
pools with 25 persons per pool and so on.
  }
}
\details{
This function simulates pooled sequencing data based on Gamma
distribution. To further increase between-pool variation and make our
simulated data closer to biological data, we assume the shape parameters
followed a Gamma distribution and the rate parameters followed a Beta distribution. 
}
\value{
A list containing the following components:
 \item{R.err}{A matrix with rows indexed by SNPs and columns by
   pools. The entries are read counts of the major allele.}  
\item{R.alt.err}{A similarly formatted matrix containing the read counts
  of the minor allele.}
}
\references{
Yang et al. "Extra-binomial variation approach for reanalysis of pooled
DNA sequencing data", under review.
}
\author{
Xin Yang
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~simulation }

