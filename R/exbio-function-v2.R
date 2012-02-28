##This funtion yields minor allele frequencies and p-values based on extra-binomial variation model for pooled sequencing data. 

##n:           number of chromosomes (twice the number of subjects) in each pooled sample
##tol:         maximum difference between coefficient values in successive glm before we can stop, the default = 1e-3.
##max.it:    maximum iterations.
##digits:     how many significant digits are to be used for allele frequency and p-value.  The default, ‘NULL’, uses ‘getOption(digits)’.
## R:         a matrix with rows indexed by SNPs and columns by pools.  the entries are counts of allele 1
## R.alt:     a similarly formatted matrix containing the counts of allele 2
## cc:         a  case/control indicator vector with length = number of pools containing 0s (control pool) and 1s (case pool).
##a.start:     an intial value for the parameter a in linear regression, the default=1.
##b.start:     an intial value for the parameter a in linear regression, the default=1.
## helper functions

calcp <- function(W,R,N,wh.cse,wh.ctl) {
  WRN <- W*R/N
  p<-matrix(nrow=nrow(WRN), ncol=2)
  p[,1] <- apply(WRN[,wh.ctl],1,sum,na.rm=TRUE) / apply(W[,wh.ctl],1,sum,na.rm=TRUE)    
  p[,2] <- apply(WRN[,wh.cse],1,sum,na.rm=TRUE) / apply(W[,wh.cse],1,sum,na.rm=TRUE)
  p[p>1] <- 1
  return(p)
}


exbio<-function(R,R.alt,cc,n,tol=1e-3,a.start=1,b.start=1,max.it=1000,digits=NULL){
  
  ## some checks
  n.pools <- length(cc)
  if(ncol(R)!=n.pools | ncol(R.alt)!=n.pools)
    stop("R, R.alt should have ncol = length(cc)")
  n.snps <- nrow(R)
  if(nrow(R.alt)!=n.snps)
    stop("R.alt should have nrow = nrow(R)")
  cc <- as.integer(cc)
  if(!identical(sort(unique(cc)), as.integer(c(0,1))))
    stop("cc should contain only 0 (control pools) and 1 (case pools)")
  
  ## Identify cases and controls
  wh.cse <- which(cc==1)
  wh.ctl <- which(cc==0)
  n.cse <- length(wh.cse)
  n.ctl <- length(wh.ctl)
  
  ##Coverage(total number of reads of both alleles of a SNP)
  N<-R+R.alt
  ##weight
  W<-matrix(nrow=n.snps, ncol=(n.ctl+n.cse))
  ##p[,1]--the allele freqency in control;p[,2]--the allele frequency in case
  p<-matrix(nrow=n.snps, ncol=2)
  ##r2--the quantities
  r2<-matrix(nrow=n.snps, ncol=n.pools)
  ab<-matrix(nrow=max.it,ncol=2)
  para <- NULL

###initial estimates
  k=0   # the times of errors
  j=1
  a <- b <- 1
  ab[j,1]<- ab[j,2]<-1
###Adopt gaussian errors to generate an initial values of a and b
  repeat{
    
  W<-1/(a/n+b/N)
  ##p[,1]stores the efficient estimate of the allele frequency in control         
  ##p[,2]stores the efficient estimate of the allele frequency in case        
  p <- calcp(W,R,N,wh.cse,wh.ctl)

  ##estimate the quantities   
  r2[,wh.ctl]<-n.ctl/(n.ctl-1)/(p[,1]*(1-p[,1]))*(R[,wh.ctl]/N[,wh.ctl]-p[,1])^2
  r2[,wh.cse]<-n.cse/(n.cse-1)/(p[,2]*(1-p[,2]))*(R[,wh.cse]/N[,wh.cse]-p[,2])^2

  ## generate variables for regression
  y <- as.numeric(r2)
  x1 <- as.numeric(1/N) # 1/N - to estimate b

  ##Adopt gaussian errors to generate an initial values of a and b
  
  para<-glm(y~x1,family=gaussian, subset=y<40 )
  a<-para$coefficients[1]*n
  b<-para$coefficients[2]
  #cat(j,"\n")
  #print(para$coefficients)
   j<-j+1
   ab[j,1]<-a
   ab[j,2]<-b
  if ((a>0 &b>0) | max(abs(ab[j,]-ab[(j-1),]))<=1e-5 |  j==max.it ) {
      break
    }  
  }
  if(a<0) { a<-1}
  if(b<0) {b<-1}
  j=1
### iterate
  repeat {
    
    W<-1/(a/n+b/N)
    p <- calcp(W,R,N,wh.cse,wh.ctl)
    r2[,wh.ctl]<-n.ctl/(n.ctl-1)/(p[,1]*(1-p[,1]))*(R[,wh.ctl]/N[,wh.ctl]-p[,1])^2
    r2[,wh.cse]<-n.cse/(n.cse-1)/(p[,2]*(1-p[,2]))*(R[,wh.cse]/N[,wh.cse]-p[,2])^2
    
    y <- as.numeric(r2)
    x1 <- as.numeric(1/N)
    
    para.old <- para
    
    ##GLM with gamma errors given a and b are positive.   
    para<-try(glm(y~x1,family=Gamma(link="identity"),start=c(a/n,b), subset=y<40,control = list(maxit = 70) ),silent=TRUE)
    if (inherits(para,"try-error")) {
      if (k<=10){
      #reset the start values
       a<-ab[1,1]
       ab[1,2]<-b<-ab[1,2]+1
       j=1
       k=k+1
       para<-para.old}
      else{
        break
        result<-data.frame(control.maf=NA,
                     case.maf=NA,
                     p.value=NA)
        return(list(result=result,
              parameters=warning("inner loop 1; cannot correct step size")))}}
        
    else{
      if(sum(is.na(para$coefficients)) ){
          #reset the start values
          a<-ab[1,1]
          ab[1,2]<-b<-ab[1,2]+1
          j=1
          para<-para.old}
      else {
    #      cat(j,"\n")
    #      print(para$coefficients)
          a<-para$coefficients[1]*n
          b<-para$coefficients[2]
            ## converged?
              if((j>2& max(abs(para.old$coefficients - para$coefficients))<=tol) ) {
                 break
               }
              else{
                if (j==max.it) {
                  break
                  result<-data.frame(control.maf=NA,
                     case.maf=NA,
                     p.value=NA)
                  return(list(result=result,
                     parameters=warning("glm is not converged")))}
               }
          j<-j+1
          ab[j,1]<-a
          ab[j,2]<-b
          }
    }
  }

      

   
  
  W<-1/(a/n+b/N)
  p <- calcp(W,R,N,wh.cse,wh.ctl)
  
  ##weighted allele counts in each pool
  count1<-matrix(nrow=n.snps, ncol=2)
  count1[,1]<-p[,1]*apply(W[,wh.ctl],1,sum,na.rm=TRUE)
  count1[,2]<-p[,2]*apply(W[,wh.cse],1,sum,na.rm=TRUE)
  count2<-matrix(nrow=n.snps, ncol=2)
  count2[,1]<-(1-p[,1])*apply(W[,wh.ctl],1,sum,na.rm=TRUE)
  count2[,2]<-(1-p[,2])*apply(W[,wh.cse],1,sum,na.rm=TRUE)
  
  chi.sq<-numeric(n.snps)
  chisq.pvalue<-numeric(n.snps)
  
#chi-square test with Yales' s correction
chi.sq<-(abs(count1[,2]*count2[,1]-count2[,2]*count1[,1])-(count1[,2]+count2[,2]+count1[,1]+count2[,1])/2)^2*(count1[,2]+count2[,2]+count1[,1]+count2[,1])/((count1[,2]+count2[,2])*(count1[,1]+count2[,1])*(count2[,2]+count2[,1])*(count1[,2]+count1[,1]))
chisq.pvalue<-pchisq(chi.sq, 1, lower.tail=FALSE)
    
                
  result<-data.frame(control.maf=format(pmin(p[,1],1-p[,1]),digits=digits),
                     case.maf=format(pmin(p[,2],1-p[,2]),digits=digits),
                     p.value=format(as.numeric(chisq.pvalue),digits=digits))
  return(list(result=result,
              parameters=c(a=a,b=b,n.iterations=format(j,nsmall=0))))
}

