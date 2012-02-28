# #This function simulates pooled sequencing data based on Gamma distribution.

##n: number of persons per pool
##m: number of pools
##maf.ctl: minor allele frequency in controls
##maf.cas: minor allele frequency in cases
##n.snp: number of SNPs
##cc: a  case/control indicator vector with length = number of pools containing 0s (control pool) and 1s (case pool).
##e: sequencing error rate
##rr: relative risk ratio
##z: a scale parameter to keep the expected total reads per pool stable. z=1 when there are 20 pools with 50 persons per pool.

mysim <- function(n,cc, n.snp,maf.ctl, e, rr,z) {
                                        
#number of pools
  m<-length(cc)
                    
#calculate minor allele frequencies in cases
  a<-rep(1,n.snp); b<-rep(-1,n.snp); c<-rr*maf.ctl*(1-maf.ctl)  
if(rr==1) {
  maf.cse <- maf.ctl
} else {
  maf.cse<-apply(rbind((-b + sqrt(b^2 - 4*a*c))/(2*a),(-b - sqrt(b^2 - 4*a*c))/(2*a)),2,min)
}

R.err<-matrix(0, nrow=n.snp, ncol=m) #major reads per pool considering error rates
R.alt.err<-matrix(0, nrow=n.snp, ncol=m)   #minor reads per pool considering error rates

m.case<-sum(cc==1)
m.control<-sum(cc==0)
  
  #we assume shape ~ Gamma(shape=1.5773,rate=0.1846) given 20 pools with 50 persons per pool.
  shape<-matrix(rgamma(n.snp*m,shape=1.5773*z,rate=0.1846),nrow=n.snp, ncol=m)
  #we assume rate ~ Beta(2.094,397.611) given our real data (20 pools with 50 persons per pool).
  rate<-matrix(rbeta(n.snp*m,2.094,397.611),nrow=n.snp, ncol=m)
 
  for (i in 1:n.snp) {
  #control pools
    for (j in 1:m.control){
      ##parameters keep the same within each pool

    ##controls
    #major allele reads 
    major.reads.person.ctl<-0.01* rgamma(n,shape=(1-maf.ctl[i])*shape[i,j],rate=rate[i,j])
    ## R:         counts of allele 1 in this pool
    R<-major.reads.pool.ctl<-sum(major.reads.person.ctl)
                            
    #minor allele reads:
    minor.reads.person.ctl<-0.01* rgamma(n,shape=maf.ctl[i]*shape[i,j],rate=rate[i,j])
    ## R.alt:     counts of allele 2 in this pool
    R.alt<-minor.reads.pool.ctl<-sum(minor.reads.person.ctl)
                      
    
    ##Considering error rates
    R.err[i,j]<-round(R*(1-e)+R.alt*e)
    R.alt.err[i,j]<-round(R.alt*(1-e)+R*e)
    
  }

  #case pools:
  for (j in (m.control+1):m){
    #parameters keep the same within each pool

     ##cases:
     #major allele reads 
    major.reads.person.cse<-0.01* rgamma(n,shape=(1-maf.cse[i])*shape[i,j],rate=rate[i,j])
   R<-major.reads.pool.cse<-sum(major.reads.person.cse)
                            
    #minor allele reads:
    minor.reads.person.cse<-0.01* rgamma(n,shape=maf.cse[i]*shape[i,j],rate=rate[i,j])
    R.alt<-minor.reads.pool.cse<-sum(minor.reads.person.cse)
                          
    ##Considering error rates
    R.err[i,j]<-round(R*(1-e)+R.alt*e)
    R.alt.err[i,j]<-round(R.alt*(1-e)+R*e)
    
  }
}

 return(list(R.err=R.err,
             R.alt.err=R.alt.err
              ))
}


  
