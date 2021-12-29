# Function for computing 
# rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties 
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance

smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    library(MASS)
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }

# Function for adding a propensity score caliper to a distance matrix dmat
# calipersd is the caliper in terms of standard deviation of the logit propensity scoe
addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

##### matching: does matching (a wrapper on top of optmatch...you should use this after understanding each step below)
### FUNCTION: creates a vector that matches treatment to control
### INPUT: 1) A: a n-length vector that lists treatment (as 1) and control (as 0)
###        2) X: a n-by-p matrix of covariates
###        3) type: type of matching, defaults to full matching
###        4) ratio: if full matching or 1:m matching, checks to see the control/treatment ratio
###        5) matchInfo: should the function print out the matching information?
### OUTPUT: 1) a n-vector vector containing match numbers.
### NOTE: The matching function (smahal, prop.score, caliper, and fullmatch) are
###       very sensitive to how the observations are ordered. Currentyl, smahal
###       orders the observations by how Z is ordered.
matching = function(A,X,type=c("full","pair","1:m"),ratio,matchInfo=FALSE) {
  type = match.arg(type)
  
  # Pre-processing
  Z = as.logical(Z);
  
  # Rank-based distance matrix
  rank.mat = smahal(Z,X)
  
  # Propensity scores
  model = glm(Z ~ X,family="binomial"); prop.score.logit = predict(model);
  
  # Fix rank-based distance matrix by propensity score calipers
  rank.mat.caliper = addcaliper(rank.mat,Z,prop.score.logit)
  
  # Actual matching
  if(type == "pair") pairmatchvec = pairmatch(rank.mat.caliper,data=1:length(Z))
  if(type == "full") {
    if(missing(ratio)) pairmatchvec = fullmatch(rank.mat.caliper,data=1:length(Z))
    else pairmatchvec = fullmatch(rank.mat.caliper,data=1:length(Z),min.controls=1/ratio,max.controls=ratio)
  }
  if(type == "1:m") pairmatchvec = pairmatch(rank.mat.caliper,data=1:length(Z),controls=ratio)
  if(matchInfo) print(summary(pairmatchvec))
  
  # Post-processing
  pairmatchvec_numeric = as.numeric(substr(pairmatchvec,start=3,stop=20))  
  return(pairmatchvec_numeric) 
}