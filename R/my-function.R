#'Fitting promotion time cure model for clustered failure time data under Case-cohort study.
#'
#'This function fits a promotion time cure model for clustered failure time data under Case-cohort study, where the unknown baseline distribution function is
#'approximated by Bernstein polynomial. In this case-cohort design, the covariates can be considered as missing at random. This function employs the inverse probability
#'weighting to construct likelihood function.
#'
#'
#'@param K  numeric, the number of cluster.
#'@param n  numeric, the cluster size.
#'@param Data.1 a data.frame, the column names of the data.frame are Z.1(Covariate), Z.2(Covariate), T.C(Observation time), delta(Censor indicator),
#'Delta(Cure indicator, 1 if subject is not cure and 0 otherwise), HH.ij(Whether the covariates of the subject are observable).
#'@param beta numeric vector of 3 dimension, initial value of parameters.
#'@param gam numeric, the parameter of the Gumbel model.
#'@param a numeric, the parameter of covariate Z.2.
#'@param tau numeric, the cure threshold.
#'@param c0  numeric, the parameter of the distribution function for censoring time.
#'@param P.BerC numeric, Sampling clusters using Bernoulli sample with probability P.BerC.
#'@param P.BerS numeric, Sampling subjects from the selected clusters using Bernoulli sample with probability P.BerS.
#'@param N numeric, the degree of the Bernstein basis polynomials.
#'@param boot_n numeric, the replication of the weighted bootstrap procedure.
#'@return clustcureCaseC returns a list containing the following components:
#'@return \item{coefficients}{the estimators of coefficients.}
#'@return \item{se}{the standard errors of the coefficients estimators.}
#'@name clustcureCaseC
#'@examples
#'library(maxLik)
#'## true parameter: c(1,1,-2)
#'K<-400; n<-6
#' data(sampleData)
#' beta<-c(1,1,-2); gam<-0.2; a<-0.15; tau<-1; c0<-4
#' P.BerC<-0.5; P.BerS<-0.5; N<-4; boot_n<-3
#'## estimate parameters
#'# clustcureCaseC(K,n,sampleData,beta,gam,a,tau,c0,P.BerC,P.BerS,N,boot_n)
#' @export clustcureCaseC
clustcureCaseC<-function(K,n,Data.1,beta,gam,a,tau,c0,P.BerC,P.BerS,N,boot_n)
{
  beta.0<-beta[1];beta.1<-beta[2];beta.2<-beta[3]

  #Data.1<-gene.data1(K,n,beta.0,beta.1,beta.2,gam,a,tau,c0,P.BerC,P.BerS)
  Z.1<-Data.1$Z.1; Z.2<-Data.1$Z.2; T.C<-Data.1$T.C; delta<-Data.1$delta; HH.ij<-Data.1$HH.ij
  Delta<-Data.1$Delta

  omega<-HH.ij/(delta+(1-delta)*P.BerC*P.BerS)

  F_n<-function(para,TT)
  {
    F.n<-0
    for(k in 1:(N+1))
    {
      F.n<-F.n+para[k]*choose(N,k-1) * ((TT)/tau)^(k-1) * (1-TT/tau)^(N-k+1)
    }

    return(F.n)
  }

  f_n<-function(para,TT)
  {
    f.n<-0
    for(k in 1:N)
    {
      f.n<-f.n+(N/tau)*(para[k+1]-para[k])*choose(N-1,k-1) * (TT/tau)^(k-1) * (1-TT/tau)^(N-k)
    }

    return(f.n)
  }


  Loglike<-function(parabeta)
  {
    beta0<-parabeta[1];beta1<-parabeta[2];beta2<-parabeta[3]
    phi<-parabeta[4:(N+4)]
    F.n<-F_n(phi,T.C)
    f.n<-f_n(phi,T.C)

    Like<-sum(omega*(delta*(beta0+beta1*Z.1+beta2*Z.2)+
                       delta*log(f.n)-Delta*exp(beta0+beta1*Z.1+beta2*Z.2)*F.n
                     -(1-Delta)*exp(beta0+beta1*Z.1+beta2*Z.2)
    ))
    return(Like)
  }


  C1<-matrix(0,N+1,N+1)
  diag(C1)<-1

  for(j1 in 1:N)
  {
    C1[j1+1,j1]<--1
  }
  C1<-rbind(C1,c(rep(0,N),-1),c(rep(0,N),1))
  C<-cbind(matrix(0,N+3,3),C1)
  D<-c(rep(0,N+1),1.01,-0.99)

  re_esti<-maxLik::maxLik(Loglike, start = c(c(beta.0+runif(1,0,0.05),beta.1+runif(1,0,0.05),beta.2+runif(1,-0.05,0)),seq(0.05,1,length=N+1)), method = 'NM', constraints = list(ineqA = C, ineqB = D))$estimate

  Loglikeb<-function(parabeta)
  {
    beta0<-parabeta[1];beta1<-parabeta[2];beta2<-parabeta[3]
    phi<-parabeta[4:(N+4)]

    F.n<-F_n(phi,T.C)
    f.n<-f_n(phi,T.C)

    Like<-sum(e*omega*(delta*(beta0+beta1*Z.1+beta2*Z.2)+
                         delta*log(f.n)-Delta*exp(beta0+beta1*Z.1+beta2*Z.2)*F.n
                       -(1-Delta)*exp(beta0+beta1*Z.1+beta2*Z.2)
    ))
    return(Like)
  }


  Re_bt<-NULL;Re_btFullC<-NULL
  for(j3 in 1:boot_n)
  {
    e<-rexp(K*n,rate = 1)

    re_bt<-maxLik::maxLik(Loglikeb, start = c(c(beta.0+runif(1,0,0.05),beta.1+runif(1,0,0.05),beta.2+runif(1,-0.05,0)),seq(0.05,1,length=N+1)), method = 'NM', constraints = list(ineqA = C, ineqB = D))$estimate
    Re_bt<-rbind(Re_bt,re_bt)

  }

  se_case<-apply(Re_bt[,1:3],MARGIN = 2,sd)

  result<-list(coefficients=re_esti[1:3],se=se_case)

  return(result)
}


#' @title sample data
#' @description this is a sample data
#' @format A data frame with ... rows and ... variables, including:
#' \describe{
#' \item{id}{the clustering number}
#' \item{Z.1}{covariate 1 comes from b(1,0.5) distribution}
#' \item{Z.2}{covariate 2 comes from U(a,a+1) distribution}
#' \item{T.c}{the observed failure time}
#' \item{delta}{the censoring indicator}}
#' \item{Delta}{cured indicator}
#' \item{HH.ij}{the case-cohort indicator, 1 indicating the subject from the case-cohort sample}
#' }
#' @keywords datasets
"sampleData"

