library('tbd')

expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

SACE.est <- function(Z,S,Y,X,V,rho=1,need.variance=FALSE){
  N <- length(Z)
  if (is.null(X)) X=as.matrix(rep(1,N))
  SE = NULL
  SEu = NULL
  SE.sc = NULL
  SE.wzr = NULL
  m.S <- glm(S~., family='binomial', data=data.frame(S,Z,X,V))
  m.Y <- glm(Y~., family='binomial', data=data.frame(Y,Z,X,V)[S==1,])
  fit.S.V1 <- predict(m.S, newdata=data.frame(S,Z=1,X,V=1), type='response')
  fit.S.V0 <- predict(m.S, newdata=data.frame(S,Z=1,X,V=0), type='response')
  fit.Y1.V1 <- predict(m.Y, newdata=data.frame(Y,Z=1,X,V=1), type='response')
  fit.Y1.V0 <- predict(m.Y, newdata=data.frame(Y,Z=1,X,V=0), type='response')
  fit.Y1 <- predict(m.Y, newdata=data.frame(Y,Z=1,X,V), type='response')
  fit.Y0 <- predict(m.Y, newdata=data.frame(Y,Z=0,X,V), type='response')
  mu1X <- (rho*fit.S.V1*fit.Y1.V1*(1-fit.S.V0) - fit.S.V0*fit.Y1.V0*(1-fit.S.V1))/
    (rho*fit.S.V1*(1-fit.S.V0) - fit.S.V0*(1-fit.S.V1))
  mu0X <- fit.Y0
  mu1 <- mean(mu1X)
  mu0 <- mean(mu0X)
  mu1.sc <- mean(fit.Y1)
  mu0.sc <- mu0
  RD <- mu1-mu0
  RD.sc <- mu1.sc-mu0.sc
  fit.wzr <- sace(Z,S,Y,X,V,need.variance=FALSE)
  mu1.wzr <- fit.wzr$mu_1_LL
  mu0.wzr <- fit.wzr$mu_0_LL
  RD.wzr <- fit.wzr$sace
  
  if (need.variance==TRUE){
    m.V <- glm(V~., family='binomial', data=data.frame(V,X))
    m.Z <- glm(Z~., family='binomial', data=data.frame(Z,V,X))
    fit.V <- predict(m.V, type='response')
    fit.Z <- predict(m.Z, type='response')
    fit.S0 <- predict(m.S, newdata=data.frame(S,Z=0,X,V), type='response')
    p1 = p0111 = (1-fit.V)*fit.Z*fit.S.V0*fit.Y1.V0
    p2 = p0110 = (1-fit.V)*fit.Z*fit.S.V0*(1-fit.Y1.V0)
    p3 = p1111 = fit.V*fit.Z*fit.S.V1*fit.Y1.V1
    p4 = p1110 = fit.V*fit.Z*fit.S.V1*(1-fit.Y1.V1)
    p5 = p010s = (1-fit.V)*fit.Z*(1-fit.S.V0)
    p6 = p110s = fit.V*fit.Z*(1-fit.S.V1)
    p7 = ps011 = (1-fit.Z)*fit.S0*fit.Y0
    p8 = ps010 = (1-fit.Z)*fit.S0*(1-fit.Y0)
    p12 = p011s = p1+p2
    p34 = p111s = p3+p4
    p56 = ps10s = p5+p6
    p78 = ps01s = p7+p8
    ps11s = p12+p34
    VARX = ((p010s^4*p1111*p1110*p111s
             +p010s*p110s*ps10s*(p0110*p1111-p1110*p0111)^2
             -2*p010s^3*p1111*p1110*p011s*p110s
             +p010s^2*p110s^2*(p0110*p1111*(p0110+p1111)
                               +p1110*p0111*(p1110+p0111))
             -2*p110s^3*p0111*p0110*p111s*p010s
             +p110s^4*p0111*p0110*p011s)/
              (p111s*p010s-p011s*p110s)^4) + (ps011*ps010/ps01s^3)
    plaus = (VARX<1)
    SE = sqrt(VARX[plaus]/sum(plaus))
    SEu = sqrt(VARX[plaus]/sum(plaus)+var((mu1X-mu0X)[plaus]))
  }
  
  return(list(Z=Z,S=S,Y=Y,X=X,V=V,mu1=mu1,mu0=mu0,sace=RD,
              mu1.sc=mu1.sc,mu0.sc=mu0,sace.sc=RD.sc,
              mu1.wzr=mu1.wzr,mu0.wzr=mu0.wzr,sace.wzr=RD.wzr,
              se=SE,seu=SEu))
}





## SORE, No covariates

sc <- function(Z,S,Y,V,need.variance=FALSE,alpha=0.05){
  N <- length(Z)
  p111 = mean(Z==1&S==1&Y==1)
  p110 = mean(Z==1&S==1&Y==0)
  p11s = p111 + p110
  p011 = mean(Z==0&S==1&Y==1)
  p010 = mean(Z==0&S==1&Y==0)
  p01s = p011 + p010
  mu1 = p111/p11s
  mu0 = p011/p01s
  RD = mu1 - mu0
  if (RD>1) RD=1
  if (-1>RD) RD=-1
  if (need.variance==TRUE){
    za = -qnorm(alpha/2)
    var1 = mean(p111*p110/p11s^3)/N
    var0 = mean(p011*p010/p01s^3)/N
    SE = sqrt(var1+var0)
    ci_l0 = expit(logit(mu0)-za*sqrt(mean(p01s/(p010*p011*N))))
    ci_u0 = expit(logit(mu0)+za*sqrt(mean(p01s/(p010*p011*N))))
    ci_l1 = expit(logit(mu1)-za*sqrt(mean(p11s/(p110*p111*N))))
    ci_u1 = expit(logit(mu1)+za*sqrt(mean(p11s/(p110*p111*N))))
    ci_l = RD - sqrt((ci_l1-mu1)^2+(ci_u0-mu0)^2)
    ci_u = RD + sqrt((ci_u1-mu1)^2+(ci_l0-mu0)^2)
    return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD,
                se=SE,ci_l=ci_l,ci_u=ci_u))
  }
  return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD))
}

sore <- function(Z,S,Y,V,need.variance=FALSE,nboot=0,alpha=0.05){
  N <- length(Z)
  p1 = p0111 = mean(V==0&Z==1&S==1&Y==1)
  p2 = p0110 = mean(V==0&Z==1&S==1&Y==0)
  p3 = p1111 = mean(V==1&Z==1&S==1&Y==1)
  p4 = p1110 = mean(V==1&Z==1&S==1&Y==0)
  p5 = p010s = mean(V==0&Z==1&S==0)
  p6 = p110s = mean(V==1&Z==1&S==0)
  p7 = ps011 = mean(Z==0&S==1&Y==1)
  p8 = ps010 = mean(Z==0&S==1&Y==0)
  ps010 = mean(Z==0&S==1&Y==0)
  p111s = p1111 + p1110
  ps111 = p1111 + p0111
  ps110 = p1110 + p0110
  p011s = p0111 + p0110
  ps01s = ps011 + ps010
  ps10s = p110s + p010s
  ps11s = p111s + p011s
  mu1 = (p0111*p110s-p1111*p010s)/(p011s*p110s-p111s*p010s)
  mu0 = ps011/ps01s
  if (!is.na(mu1)){
    if (mu1>1) mu1=1
    if (mu1<0) mu1=0
  }
  RD = mu1 - mu0
  singularity = FALSE
  if (!is.na(RD)){
    if (RD>1) RD=1
    if (-1>RD) RD=-1
  }
  if (is.na(RD)){
    mu1 = ps111/ps11s
    RD = mu1-mu0
    if (is.infinite(RD)) RD=0
    singularity = TRUE
  }
  if (need.variance==TRUE){
    if (singularity==TRUE){
      za = -qnorm(alpha/2)
      var1 = mean(ps111*ps110/ps11s^3)/N
      var0 = mean(ps011*ps010/ps01s^3)/N
      SE = sqrt(var1+var0)
      ci_l0 = expit(logit(mu0)-za*sqrt(mean(ps01s/(ps010*ps011*N))))
      ci_u0 = expit(logit(mu0)+za*sqrt(mean(ps01s/(ps010*ps011*N))))
      ci_l1 = expit(logit(mu1)-za*sqrt(mean(ps11s/(ps110*ps111*N))))
      ci_u1 = expit(logit(mu1)+za*sqrt(mean(ps11s/(ps110*ps111*N))))
    } else{
    za = -qnorm(alpha/2)
    vu = (p010s^4*p1111*p1110*p111s
          +p010s*p110s*ps10s*(p0110*p1111-p1110*p0111)^2
          -2*p010s^3*p1111*p1110*p011s*p110s
          +p010s^2*p110s^2*(p0110*p1111*(p0110+p1111)
                            +p1110*p0111*(p1110+p0111))
          -2*p110s^3*p0111*p0110*p111s*p010s
          +p110s^4*p0111*p0110*p011s)
    var1 = mean(vu/(p111s*p010s-p011s*p110s)^4)/N
    var0 = mean(ps011*ps010/ps01s^3)/N
    SE = sqrt(var1+var0)
    ci_l0 = expit(logit(mu0)-za*sqrt(mean(ps01s/(ps010*ps011*N))))
    ci_u0 = expit(logit(mu0)+za*sqrt(mean(ps01s/(ps010*ps011*N))))
    ci_l1 = expit(logit(mu1)-za*
                  sqrt(mean(vu/(p3*p5-p1*p6)^2/(p4*p5-p2*p6)^2)/N))
    ci_u1 = expit(logit(mu1)+za*
                  sqrt(mean(vu/(p3*p5-p1*p6)^2/(p4*p5-p2*p6)^2)/N))
    }
    if (is.na(SE)) SE=1
    if (is.infinite(SE)) SE=1
    if (SE>1) SE=1
    if (is.na(ci_l0)) ci_l0=0
    if (is.na(ci_l1)) ci_l1=0
    if (is.na(ci_u0)) ci_u0=1
    if (is.na(ci_u1)) ci_u1=1
    if (is.infinite(ci_l0)) ci_l0=0
    if (is.infinite(ci_l1)) ci_l1=0
    if (is.infinite(ci_u0)) ci_u0=1
    if (is.infinite(ci_u1)) ci_u1=1
    if (ci_l0<0) ci_l0=0
    if (ci_l1<0) ci_l1=0
    if (1<ci_u0) ci_u0=1
    if (1<ci_u1) ci_u1=1
    ci_l = RD - sqrt((ci_l1-mu1)^2+(ci_u0-mu0)^2)
    ci_u = RD + sqrt((ci_u1-mu1)^2+(ci_l0-mu0)^2)
    if (nboot!=0){
      bt = rep(NA, nboot)
      for (b in 1:nboot){
      s = sample(N,replace=TRUE)
      Zs = Z[s]
      Ss = S[s]
      Ys = Y[s]
      Vs = V[s]
      bt[b] = sore(Zs,Ss,Ys,Vs,need.variance=FALSE,nboot=FALSE,alpha=0.05)$sace
      }
      bt = sort(na.omit(bt))
      nboot = length(bt)
      bt_l = bt[alpha/2*nboot]
      bt_u = bt[(1-alpha/2)*nboot+1]
      return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD,
                  se=SE,ci_l=ci_l,ci_u=ci_u,bt_l=bt_l,bt_u=bt_u))
    }
    return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD,
                se=SE,ci_l=ci_l,ci_u=ci_u))
  }
  return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD))
}

wzr <- function(Z,S,Y,V,need.variance=FALSE,alpha=0.05){
  Y[is.na(Y)] = 0
  YS = Y*S
  PLL0 = mean(S[Z==0&V==0])
  PLL1 = mean(S[Z==0&V==1])
  PS0 = mean(S[Z==1&V==0])
  PS1 = mean(S[Z==1&V==1])
  EYS0 = mean(YS[Z==1&V==0])
  EYS1 = mean(YS[Z==1&V==1])
  mu0 = mean(Y[Z==0&S==1])
  mu1 = (EYS1*(PLL0-PS0)-EYS0*(PLL1-PS1))/(PLL0*PS1-PLL1*PS0)
  RD = mu1 - mu0
  if (RD>1) RD=1
  if (-1>RD) RD=-1
  return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD))
}

sace <- function(Z,S,Y,V){
  sace.sc = sc(Z,S,Y,V)$sace
  sace.sore = sore(Z,S,Y,V)$sace
  sace.wzr = wzr(Z,S,Y,V)$sace
  return(list(sace.sc=sace.sc,sace.sore=sace.sore,sace.wzr=sace.wzr))
}
