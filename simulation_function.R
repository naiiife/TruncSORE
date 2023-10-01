expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

sc <- function(Z,S,Y,V,need.variance=FALSE,nboot=0,alpha=0.05){
  N = length(Z)
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
    if (nboot!=0){
      bt = rep(NA, nboot)
      for (b in 1:nboot){
        s = sample(N,replace=TRUE)
        Zs = Z[s]
        Ss = S[s]
        Ys = Y[s]
        Vs = V[s]
        bt[b] = sc(Zs,Ss,Ys,Vs,need.variance=FALSE,nboot=0,alpha=0.05)$sace
      }
      bt = sort(na.omit(bt))
      nboot = length(bt)
      bt_l = bt[alpha/2*nboot]
      bt_u = bt[(1-alpha/2)*nboot+1]
      return(list(mu1=mu1,mu0=mu0,sace=RD,
                  se=SE,ci_l=ci_l,ci_u=ci_u,bt_l=bt_l,bt_u=bt_u))
    }
    return(list(mu1=mu1,mu0=mu0,sace=RD,
                se=SE,ci_l=ci_l,ci_u=ci_u))
  }
  return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD))
}

sore <- function(Z,S,Y,V,need.variance=FALSE,nboot=0,alpha=0.05){
  N = length(Z)
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
      bt[b] = sore(Zs,Ss,Ys,Vs,need.variance=FALSE,nboot=0,alpha=0.05)$sace
      }
      bt = sort(na.omit(bt))
      nboot = length(bt)
      bt_l = bt[alpha/2*nboot]
      bt_u = bt[(1-alpha/2)*nboot+1]
      return(list(mu1=mu1,mu0=mu0,sace=RD,
                  se=SE,ci_l=ci_l,ci_u=ci_u,bt_l=bt_l,bt_u=bt_u))
    }
    return(list(mu1=mu1,mu0=mu0,sace=RD,
                se=SE,ci_l=ci_l,ci_u=ci_u))
  }
  return(list(Z=Z,S=S,Y=Y,V=V,mu1=mu1,mu0=mu0,sace=RD))
}

wzr <- function(Z,S,Y,V,need.variance=FALSE,nboot=0,alpha=0.05){
  N = length(Z)
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
  if (nboot!=0){
    bt = rep(NA, nboot)
    for (b in 1:nboot){
      s = sample(N,replace=TRUE)
      Zs = Z[s]
      Ss = S[s]
      Ys = Y[s]
      Vs = V[s]
      bt[b] = wzr(Zs,Ss,Ys,Vs,need.variance=FALSE,nboot=FALSE,alpha=0.05)$sace
    }
    bt = sort(na.omit(bt))
    nboot = length(bt)
    bt_l = bt[alpha/2*nboot]
    bt_u = bt[(1-alpha/2)*nboot+1]
    return(list(mu1=mu1,mu0=mu0,sace=RD,bt_l=bt_l,bt_u=bt_u))
  }
  return(list(mu1=mu1,mu0=mu0,sace=RD))
}

sace <- function(Z,S,Y,V){
  sace.sc = sc(Z,S,Y,V)$sace
  sace.sore = sore(Z,S,Y,V)$sace
  sace.wzr = wzr(Z,S,Y,V)$sace
  return(list(sace.sc=sace.sc,sace.sore=sace.sore,sace.wzr=sace.wzr))
}
