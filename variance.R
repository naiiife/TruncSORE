
iter <- 10000
sizel = c(200,800,3000)

a1l <- c(.4, .4, .4, .4, .4)  #strength of V
a2l <- c(.0, .0, .0, .0, .0)  #sensitivity
b1l <- c(.0, .4, .0, .4, .4)  #extent of nonrandomization
b2l <- c(.0, .0, .0, .0, .2)  #sensitivity
c1l <- c(.4, .4, .4, .4, .4)
c2l <- c(.2, .4, .4, .2, .2)
d1 = d0 = 0

CV = NULL
WD = NULL
for (size in sizel){
  for (setting in 1:length(a1l)){
    a1 = a1l[setting]
    a2 = a2l[setting]
    b1 = b1l[setting]
    b2 = b2l[setting]
    c1 = c1l[setting]
    c2 = c2l[setting]
    
    cover.z = cover.t = cover.b = 0
    width.z = width.t = width.b = rep(NA,iter)
    for (i in 1:iter){
      dat = generate(size,a1,a2,b1,b2,c1,c2,seeds=setting*10000+i,d1,d0)
      #TrueRD = dat$TrueRD
      TrueRD = 0.2
      Z = dat$Z
      S = dat$S
      Y = dat$Y
      V = dat$V
      res = sore(Z,S,Y,V,need.variance=TRUE,nboot=1000,alpha=0.05)
      if (res$sace-1.96*res$se<=TrueRD & res$sace+1.96*res$se>=TrueRD) cover.z = cover.z+1
      if (res$ci_l<=TrueRD & res$ci_u>=TrueRD) cover.t = cover.t+1
      if (res$bt_l<=TrueRD & res$bt_u>=TrueRD) cover.b = cover.b+1
      width.z[i] = min(1,res$sace+1.96*res$se)-max(-1,res$sace-1.96*res$se)
      width.t[i] = res$ci_u-res$ci_l
      width.b[i] = res$bt_u-res$bt_l
    }
    CVnew = c(cover.z/iter, cover.t/iter, cover.b/iter)
    WDnew = c(mean(na.omit(width.z)), mean(na.omit(width.t)), mean(na.omit(width.b)))
    CV = rbind(CV, CVnew)
    WD = rbind(WD, WDnew)
    print(CV)
    print(WD)
  }
}

print(CV)
print(WD)
rbind(CV,WD)
