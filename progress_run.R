library(ggplot2)

iter <- 1000
sizel = c(200,800,3000)

# setting1 WZR and SORE correct
# setting2 SC and SORE correct
# setting3 SC,WZR,SORE correct
# setting4 SORE correct
# setting5 SORE correct (remark)
# SC correct <= c1=c2
# WZR correct <= b1=b2=0
# SORE correct <= a2=b2=0

a1l <- c(.4, .4, .4, .4, .4)  #strength of V
a2l <- c(.0, .0, .0, .0, .0)  #sensitivity
b1l <- c(.0, .4, .0, .4, .4)  #extent of nonrandomization
b2l <- c(.0, .0, .0, .0, .2)  #sensitivity
c1l <- c(.4, .4, .4, .4, .4)
c2l <- c(.2, .4, .4, .2, .2)
d1 = d0 = 0

## In sensitivity analysis I, change
## (1) a2l = c(.4, .4, .4, .4, .4)
## (2) a2l = c(.2, .2, .2, .2, .2)
## (3) a1l = c(.0, .0, .0, .0, .0)
## In sensitivity analysis II, change
## (1) d1 = d0 = 0.1
## (2) d1 = 0.1; d0 = 0
## (3) d1 = 0; d0 = 0.1

#for (size in sizel){
#size = sizel[1]
size = sizel[2]
#size = sizel[3]
  restable <- data.frame(Bias=double(),
                         Setting=integer(),
                         Method=character(),stringsAsFactors=FALSE)
  
  for (setting in 1:length(a1l)){
    a1 = a1l[setting]
    a2 = a2l[setting]
    b1 = b1l[setting]
    b2 = b2l[setting]
    c1 = c1l[setting]
    c2 = c2l[setting]
    
    for (i in 1:iter){
      dat = generate(size,a1,a2,b1,b2,c1,c2,seeds=setting+i,d1,d0)
      #TrueRD = dat$TrueRD
      TrueRD = 0.2
      Z = dat$Z
      S = dat$S
      Y = dat$Y
      #X = dat$X
      V = dat$V
      res = sace(Z,S,Y,V)
      
      restable[nrow(restable)+1,] <- c(res$sace.sc - TrueRD, setting, 'SC')
      restable[nrow(restable)+1,] <- c(res$sace.wzr - TrueRD, setting, 'WZR')
      restable[nrow(restable)+1,] <- c(res$sace.sore - TrueRD, setting, 'SORE')
      }
  }
  # Change bar order grouped box-plot (fill-variable), prepare for plot
  restable <- transform(restable, Bias = as.numeric(Bias), 
                Setting = as.factor(Setting), Method = ordered(as.factor(Method), levels = c('SC', 'WZR', 'SORE')))
  restable <- na.omit(restable)

  # boxplot + ylim + yintercept
  bp = ggplot(restable, aes(x=Setting, y=Bias, fill=Method)) + 
    geom_boxplot() + coord_cartesian(ylim = c(-0.5, 0.5)) + geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
  bp
  # if you need grey:
  #bp + scale_fill_grey(start=0.2, end=0.8) + theme_classic()
  
#}
