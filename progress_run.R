library(ggplot2)

iter <- 1000
sizel = c(200,800,3000)
# size = 800

# setting1 WZR andSORE correct
# setting2 SC and SORE correct
# setting3 SC, WZR, SORE correct
# setting4 SORE correct
# setting5 SORE correct, and sensitivity
# SC correct <=> c1=c2
# WZR correct <=> b1=b2=0
# SORE correct <=> a2=b2=0

a1l <- c(.4, .4, .4, .4, .4)  #strength of V
a2l <- c(.0, .0, .0, .0, .0)  #sensitivity
b1l <- c(.0, .4, .0, .4, .4)  #extent of nonrandomization
b2l <- c(.0, .0, .0, .0, .2)  #sensitivity
c1l <- c(.4, .4, .4, .4, .4)
c2l <- c(.2, .4, .4, .2, .2)
# 
# a1l <- c(5,5,5,5,5)  #strength of V
# a2l <- c(0,0,0,0,0)  #sensitivity
# b1l <- c(0,0,0,0,0)  #extent of nonrandomization
# b2l <- c(0,0,0,0,0)  #sensitivity
# c1l <- c(4,3.5,3,2.5,2)
# c2l <- c(4,3.5,3,2.5,2)

for (size in sizel){
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
      dat = generate(size,a1,a2,b1,b2,c1,c2,seeds=setting+i)
      #TrueRD = dat$TrueRD
      TrueRD = 0.2
      Z = dat$Z
      S = dat$S
      Y = dat$Y
      #X = dat$X
      V = dat$V
      #res = SACE.est(Z,S,Y,X[,-1],V)
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
  # bp = ggplot(restable, aes(x=Setting, y=Bias, fill=Method)) + 
  #   geom_boxplot() + coord_cartesian(ylim = c(-0.25, 0.4)) + geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
  bp = ggplot(restable, aes(x=Setting, y=Bias, fill=Method)) + 
    geom_boxplot() + coord_cartesian(ylim = c(-0.5, 0.5)) + geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
  bp
  # if you need grey:
  #bp + scale_fill_grey(start=0.2, end=0.8) + theme_classic()
  
}
