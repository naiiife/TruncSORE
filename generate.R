generate <- function(size,a1,a2,b1,b2,c1,c2,seeds,d1=0,d0=0){
  set.seed(seeds)
  
  G1 <- rbinom(size, 1, 0.5)
  G2 = (1-G1) * rbinom(size, 1, 0.6)
  S0 = G1
  S1 = G1 + G2
  
  V <- rbinom(size, 1, 0.3 + a1*G1 + a2*G2)
  
  Z = rbinom(size, 1, 0.3 + b1*G1 + b2*G2 - 0.1*V)
  S = Z*S1 + (1-Z)*S0
  
  Y1 = rbinom(size, 1, 0.3 + c1*G1 + c2*G2 + d1*V)
  Y0 = rbinom(size, 1, 0.1 + c1*G1 + c2*G2 + d0*V)
  Y = Z*Y1 + (1-Z)*Y0
  Y[S==0] = NA
  TrueRD = mean((Y1-Y0)[G1==1])
  
  return(list(Z=Z,S=S,Y=Y,V=V,S1=S1,S0=S0,Y1=Y1,Y0=Y0,TrueRD=TrueRD))
}
