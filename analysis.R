dat = read.csv('leukemiaPKU.csv')
Z = 1-dat$TRANSPLANT
V = as.numeric(dat$AGE>27.15)
S = as.numeric(dat$TRM==0|dat$TRMT>730)
Y = 1-(dat$RELAPSE==0|dat$RELAPSET>730)
M = dat$MRD

sc(Z, S, Y, V, need.variance=TRUE, nboot=1000)
sore(Z, S, Y, V, need.variance=TRUE, nboot=1000)

sc(Z[M==0], S[M==0], Y[M==0], V[M==0], need.variance=TRUE, nboot=1000)
sore(Z[M==0], S[M==0], Y[M==0], V[M==0], need.variance=TRUE, nboot=1000)

sc(Z[M==1], S[M==1], Y[M==1], V[M==1], need.variance=TRUE, nboot=1000)
sore(Z[M==1], S[M==1], Y[M==1], V[M==1], need.variance=TRUE, nboot=1000)
