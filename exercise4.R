rm(list=ls())

library(microbenchmark)
library(Rcpp)
library(RcppEigen)
library(Matrix)
library(permute)

sourceCpp(file= "C:/Users/weizhang/Desktop/R/improvingSGD.cpp")

# First, we try small sample data
wdbc = read.csv('C:/Users/weizhang/Desktop/R/wdbc.csv', header=FALSE)
y = wdbc[,2]
Y = rep(0,length(y)); Y[y=='M']=1
X = as.matrix(wdbc[,-c(1,2)])
scrub = which(1:ncol(X) %% 3 == 0)
scrub = 11:30
X = X[,-scrub]
X = scale(X) #Normalize design matrix features.
X = cbind(rep(1,nrow(X)),X)
X = Matrix(X,sparse=T)
Xt = t(X)
m = rep(1,nrow(X))	

#without penalty
testset1 = improvingSGD(Xt,Y,m,step=5,beta0=rep(0,11),lambda=0,npass=5000)
testset1$beta

#[1]  0.3394583 -3.9362003  1.6120498 -4.8625655 13.7152620  1.0397486  0.1588215  0.7129504  2.6758299  0.4568794
#[11] -0.5093412

#with penalty
testset2 = improvingSGD(Xt,Y,m,step=5,beta0=rep(0,11),lambda=.1,npass=5000)
testset2$beta

#[1] -1.259868e+68 -2.116988e-05 -2.523405e-05 -2.296332e-05 -3.845077e-05  5.146706e+84 -3.454989e-06 -2.106173e-05
#[9] -4.594802e-05 -1.586807e-06 -3.025312e+40

#We can see the algorithm works well with the small sample, and it is also fast.

# Now we work on the full size data

Xt_URL = readRDS(file='C:/Users/weizhang/Desktop/R/url_tX.rds')
Y_URL = readRDS(file='C:/Users/weizhang/Desktop/R/url_y.rds')
m_URL <- rep(1,ncol(Xt_URL))

#without penalty
fullsize1 = improvingSGD(Xt_URL,Y_URL,m_URL,step=1,beta0=rep(0,nrow(Xt_URL)),lambda=0,npass=5) #took 59 sec

#With penalty
fullsize2 = improvingSGD(Xt_URL,Y_URL,m_URL,step=1,beta0=rep(0,nrow(Xt_URL)),lambda=.1,npass=5) #took 103 sec
