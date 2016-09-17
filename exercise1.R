library(Matrix)
##########################################################################################
# Linear regression
##########################################################################################
# C)

invmethod=function(X,y,W)
{
  beta=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y
  return(beta)
}

mymethod=function(X,y,W) # LU decomposition
{
  A=t(X)%*%W%*%X
  b=t(X)%*%W%*%y
  ch=chol(A)
  D=diag(ch)
  L=t(ch/D)
  c=forwardsolve(L,b) # solve c
  d=solve(diag(D))%*%b 
  beta=backsolve(d,c) # solve for beta
  return(beta)
}

mymethod2=function(A,b) # LU decomposition
{
  L=chol(A)
  c=solve(t(L))%*%b
  d=solve(L)%*%c
  return(d)
}
##########################################################################################
#D)

N = 2000
P = 500
spar=0.05
X = matrix(rnorm(N*P), nrow=N)
mask1 = matrix(rbinom(N*P,1,spar), nrow=N)
X = mask1*X
y = matrix(rnorm(N*P), nrow=N)
mask2 = matrix(rbinom(N*P,1,spar), nrow=N)
y = mask1*y
W=diag(N)
#X[1:10, 1:10]
beta1=invmethod(X,y,W)
beta2=mymethod(X,y,W)
##########################################################################################
# Generalized linear models
##########################################################################################
# B)

# Functions:
# Calculate weights:
wts = function(X,B)
{
  1/(1+exp(-X%*%B))
}
# calculate differences:
dist=function(mob)
{
  sqrt(sum(mob^2))
}
# Gradient:
grad = function(y,X,w,m)
{
  -t(X)%*%(y-m*w)	
}
# Descent:
descent = function(beta0,y,X,alpha,tol,iter)
{
  m=dim(X)[1]
  n=dim(X)[2]
  mob=matrix(0,iter,n)
  mob[1,]=beta0
  logl=rep(0,iter)
  dist=rep(0,iter)
  gradient=rep(0,iter)
  mvec=rep(1,m)
  
  for(ii in 2:iter)
  {
    w=wts(X,mob[ii-1,])
    gradient[ii-1]=grad(y,X,w,mvec)
    mob[ii,]=mob[ii-1,]-alpha*gradient[ii-1]
    dist[ii]=sqrt(sum((mob[ii,]-mob[ii-1,])^2))
    if (dist[ii]<tol) 
    {break}
    logl[ii]=-sum(y*log(w+1e-6)+(1-y)*log(1-w+1e-6))
  }
  return(list(mob=mob,logl=logl,dist=dist))
}

# Input data:
wdbc = read.csv('C:/Users/Sputnik/Desktop/Statistical Models for Big Data/wdbc.csv', header=FALSE)
ya = as.character(wdbc[,2])
y = rep(0,length(ya))
y[which(ya=='M')] = 1
X = as.matrix(data[,2:11])
X = scale(X)
X = cbind(rep(1,length(ya)),X)


#glm
glm1 = glm(y~X-1, family='binomial')
coef=glm1$coefficients

#
beta0=matrix(rnorm(11),1)
fit=descent(beta0,y,X,alpha=1e-2,tol=1e-6,iter=10000)
fit

plot(fit$logl,type='l',log='xy')
plot(fit$dist,type='l',log='xy')

## The algorithm converges very fast.

# D) Newton's method:

hessian = function(X,m,w)
{
  t(X) %*% diag(m*(w+1e-6)*(1-w+1e-6)) %*% X
}

newton = function(y,X,B0,m=1,tol,iter,alpha)
{
  # defining relevant variables and Bmat
  
  m=dim(X)[1]
  n=dim(X)[2]
  mob=matrix(0,iter,n)
  mob[1,]=beta0
  dist=rep(0,iter)
  logl=rep(0,iter)
  mvec=rep(1,m)
  # iteration loop
  for(ii in 2:iter)
  {
    w=as.numeric(wts(X,mob[ii-1,]))
    Hess=hessian(X,mvec,w)
    Grad=-grad(y,X,w,mvec)
    delB=mymethod2(Hess,Grad)
    mob[ii,]=mob[ii-1,]+delB
    dist[ii]=sqrt(sum((mob[ii,]-mob[ii-1,])^2))
    if(dist[ii] <= tol){ break }
    logl[ii]=-sum(y*log(w+1e-6)+(1-y)*log(1-w+1e-6))
  }
  return(list(mob=mob,logl=logl,dist=dist))
}

fit2 = newton(y,X,B0,tol=1e-2,iter=1000,alpha=1)
plot(fit2$loglik,type='l',log='xy')
plot(fit2$dist,type='l',log='xy')


# Form the plot, we can see that it vibrates around a certain level and does not converge at the end.




