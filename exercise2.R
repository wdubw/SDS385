# Online learning
# Input data
wdbc = read.csv('C:/Users/Sputnik/Desktop/Statistical Models for Big Data/wdbc.csv', header=FALSE, row.names = 1)
ya = as.character(wdbc[,1])
y = rep(0,length(ya))
y[which(ya=='M')] = 1
X = as.matrix(data[,2:11])
X = scale(X)
X = cbind(rep(1,length(ya)),X)

#####################################################################
# C)
# calculate difference
dist=function(mob)
{
  sqrt(sum(mob^2))
}
# calculate the gradient
gradnorm=function(y,X,B)
{
  -sum((y-X*B)*X)
}
# calculate the loglike
logliknorm=function(y,X,B)
{
  -sum(dnorm(y,X*B,sd=0.5,log=TRUE))
}

sgd=function(y,X,beta0,m=1,iter,alpha)
{
  p=1
  N=500
  mob=matrix(0,iter,p)
  mob[1,]=beta0
  loglik=rep(0,iter)
  distance=rep(0,iter)
  
  for(ii in 2:iter)
  {
    # randomize
    ind=sample(1:N,1)
    yran=y[ind]
    Xran=X[ind]
    mob[ii,]=mob[ii-1,]-alpha*gradnorm(yran,Xran,mob[ii-1,])
    distance[ii] = dist(mob[ii,]-mob[ii-1,])
    loglik[ii] = logliknorm(y,X,mob[ii,])
  }
  return(list(mob=mob,loglik=loglik,dist=distance))
}

beta0=-10
m=1
iter=1000
alpha=0.1
fit2=sgd(y,X,beta0,m,iter,alpha)
fit2$dist
plot(fit2$dist[1:iter],type='l')
abline(h=1,lty=2,col=2)
plot(fit2$loglik[1:iter],type='l')

#####################################################################
# D)

sgd2=function(y,X,beta0,m=1,iter,alpha,C)
{
  p=1
  N=500
  mob=matrix(0,iter,p)
  mob[1,]=beta0
  loglik=rep(0,iter)
  distance=rep(0,iter)
  
  for(ii in 2:iter)
  {
    # randomize
    ind=sample(1:N,1)
    yran=y[ind]
    Xran=X[ind]
    CC=C
    mob[ii,]=mob[ii-1,]-stepsizeRM(CC,ii,alpha)*gradnorm(yran,Xran,mob[ii-1,])
    distance[ii] = dist(mob[ii,]-mob[ii-1,])
    loglik[ii] = logliknorm(y,X,mob[ii,])
  }
  return(list(mob=mob,loglik=loglik,dist=distance))
}

beta0=-2
m=1
iter=200
alpha=1
C=1

fit3=sgd2(y,X,beta0,m,iter,alpha,C)
fit3$dist
plot(fit3$dist[1:iter],type='l')
abline(h=1,lty=2,col=2)
plot(fit3$loglik[1:iter],type='l')

#####################################################################
# E)
sgd3=function(y,X,beta0,m=1,iter,alpha)
{
  p=1
  N=500
  mob=matrix(0,iter,p)
  mob[1,]=beta0
  loglik=rep(0,iter)
  distance=rep(0,iter)
  sumb=matrix(0,iter,p)
  sumb[1,]=mob[1,]
  for(ii in 2:iter)
  {
    # randomize
    ind=sample(1:N,1)
    yran=y[ind]
    Xran=X[ind]
    mob[ii,]=mob[ii-1,]-alpha*gradnorm(yran,Xran,mob[ii-1,])
    sumb[ii,]=(sumb[ii-1,]*(ii-1)+mob[ii,])/ii # Update the average
    distance[ii] = dist(sumb[ii,]-sumb[ii-1,])
    loglik[ii] = logliknorm(y,X,mob[ii,])
  }
  return(list(mob=mob,loglik=loglik,dist=distance))
}
beta0=-10
m=1
iter=1000
alpha=0.1
fit4=sgd3(y,X,beta0,m,iter,alpha)
fit4$dist
plot(fit4$dist[1:iter],type='l')
abline(h=1,lty=2,col=2)
plot(fit4$loglik[1:iter],type='l')