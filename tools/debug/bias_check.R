# Setup------------------------------------------------------------------------
library(searchPTA)
# DGP functions----------------------------------------------------------------
gamma_fun <- function(xvals,g){
  if (length(g)==1) g <- rep(g,nrow(xvals))
  value <- g*(2*(xvals[,1]-0.25)^2 + xvals[,2]) + (1-g)*(2*(xvals[,2]-0.5)^2 + xvals[,1])
  return(value)
}

fx <- function(nsamp,g,x1.shape1.g1,x1.shape2.g1,x1.shape1.g0,x1.shape2.g0,x2.shape1.g1,x2.shape2.g1,x2.shape1.g0,x2.shape2.g0){
  if (length(g)==1) g <- rep(g,nsamp)
  x1 <- g*rbeta(nsamp,x1.shape1.g1,x1.shape2.g1) + (1-g)*rbeta(nsamp,x1.shape1.g0,x1.shape2.g0)
  x2 <- g*rbeta(nsamp,x2.shape1.g1,x2.shape2.g1) + (1-g)*rbeta(nsamp,x2.shape1.g0,x2.shape2.g0)
  x <- cbind(x1,x2)
  return(x)
}
# Sample data------------------------------------------------------------------
n <- 20000
g <- sample(c(1,0),size=n,replace=TRUE,prob=c(0.5,0.5))
x1.shape1.g1 <- 1
x1.shape2.g1 <- 1
x1.shape1.g0 <- 1
x1.shape2.g0 <- 1
x2.shape1.g1 <- 1
x2.shape2.g1 <- 1
x2.shape1.g0 <- 1
x2.shape2.g0 <- 1
x <- fx(n,g,x1.shape1.g1,x1.shape2.g1,x1.shape1.g0,x1.shape2.g0,x2.shape1.g1,x2.shape2.g1,x2.shape1.g0,x2.shape2.g0)
X <- data.frame(x = x,g = g)
gamma1 <- gamma_fun(x,1)
gamma0 <- gamma_fun(x,0)
x1 <- subset(X, g==1, select=1:2)
x0 <- subset(X, g==0, select=1:2)
# Search for PTA regions-------------------------------------------------------
raf <- searchPTA::searchPTA(x=X,gamma1=gamma1,gamma0=gamma0,bta1=gamma1,beta0=gamma0,epsilon = 0.05,control = list(minbucket=50,maxdepth = 10,cp = 0))
preds <- raf$beta.diff
preds[is.na(preds)] <- -1 ## Substitute NAs for unidentified points by -1; necessary because of the behavior of "unique" function
# Plot bias--------------------------------------------------------------------
bias.per.region <- sort(unique(preds)) ## Get bias per PTA region and sort in ascending order; used to color by bias size
## Point colors
r <- ifelse(bias.per.region==-1,1,0)
g <- seq(0,1,length.out=length(bias.per.region)+1)[-1]
b <- ifelse(bias.per.region==-1,0,1)
##
par(mfrow=c(2,2))
## Plot per point
# cols <- seq(0,1,length.out=n)
# bias <- gamma1-gamma0
# plot(X[,1],X[,2],type="p",pch=15,xlim=c(0,1),ylim=c(0,1),xlab=bquote(X[1]),ylab=bquote(X[2]),col=rgb(cols[order(bias)],cols[order(bias)],cols[order(bias)]),main="Pointwise bias")

# x1 <- x2 <- seq(0,1,length.out=1000)
# X.plot <- expand.grid(x1,x2)
# g1 <- gamma_fun(X.plot,1)
# g0 <- gamma_fun(X.plot,0)
# bias <- matrix(g1-g0, nrow=length(unique(X.plot[,1])), ncol=length(unique(X.plot[,2])))
# cols <- seq(0,1,length.out=nrow(bias))
#
# fields::image.plot(x=sort(unique(X.plot[,1])),y=sort(unique(X.plot[,2])),z=bias,xlab=bquote(X[1]),ylab=bquote(X[2]),main="Pointwise bias")

x1.plot <- x2.plot <- seq(0,1,length.out=1000)
X.plot <- expand.grid(x1.plot,x2.plot)
g1 <- gamma_fun(X.plot,1)
g0 <- gamma_fun(X.plot,0)
bias <- matrix(g1-g0, nrow=1000, ncol=1000)
colMap <- colorRampPalette(c("red","white","blue" ))(100)

image(x1.plot,x2.plot,z=bias,main="Pointwise bias",col=colMap)
legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
       c(round(min(bias),2),0,round(max(bias),2)), fill = colMap[c(1, 50, 100)], xpd = NA)
## Plot per PTA region
plot(0.5,0.5,type="n",xlim=c(0,1),ylim=c(0,1),xlab=bquote(X[1]),ylab=bquote(X[2]),main="Bias per id. region")
for (j in 1:length(bias.per.region)){
  id1 <- preds[X$g==1]==bias.per.region[j]
  id0 <- preds[X$g==0]==bias.per.region[j]
  points(x1[id1,1],x1[id1,2],pch=15,col=rgb(r[j],g[j],b[j],alpha = 1))
  points(x0[id0,1],x0[id0,2],pch=15,col=rgb(r[j],g[j],b[j],alpha = 1))
}
## Plot X densities
X.plot <- expand.grid(seq(0,1,length.out=100),seq(0,1,length.out=100))
### G=1
P1 <- apply(X.plot,1,function(i) dbeta(i[1],x1.shape1.g1,x1.shape2.g1)*dbeta(i[2],x2.shape1.g1,x2.shape2.g1))
P1 <- matrix(P1,nrow=length(unique(X.plot[,1])),ncol=length(unique(X.plot[,2])),byrow=FALSE)

contour(unique(X.plot[,1]),unique(X.plot[,2]),P1,xlim=c(0,1),ylim=c(0,1),xlab=bquote(X[1]),ylab=bquote(X[2]),main=bquote(P(X~"|"~g==1)))
### G=0
P0 <- apply(X.plot,1,function(i) dbeta(i[1],x1.shape1.g0,x1.shape2.g0)*dbeta(i[2],x2.shape1.g0,x2.shape2.g0))
P0 <- matrix(P0,nrow=length(unique(X.plot[,1])),ncol=length(unique(X.plot[,2])),byrow=FALSE)

contour(unique(X.plot[,1]),unique(X.plot[,2]),P0,xlim=c(0,1),ylim=c(0,1),xlab=bquote(X[1]),ylab=bquote(X[2]),main=bquote(P(X~"|"~g==0)))
##
par(mfrow=c(1,1))
