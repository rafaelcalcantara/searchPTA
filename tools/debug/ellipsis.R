###############################################################################
## Benchmarking functions time and memory-wise
###############################################################################
devtools::load_all()
epsilon <- 0.01
set.seed(0)
# We check the functions with continuous data (no time issues for discrete data)
## Generate data
### Define parameters for each group
mu_1 <- c(1,3,1,3)
mu_0 <- c(1,2,3,4)
mu <- cbind(mu_0,mu_1)
###
beta_1 <- c(1,1,1,1)
beta_0 <- c(7,1,2,0)
beta <- cbind(beta_0,beta_1)
###
# alpha_1 <- seq(0,1,length.out = 4)
# alpha_0 <- seq(1,0,length.out = 4)
# alpha <- cbind(alpha_0, alpha_1)
# ###
# tau_1 <- c(1,2,3,4)
# tau_0 <- c(2,1,4,3)
# tau <- cbind(tau_0, tau_1)
### Generate features
n <- 900
t <- c(rep(-1,n/3),rep(0,n/3),rep(1,n/3))
g <- rep(c(0,1),n/2)
z <- g*(t==1)
#########
x1 <- runif(n/3)
x2 <- runif(n/3)
x <- rep(1,n/3)
for (i in 1:(n/3))
{
  if (x1[i]>0.5 & x2[i]<0.5) x[i] <- 2
  if (x1[i]>0.5 & x2[i]>0.5) x[i] <- 3
  if (x1[i]<0.5 & x2[i]>0.5) x[i] <- 4
}
x <- c(x,x,x)
#########
alpha_1 <- 0.5*pnorm(x1,0.5,0.2)
alpha_0 <- pnorm(x1,0.5,1)
alpha <- c(alpha_1,alpha_1,alpha_1)*g + c(alpha_0,alpha_0,alpha_0)*(1-g)
tau_1 <- 0.5*pnorm(-x2,-0.5,0.1)
tau_0 <- pnorm(-x2,-0.5,1)
tau <- c(tau_1,tau_1,tau_1)*g + c(tau_0,tau_0,tau_0)*(1-g)
### Generate outcome
idx <- cbind(x,g+1)
Ey <- mu[idx] + beta[idx]*t + tau*z + alpha*z*t
y <- Ey + 0.5*sd(Ey)*rnorm(n)
df <- data.frame(y,t,g,x1=c(x1,x1,x1),x2=c(x2,x2,x2))
df_main <- subset(df,t>=0)
df_placebo <- subset(df,t<=0)
### Checking the functions for true beta_1 and beta_0
b1 <- beta_1[x]
b0 <- beta_0[x]
tau1 <- c(tau_1,tau_1,tau_1)
alpha1 <- c(alpha_1,alpha_1,alpha_1)
######
test <- searchPTA(x=subset(df_main,select=4:5),beta1_placebo=b1[t<=0],beta0_placebo=b0[t<=0],beta1_main=b1[t>=0]+alpha1[t>=0]+tau1[t>=0],beta0_main=b0[t>=0],epsilon=epsilon,saveCART = FALSE,cp=0)
######
col <- c("black","blue","red")
bg <- scales::alpha(col,c(0,0.4,0.4))
layout(matrix(c(1,3,2,3),ncol=2),heights = c(3,1))
par(mar=c(4,3,3,1))
matplot(df_main$x1, cbind((tau1+alpha1)[t>=0],test$catt),bty="n",xlab=bquote(X[1]),ylab="CATT",pch=21,col=col[1:2],bg=bg[1:2],cex=0.7)
points(cbind(df_main$x1[x[t>=0]==2],0.02+mean((tau1+alpha1)[x==2 & t>=0])),pch=21,col=col[3],bg=bg[3],cex=0.7)
points(cbind(df_main$x1[x[t>=0] %in% 3:4],0.025+mean((tau1+alpha1)[x %in% 3:4 & t>=0])),pch=21,col=col[3],bg=bg[3],cex=0.7)
# text(x=mean(df_main$x1[x[t>=0]==1]),y=unique((tau1+alpha1)+2*pnorm(x1+x2,x1+x2,x1+x2))[x==1],labels=1,pos=3)
# text(x=mean(df_main$x1[x[t>=0]==2]),y=unique((tau1+alpha1))[x==2],labels=2,pos=3)
# text(x=mean(df_main$x1[x[t>=0]==3]),y=unique((tau1+alpha1))[x==3],labels=3,pos=3)
# text(x=mean(df_main$x1[x[t>=0]==4]),y=unique((tau1+alpha1))[x==4],labels=4,pos=1)
matplot(df_main$x2[t>=0], cbind((tau1+alpha1)[t>=0],test$catt),bty="n",xlab=bquote(X[2]),ylab="CATT",pch=21,col=col[1:2],bg=bg[1:2],cex=0.7)
points(cbind(df_main$x2[x[t>=0]==2],0.02+mean((tau1+alpha1)[x==2 & t>=0])),pch=21,col=col[3],bg=bg[3],cex=0.7)
points(cbind(df_main$x2[x[t>=0] %in% 3:4],0.035+mean((tau1+alpha1)[x %in% 3:4 & t>=0])),pch=21,col=col[3],bg=bg[3],cex=0.7)
# text(x=mean(df_main$x2[x[t>=0]==1]),y=unique((tau1+alpha1))[x==1],labels=1,pos=3)
# text(x=mean(df_main$x2[x[t>=0]==2]),y=unique((tau1+alpha1))[x==2],labels=2,pos=3)
# text(x=mean(df_main$x2[x[t>=0]==3]),y=unique((tau1+alpha1))[x==3],labels=3,pos=3)
# text(x=mean(df_main$x2[x[t>=0]==4]),y=unique((tau1+alpha1))[x==4],labels=4,pos=1)
plot.new()
legend("center",legend=c(bquote(CATT(x[i])),bquote(CATT(x %in% S)),"CART"),col=col[c(1,3,2)],pt.bg=bg[c(1,3,2)],pch=21,ncol=3,cex=0.9)
######
dev.off()

library(scatterplot3d)
par(mfrow=c(1,2))
####
temp <- expand.grid(x1,x2)
alpha1 <- 0.5*pnorm(temp[,1],0.5,0.2)
tau1 <- 0.5*pnorm(-temp[,1],-0.5,0.1)
x3d <- temp[,1]
y3d <- temp[,2]
z3d <- tau1+alpha1
zlim <- c(min(z3d),max(z3d))
scatterplot3d(x3d,y3d,z3d,zlim=zlim,pch=20,box = FALSE,highlight.3d = TRUE,type='p')
####
temp <- expand.grid(seq(0,1,0.01),seq(0,1,0.01))
alpha1 <- 0.5*pnorm(temp[,1],0.5,0.2)
tau1 <- 0.5*pnorm(-temp[,1],-0.5,0.1)
x3d <- temp[,1]
y3d <- temp[,2]
z3d.temp <- tau1+alpha1
z3d <- rep(NA,length(z3d.temp))
for (i in 1:length(z3d))
{
  if (x3d[i] >= 0.5 & y3d[i] <= 0.5) z3d[i] <- mean(z3d.temp[x3d >= 0.5 & y3d <= 0.5])
  if (y3d[i] >= 0.5) z3d[i] <- mean(z3d.temp[y3d >= 0.5])
}
scatterplot3d(x3d,y3d,z3d,pch=20,box = FALSE,highlight.3d = TRUE,type='l')
