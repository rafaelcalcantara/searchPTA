###############################################################################
## Benchmarking functions time and memory-wise
###############################################################################
devtools::load_all()
library(searchPTA)
epsilon <- 0.01
m <- 1000 ## Times to run function
set.seed(0)
cp <- seq(0,0.0025,length=6) # grid of cp parameters to use for the rpart function
# 1) First, we check the functions with discrete data--------------------------
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
alpha_1 <- seq(0,1,length.out = 4)
alpha_0 <- seq(1,0,length.out = 4)
alpha <- cbind(alpha_0, alpha_1)
###
tau_1 <- c(1,2,3,4)
tau_0 <- c(2,1,4,3)
tau <- cbind(tau_0, tau_1)
### Generate features
n <- 4500
t <- c(rep(-1,n/3),rep(0,n/3),rep(1,n/3))
g <- rep(c(0,1),n/2)
z <- g*(t==1)
x <- sample(1:4,n/3,replace=TRUE)
x <- c(x,x,x)
### Generate outcome
idx <- cbind(x,g+1)
Ey <- mu[idx] + beta[idx]*t + tau[idx]*z + alpha[idx]*z*t
y <- Ey + 0.1*sd(Ey)*rnorm(n)
df <- data.frame(y,t,g,x=as.factor(x))
df_main <- subset(df,t>=0)
df_placebo <- subset(df,t<=0)
### Checking the functions for true beta_1 and beta_0
b1 <- beta_1[x]
b0 <- beta_0[x]
tau1 <- tau_1[x]
alpha1 <- alpha_1[x]
time <- matrix(0,nrow=m,ncol=length(cp))
for (i in 1:length(cp))
{
  for (j in 1:m)
  {
    t0 <- Sys.time()
    test <- searchPTA(subset(df_main,select=4),b1[t<=0],b0[t<=0],b1[t>=0]+alpha1[t>=0]+tau1[t>=0],b0[t>=0],saveCART = FALSE,epsilon,cp=cp[i])
    t1 <- Sys.time()
    time[j,i] <- as.numeric(t1-t0)
  }
}
## Check times
hist(time)
## Check if results are right
summary(test$catt[x[t>=0] %in% 2:4])
mean(tau_1[2:4] + alpha_1[2:4])
summary(test$catt[x[t>=0]==2])
tau_1[2] + alpha_1[2]
summary(test$catt[x[t>=0] %in% 3:4])
mean(tau_1[3:4] + alpha_1[3:4])
# 2) Now, we check the functions with continuous data--------------------------
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
alpha_1 <- seq(0,1,length.out = 4)
alpha_0 <- seq(1,0,length.out = 4)
alpha <- cbind(alpha_0, alpha_1)
###
tau_1 <- c(1,2,3,4)
tau_0 <- c(2,1,4,3)
tau <- cbind(tau_0, tau_1)
### Generate features
n <- 4500
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
### Generate outcome
idx <- cbind(x,g+1)
Ey <- mu[idx] + beta[idx]*t + tau[idx]*z + alpha[idx]*z*t
y <- Ey + 0.5*sd(Ey)*rnorm(n)
df <- data.frame(y,t,g,x1=c(x1,x1,x1),x2=c(x2,x2,x2))
df_main <- subset(df,t>=0)
df_placebo <- subset(df,t<=0)
### Checking the functions for true beta_1 and beta_0
b1 <- beta_1[x]
b0 <- beta_0[x]
tau1 <- tau_1[x]
alpha1 <- alpha_1[x]
time <- matrix(0,nrow=m,ncol=length(cp))
for (i in 1:length(cp))
{
  for (j in 1:m)
  {
    t0 <- Sys.time()
    test <- searchPTA(subset(df_main,select=4),b1[t<=0],b0[t<=0],b1[t>=0]+alpha1[t>=0]+tau1[t>=0],b0[t>=0],saveCART = FALSE,epsilon,cp=cp[i])
    t1 <- Sys.time()
    time[j,i] <- as.numeric(t1-t0)
  }
}
######
profvis::profvis({
  test <- searchPTA(subset(df_main,select=4),b1[t<=0],b0[t<=0],b1[t>=0]+alpha1[t>=0]+tau1[t>=0],b0[t>=0],saveCART = FALSE,epsilon,cp=cp[i])
})
## Check times
summary(time)
hist(time)
## Check if results are right
summary(test$catt[x[t>=0] %in% 2:4])
mean(tau_1[2:4] + alpha_1[2:4])
summary(test$catt[x[t>=0]==2])
tau_1[2] + alpha_1[2]
summary(test$catt[x[t>=0] %in% 3:4])
mean(tau_1[3:4] + alpha_1[3:4])
