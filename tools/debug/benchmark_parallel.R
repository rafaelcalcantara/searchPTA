###############################################################################
## Benchmarking functions for multiple runs (e.g. multiple BART draws)
###############################################################################
devtools::load_all()
epsilon <- 0.01
m <- 1000 ## Times to run function
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
alpha_1 <- seq(0,1,length.out = 4)
alpha_0 <- seq(1,0,length.out = 4)
alpha <- cbind(alpha_0, alpha_1)
###
tau_1 <- c(1,2,3,4)
tau_0 <- c(2,1,4,3)
tau <- cbind(tau_0, tau_1)
### Generate features
n <- 1500
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
## For loop
t0 <- Sys.time()
for (i in 1:m)
{
  searchPTA(subset(df_main,select=4),b1[t<=0],b0[t<=0],b1[t>=0]+alpha1[t>=0]+tau1[t>=0],b0[t>=0],saveCART = FALSE,epsilon)
}
t1 <- Sys.time()
time.for <- t1-t0
## Lapply
t0 <- Sys.time()
lapply(1:m, function(i) searchPTA(subset(df_main,select=4),b1[t<=0],b0[t<=0],b1[t>=0]+alpha1[t>=0]+tau1[t>=0],b0[t>=0],saveCART = FALSE,epsilon))
t1 <- Sys.time()
time.lapply <- t1-t0
## Check times
print("Time for loop")
print(time.for)
print("Time lapply")
print(time.lapply)
