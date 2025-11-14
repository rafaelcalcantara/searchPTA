## Setup-----------------------------------------------------------------------
seed <- 007
set.seed(seed)
n <- 5000
## Generate data---------------------------------------------------------------
### Fixed features
t <- c(rep(-1,n),rep(0,n),rep(1,n))
n1 <- n %/% 2
n0 <- n-n1
g <- c(rep(1,n1),rep(0,n0))
z <- c(g,g,g)*(t==1)
### X1: capital stock
x1 <- g*rgamma(n,3,3) + (1-g)*rgamma(n,6,3)
# x1 <- rgamma(n,8,1)
### X2: number of employees
x2 <- ifelse(g==1,sample(c("small","large"),n,replace=TRUE,prob=c(0.6,0.4)),sample(c("small","large"),n,replace=TRUE,prob=c(0.4,0.6)))
# x2 <- sample(1:5,n,replace=TRUE,prob=c(0.1,0.2,0.3,0.3,0.1))
### X3: sector
x3 <- ifelse(g==1,sample(c("a","b","c"),n,replace=TRUE,prob=c(0.4,0.25,0.35)),sample(c("a","b","c"),n,replace=TRUE,prob=rep(1/3,3)))
# x3 <- sample(c("a","c"),n,replace=TRUE,prob=rep(1/2,2))
#### Visualize feature distribution per g
par(mfrow=c(3,2),bty="n",cex.axis=0.7)
hist(x1[g==1])
hist(x1[g==0])
barplot(table(x2[g==1]))
barplot(table(x2[g==0]))
barplot(table(x3[g==1]))
barplot(table(x3[g==0]))
### Construct DGP functions
mu_fun <- function(x1,x2,g)
{
  ## Intercept
  a0 <- ifelse(x2=="large",0.5,0.2)
  a1 <- ifelse(g==1,0,0.2)
  a <- a0+a1
  ## x1 coefficient
  b <- ifelse(g==1,0.3,0.6)
  ## Return mu
  return(a+b*x1+sin(0.75*pi*x1))
}
beta_fun <- function(x3,g)
{
  ## Intercept
  a0 <- ifelse(x2=="large",0.4,0.1)
  a1 <- ifelse(x3=="a",0.6, ifelse(x3=="b",0.45,0.15))
  a2 <- ifelse(g==1,0.1,0.3)
  a <- a0+a1+a2
  ## x1 coefficient
  b <- ifelse(g==1,0.4,0.1)
  ## Return mu
  return(a+b*x1-exp(-5*x1))
}
gamma_fun <- function(x3,g) 2*beta_fun(x3,g)
tau_fun <- function(x2,x3)
{
  a0 <- ifelse(x2=="large",0.1,0.8)
  a1 <- ifelse(x3=="a",0.3, 0.1)
  return(a0+a1)
}
alpha_fun <- function(x2)
{
  a0 <- ifelse(x2=="large",0,0.4)
  a1 <- ifelse(x3=="c",0.5,0)
  return(a0+a1)
}
### Calculate function values
mu <- mu_fun(x1,x2,g)
beta <- beta_fun(x3,g)
gamma <- gamma_fun(x3,g)
tau <- tau_fun(x2,x3)
alpha <- alpha_fun(x2)
### Define PTA regions
S <- rep(3,n)
for (i in 1:n)
{
  if (x2[i] == "large") S[i] <- 1
  if (x2[i] == "small" & x3[i] %in% c("a","b")) S[i] <- 2
}
### Establish PTA in regions 1 and 2
beta_avg_1_1 <- mean(beta[S==1 & g==1])
beta_avg_1_0 <- mean(beta[S==1 & g==0])
beta_avg_2_1 <- mean(beta[S==2 & g==1])
beta_avg_2_0 <- mean(beta[S==2 & g==0])
beta_avg_3_1 <- mean(beta[S==3 & g==1])
beta_avg_3_0 <- mean(beta[S==3 & g==0])
beta <- beta - (S==1 & g==1)*(beta_avg_1_1-beta_avg_1_0) - (S==2 & g==1)*(beta_avg_2_1-beta_avg_2_0) #- (S==3 & g==1)*(beta_avg_3_1-beta_avg_3_0 - 0.2)
gamma_avg_1_1 <- mean(gamma[S==1 & g==1])
gamma_avg_1_0 <- mean(gamma[S==1 & g==0])
gamma_avg_2_1 <- mean(gamma[S==2 & g==1])
gamma_avg_2_0 <- mean(gamma[S==2 & g==0])
gamma_avg_3_1 <- mean(gamma[S==3 & g==1])
gamma_avg_3_0 <- mean(gamma[S==3 & g==0])
gamma <- gamma - (S==1 & g==1)*(gamma_avg_1_1-gamma_avg_1_0) - (S==2 & g==1)*(gamma_avg_2_1-gamma_avg_2_0) #- (S==3 & g==1)*(gamma_avg_3_1-gamma_avg_3_0 - 0.2)
### Generate outcome
Ey <- rep(mu,3) - rep(gamma,3)*(t==-1) + rep(beta,3)*(t==1) +  rep(tau,3)*z + rep(alpha,3)*z*(t==1)
error.sd <- 0.2*sd(Ey)
y <- Ey + rnorm(n,0,error.sd)
### Check that PTA holds in stipulated regions
mean(beta[S==1 & g==1])-mean(beta[S==1 & g==0])
mean(beta[S==2 & g==1])-mean(beta[S==2 & g==0])
mean(beta[S==3 & g==1])-mean(beta[S==3 & g==0])
mean(gamma[S==1 & g==1])-mean(gamma[S==1 & g==0])
mean(gamma[S==2 & g==1])-mean(gamma[S==2 & g==0])
mean(gamma[S==3 & g==1])-mean(gamma[S==3 & g==0])
mean(Ey[S==1 & g==1 & t==0])-mean(Ey[S==1 & g==0 & t==0]) - (mean(Ey[S==1 & g==1 & t==-1])-mean(Ey[S==1 & g==0 & t==-1]))
mean(Ey[S==2 & g==1 & t==0])-mean(Ey[S==2 & g==0 & t==0]) - (mean(Ey[S==2 & g==1 & t==-1])-mean(Ey[S==2 & g==0 & t==-1]))
mean(Ey[S==3 & g==1 & t==0])-mean(Ey[S==3 & g==0 & t==0]) - (mean(Ey[S==3 & g==1 & t==-1])-mean(Ey[S==3 & g==0 & t==-1]))
mean(y[S==1 & g==1 & t==0])-mean(y[S==1 & g==0 & t==0]) - (mean(y[S==1 & g==1 & t==-1])-mean(y[S==1 & g==0 & t==-1]))
mean(y[S==2 & g==1 & t==0])-mean(y[S==2 & g==0 & t==0]) - (mean(y[S==2 & g==1 & t==-1])-mean(y[S==2 & g==0 & t==-1]))
mean(y[S==3 & g==1 & t==0])-mean(y[S==3 & g==0 & t==0]) - (mean(y[S==3 & g==1 & t==-1])-mean(y[S==3 & g==0 & t==-1]))
## Estimation------------------------------------------------------------------
par(mfrow=c(1,2))

eps <- 0.01

df <- data.frame(X1=x1,X2=factor(x2,ordered=TRUE),X3=factor(x3,ordered=FALSE),g=g)
test <- searchPTA::searchPTA(df,gamma,gamma,beta+tau+alpha,beta,epsilon=eps,minsplit=1,minbucket=1,cp=0.001,maxdepth=30)
rpart.plot::rpart.plot(test$cart)
table(test$beta.diff,S)

df <- data.frame(X=S,g=g)
test <- searchPTA::searchPTA(df,gamma,gamma,beta+tau+alpha,beta,epsilon=eps,minsplit=1,minbucket=1,cp=0,maxdepth=30)
rpart.plot::rpart.plot(test$cart)
table(test$beta.diff,S)
