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
### Continuous functions
mu.cont.fun <- function(x1,g)
{
  k <- ifelse(g==1,1,0.75)
  return(1/(10+exp(-x1))+0.1*sin(k*x1))
}
gamma.cont.fun <- function(x1,g)
{
  k <- ifelse(g==1,2,1.5)
  return((0.3*x1)^k)
}
tau.cont.fun <- function(x1,g)
{
  k <- ifelse(g==1,3,6)
  return(0.5*pgamma(x1,k,3))
}
alpha.cont.fun <- function(x1,g)
{
  k <- ifelse(g==1,0.7,0.3)
  return(exp(-0.2*x1))
}
par(mfrow=c(1,1))
for (i in c(0,1))
{
  plot(x1,mu.cont.fun(x1,i))
  plot(x1,gamma.cont.fun(x1,i))
  plot(x1,tau.cont.fun(x1,i))
  plot(x1,alpha.cont.fun(x1,i))
  plot(x1,tau.cont.fun(x1,i)+alpha.cont.fun(x1,i))
}
### Discrete functions
mu.disc.fun <- function(x2,x3,g)
{
  out <- rep(NA,length(g))
  for (i in 1:length(out))
  {
    if (x3[i] == "a")
    {
      out[i] <- 0.6
    } else if (x3[i] == "b")
    {
      out[i] <- 1
    } else if (x3[i] == "c")
    {
      out[i] <- 0.8
    }
  }
  out <- out + (x2=="large")*0.5 + (g==1)*0.2
  return(out)
}
gamma.disc.fun <- function(x2,x3,g)
{
  out <- rep(NA,length(g))
  for (i in 1:length(out))
  {
    if (x3[i] == "a")
    {
      out[i] <- 0.3
    } else if (x3[i] == "b")
    {
      out[i] <- 0.5
    } else if (x3[i] == "c")
    {
      out[i] <- 0.4
    }
  }
  out <- out + (x2=="large")*0.2 + (g==1)*0.4
  return(out)
}
tau.disc.fun <- function(x2,x3,g)
{
  out <- rep(NA,length(g))
  for (i in 1:length(out))
  {
    if (x3[i] == "a")
    {
      out[i] <- 0.1
    } else if (x3[i] == "b")
    {
      out[i] <- 0.15
    } else if (x3[i] == "c")
    {
      out[i] <- 0.3
    }
  }
  out <- out + (x2=="large")*0.4 + (g==1)*0.3
  return(out)
}
alpha.disc.fun <- function(x2,x3,g)
{
  out <- rep(NA,length(g))
  for (i in 1:length(out))
  {
    if (x3[i] == "a")
    {
      out[i] <- 0.2
    } else if (x3[i] == "b")
    {
      out[i] <- 0.1
    } else if (x3[i] == "c")
    {
      out[i] <- 0.45
    }
  }
  out <- out + (x2=="large")*0.3 + (g==1)*0.1
  return(out)
}
### Construct DGP functions
#### g=1
mu_1 <- mu.disc.fun(x2,x3,g=1) + mu.cont.fun(x1,g=1)
gamma_1 <- gamma.disc.fun(x2,x3,g=1) + gamma.cont.fun(x1,g=1)
beta_1 <- gamma_1
tau_1 <- tau.disc.fun(x2,x3,g=1) + tau.cont.fun(x1,g=1)
alpha_1 <- alpha.disc.fun(x2,x3,g=1) + alpha.cont.fun(x1,g=1)
#### g=0
mu_0 <- mu.disc.fun(x2,x3,g=0) + mu.cont.fun(x1,g=0)
gamma_0 <- gamma.disc.fun(x2,x3,g=0) + gamma.cont.fun(x1,g=0)
beta_0 <- gamma_0
tau_0 <- tau.disc.fun(x2,x3,g=0) + tau.cont.fun(x1,g=0)
alpha_0 <- alpha.disc.fun(x2,x3,g=0) + alpha.cont.fun(x1,g=0)
### Define PTA regions
S <- rep(3,n)
for (i in 1:n)
{
  if (x3[i] == "a" & x1[i] > median(x1)) S[i] <- 1
  if (x3[i] == "c" & x2[i] == "large") S[i] <- 2
}
### Establish PTA in regions 1 and 2
beta_avg_1_1 <- mean(beta_1[S==1 & g==1])
beta_avg_1_0 <- mean(beta_0[S==1 & g==0])
beta_avg_2_1 <- mean(beta_1[S==2 & g==1])
beta_avg_2_0 <- mean(beta_0[S==2 & g==0])
beta_1 <- beta_1 - (S==1 & g==1)*(beta_avg_1_1-beta_avg_1_0) - (S==2 & g==1)*(beta_avg_2_1-beta_avg_2_0)
gamma_avg_1_1 <- mean(gamma_1[S==1 & g==1])
gamma_avg_1_0 <- mean(gamma_0[S==1 & g==0])
gamma_avg_2_1 <- mean(gamma_1[S==2 & g==1])
gamma_avg_2_0 <- mean(gamma_0[S==2 & g==0])
gamma_1 <- gamma_1 - (S==1 & g==1)*(gamma_avg_1_1-gamma_avg_1_0) - (S==2 & g==1)*(gamma_avg_2_1-gamma_avg_2_0)
####
mu <- mu_1*g + mu_0*(1-g)
gamma <- gamma_1*g + gamma_0*(1-g)
beta <- beta_1*g + beta_0*(1-g)
alpha <- alpha_1*g + alpha_0*(1-g)
tau <- tau_1*g + tau_0*(1-g)
### Generate outcome
Ey <- rep(mu,3) - rep(gamma,3)*(t==-1) + rep(beta,3)*(t==1) +  rep(tau,3)*z + rep(alpha,3)*z*(t==1)
error.sd <- 0.5*sd(Ey)
y <- Ey + rnorm(n,0,error.sd)
### Compare y and Ey
par(mfrow=c(1,1),bty="n",cex.axis=0.7)
plot(y,Ey,cex=0.7)
abline(a=0,b=1,lty=2)
### Check that PTA holds in stipulated regions
mean(beta[S==1 & g==1])-mean(beta[S==1 & g==0])
mean(beta[S==2 & g==1])-mean(beta[S==2 & g==0])
mean(gamma[S==1 & g==1])-mean(gamma[S==1 & g==0])
mean(gamma[S==2 & g==1])-mean(gamma[S==2 & g==0])
mean(Ey[S==1 & g==1 & t==0])-mean(Ey[S==1 & g==0 & t==0]) - (mean(Ey[S==1 & g==1 & t==-1])-mean(Ey[S==1 & g==0 & t==-1]))
mean(Ey[S==2 & g==1 & t==0])-mean(Ey[S==2 & g==0 & t==0]) - (mean(Ey[S==2 & g==1 & t==-1])-mean(Ey[S==2 & g==0 & t==-1]))
mean(y[S==1 & g==1 & t==0])-mean(y[S==1 & g==0 & t==0]) - (mean(y[S==1 & g==1 & t==-1])-mean(y[S==1 & g==0 & t==-1]))
mean(y[S==2 & g==1 & t==0])-mean(y[S==2 & g==0 & t==0]) - (mean(y[S==2 & g==1 & t==-1])-mean(y[S==2 & g==0 & t==-1]))
## Estimation------------------------------------------------------------------
### Apply CART to true data to check that algorithm does what is expected
bta1 <- beta_1 + tau_1 + alpha_1
epsilon <- 0.05
df <- data.frame(X1=x1,X2=factor(x2,ordered=TRUE),X3=factor(x3,ordered=FALSE),g=g)
catt_region <- searchPTA::searchPTA(x=df,delta1_aux=gamma_1,delta0_aux=gamma_0,delta1_main=bta1,delta0_main=beta_0,epsilon=epsilon,minsplit=0,minbucket=0,cp=0)
#### Check regions
lapply(1:ncol(catt_region$regions), function(i) table(S,catt_region$regions[,i]))
lapply(1:ncol(catt_region$regions), function(i) mean(gamma_1[g==1 & catt_region$regions[,i]])-mean(gamma_0[g==0 & catt_region$regions[,i]]))
lapply(1:ncol(catt_region$regions), function(i) mean(beta_1[g==1 & catt_region$regions[,i]])-mean(beta_0[g==0 & catt_region$regions[,i]]))
#####
estimated.regions <- sapply(1:nrow(catt_region$regions), function(i) ifelse(length(which(catt_region$regions[i,]))==0,"unidentified",letters[which(catt_region$regions[i,])]))
table(S,estimated.regions)
### Fit model to data and try to recover regions from model fit
library(stochtree)
library(foreach)
library(doParallel)
num_chains <- 20
draws <- 1000
burnin <- 50
ncores <- 10
num_trees <- 100
max_depth <- 20
model.df <- data.frame(y,t,g=c(g,g,g),x1=c(x1,x1,x1),x2=c(x2,x2,x2),x3=c(x3,x3,x3))
model.df$x2 <- factor(model.df$x2,ordered=TRUE)
model.df$x3 <- factor(model.df$x3,ordered=FALSE)
out <- list(gamma1=NULL,gamma0=NULL,beta1=NULL,beta0=NULL)
df.x <- subset(model.df,select=4:ncol(model.df))
keep.draws <- (burnin+1):(draws+burnin)
## Generate time dummies
model.df$pre <- -as.numeric(model.df$t==-1)
model.df$post <- as.numeric(model.df$t==1)
## Lists of parameters for the Stochtree BART function
bart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
bart.mean.parmlist <- list(num_trees=num_trees, min_samples_leaf=20, alpha=0.95, beta=2,
                           max_depth=max_depth, sample_sigma2_leaf=FALSE,
                           sigma2_leaf_init = diag(rep(0.1/num_trees,6)))
## Set basis vector for leaf regressions
Psi <- cbind(1-model.df$g,model.df$g,(1-model.df$g)*model.df$pre,model.df$g*model.df$pre,(1-model.df$g)*model.df$post,model.df$g*model.df$post)
## Model fit
bart.fit = stochtree::bart(X_train=df.x, y_train=model.df$y,
                           leaf_basis_train = Psi, mean_forest_params=bart.mean.parmlist,
                           general_params=bart.global.parmlist,
                           num_mcmc=draws+burnin,num_gfr=0)
## Extract parameters for g=1 and g=0
### We use the 'predict_raw' function from the bart.fit object
### 1) Pre-process test data to the format expected by the predict_raw function
X_test <- stochtree::preprocessPredictionData(df.x,bart.fit$train_set_metadata)
### 2) Create 'ForestDataset' structure to be read by the predict_raw function
X_test <- stochtree::createForestDataset(X_test,basis=Psi)
### 3) Use predict_raw to extract betas for each group. Because the basis is (1,g,t,g*t), we need the 3rd and 4th coefficients
coefs <- bart.fit$mean_forests$predict_raw(X_test)
### 4) Save results; the predict_raw output is a (n,p,m) array, where n=num_obs, p is the length of the basis vector, m is the number of MCMC samples
#### Extracting coefs[,3,] gives us a n*m matrix of individual draws for the coefficient on (1-g)*D_{-1}, and similarly coefs[,4,] for g*D_{-1}
#### Extracting coefs[,5,] gives us a n*m matrix of individual draws for the coefficient on (1-g)*D_{1}, and similarly coefs[,6,] for g*D_{1}
#### The raw predictions need to be scaled back to match the actual data (i.e. the raw predictions are based on (Y-Y_bar)/sd(Y))
#### Note that, because those two coefficients represent differences, the Y_bar term cancels out and we only need to multiply by sd(Y):
#### E.g. the coefficient on t is equal to E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=1]-E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=0] = (E[Y|X=x,G=0,T=1]-E[Y|X=x,G=0,T=0])/sd(Y)
out$gamma1 <- coefs[,4,keep.draws]*sd(model.df$y)
out$gamma0 <- coefs[,3,keep.draws]*sd(model.df$y)
out$beta1_tau1_alpha1 <- coefs[,6,keep.draws]*sd(model.df$y)
out$beta0 <- coefs[,5,keep.draws]*sd(model.df$y)
# cl <- makeCluster(ncores)
# registerDoParallel(cl)
# bart_model_outputs <- foreach (i = 1:num_chains) %dopar% {
#   out <- list(gamma1=NULL,gamma0=NULL,beta1=NULL,beta0=NULL)
#   df.x <- subset(model.df,select=4:ncol(model.df))
#   keep.draws <- (burnin+1):(draws+burnin)
#   ## Generate time dummies
#   model.df$pre <- -as.numeric(model.df$t==-1)
#   model.df$post <- as.numeric(model.df$t==1)
#   ## Lists of parameters for the Stochtree BART function
#   bart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
#   bart.mean.parmlist <- list(num_trees=num_trees, min_samples_leaf=20, alpha=0.95, beta=2,
#                              max_depth=max_depth, sample_sigma2_leaf=FALSE,
#                              sigma2_leaf_init = diag(rep(0.1/num_trees,6)))
#   ## Set basis vector for leaf regressions
#   Psi <- cbind(1-model.df$g,model.df$g,(1-model.df$g)*model.df$pre,model.df$g*model.df$pre,(1-model.df$g)*model.df$post,model.df$g*model.df$post)
#   ## Model fit
#   bart.fit = stochtree::bart(X_train=df.x, y_train=model.df$y,
#                              leaf_basis_train = Psi, mean_forest_params=bart.mean.parmlist,
#                              general_params=bart.global.parmlist,
#                              num_mcmc=draws+burnin,num_gfr=0)
#   ## Extract parameters for g=1 and g=0
#   ### We use the 'predict_raw' function from the bart.fit object
#   ### 1) Pre-process test data to the format expected by the predict_raw function
#   X_test <- stochtree::preprocessPredictionData(df.x,bart.fit$train_set_metadata)
#   ### 2) Create 'ForestDataset' structure to be read by the predict_raw function
#   X_test <- stochtree::createForestDataset(X_test,basis=Psi)
#   ### 3) Use predict_raw to extract betas for each group. Because the basis is (1,g,t,g*t), we need the 3rd and 4th coefficients
#   coefs <- bart.fit$mean_forests$predict_raw(X_test)
#   ### 4) Save results; the predict_raw output is a (n,p,m) array, where n=num_obs, p is the length of the basis vector, m is the number of MCMC samples
#   #### Extracting coefs[,3,] gives us a n*m matrix of individual draws for the coefficient on (1-g)*D_{-1}, and similarly coefs[,4,] for g*D_{-1}
#   #### Extracting coefs[,5,] gives us a n*m matrix of individual draws for the coefficient on (1-g)*D_{1}, and similarly coefs[,6,] for g*D_{1}
#   #### The raw predictions need to be scaled back to match the actual data (i.e. the raw predictions are based on (Y-Y_bar)/sd(Y))
#   #### Note that, because those two coefficients represent differences, the Y_bar term cancels out and we only need to multiply by sd(Y):
#   #### E.g. the coefficient on t is equal to E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=1]-E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=0] = (E[Y|X=x,G=0,T=1]-E[Y|X=x,G=0,T=0])/sd(Y)
#   out$gamma1 <- coefs[,4,keep.draws]*sd(model.df$y)
#   out$gamma0 <- coefs[,3,keep.draws]*sd(model.df$y)
#   out$beta1_tau1_alpha1 <- coefs[,6,keep.draws]*sd(model.df$y)
#   out$beta0 <- coefs[,5,keep.draws]*sd(model.df$y)
#   out
# }
# stopCluster(cl)
#### Obtain CATT draws per i
# bart_catt_i <- do.call("cbind",lapply(bart_model_outputs, function(i) (i$beta1_tau1_alpha1-i$beta0)[1:n,]))
# bta1 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$beta1_tau1_alpha1[1:n,]))
# b0 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$beta0[1:n,]))
# g1 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$gamma1[1:n,]))
# g0 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$gamma0[1:n,]))
# #### We check the quality of the BART estimates by verifying the CATT posterior for the true PTA regions
# par(mfrow=c(2,2),bty="n")
# plot((beta_1+tau_1+alpha_1)[g==1 & S==1],Ey[g==1 & t==1 & S==1]-Ey[g==1 & t==0 & S==1])
# abline(a=0,b=1)
# plot((beta_1+tau_1+alpha_1)[g==1 & S==2],Ey[g==1 & t==1 & S==2]-Ey[g==1 & t==0 & S==2])
# abline(a=0,b=1)
# plot((beta_0)[g==0 & S==1],Ey[g==0 & t==1 & S==1]-Ey[g==0 & t==0 & S==1])
# abline(a=0,b=1)
# plot((beta_0)[g==0 & S==2],Ey[g==0 & t==1 & S==2]-Ey[g==0 & t==0 & S==2])
# abline(a=0,b=1)
# par(mfrow=c(2,2),bty="n")
# plot(gamma_1[g==1 & S==1],Ey[g==1 & t==0 & S==1]-Ey[g==1 & t==-1 & S==1])
# abline(a=0,b=1)
# plot(gamma_1[g==1 & S==2],Ey[g==1 & t==0 & S==2]-Ey[g==1 & t==-1 & S==2])
# abline(a=0,b=1)
# plot((gamma_0)[g==0 & S==1],Ey[g==0 & t==0 & S==1]-Ey[g==0 & t==-1 & S==1])
# abline(a=0,b=1)
# plot((gamma_0)[g==0 & S==2],Ey[g==0 & t==0 & S==2]-Ey[g==0 & t==-1 & S==2])
# abline(a=0,b=1)
# #####
# par(mfrow=c(2,2),bty="n")
# plot(rowMeans(bta1[g==1 & S==1,]),(beta_1+tau_1+alpha_1)[g==1 & S==1])
# abline(a=0,b=1)
# plot(rowMeans(bta1[g==1 & S==2,]),(beta_1+tau_1+alpha_1)[g==1 & S==2])
# abline(a=0,b=1)
# plot(rowMeans(b0[g==0 & S==1,]),(beta_0)[g==0 & S==1])
# abline(a=0,b=1)
# plot(rowMeans(b0[g==0 & S==2,]),(beta_0)[g==0 & S==2])
# abline(a=0,b=1)
# par(mfrow=c(2,2),bty="n")
# plot(rowMeans(g1[g==1 & S==1,]),gamma_1[g==1 & S==1])
# abline(a=0,b=1)
# plot(rowMeans(g1[g==1 & S==2,]),gamma_1[g==1 & S==2])
# abline(a=0,b=1)
# plot(rowMeans(g0[g==0 & S==1,]),gamma_0[g==0 & S==1])
# abline(a=0,b=1)
# plot(rowMeans(g0[g==0 & S==2,]),gamma_0[g==0 & S==2])
# abline(a=0,b=1)
#####
par(mfrow=c(2,2),bty="n")
plot(rowMeans(out$beta1_tau1_alpha1)[1:n][g==1],(beta_1+tau_1+alpha_1)[g==1])
abline(a=0,b=1)
plot(rowMeans(out$beta0)[1:n][g==0],(beta_0)[g==0])
abline(a=0,b=1)
plot(rowMeans(out$gamma1)[1:n][g==1],gamma_1[g==1])
abline(a=0,b=1)
plot(rowMeans(out$gamma0)[1:n][g==0],gamma_0[g==0])
abline(a=0,b=1)
#####
cat.x2 <- "large"
cat.x3 <- "c"
par(mfrow=c(2,2),bty="n")
plot(rowMeans(out$beta1_tau1_alpha1)[1:n][x2==cat.x2 & x3==cat.x3 & g==1],(beta_1+tau_1+alpha_1)[x2==cat.x2 & x3==cat.x3 & g==1])
abline(a=0,b=1)
plot(rowMeans(out$beta0)[1:n][x2==cat.x2 & x3==cat.x3 & g==0],(beta_0)[x2==cat.x2 & x3==cat.x3 & g==0])
abline(a=0,b=1)
plot(rowMeans(out$gamma1)[1:n][x2==cat.x2 & x3==cat.x3 & g==1],gamma_1[x2==cat.x2 & x3==cat.x3 & g==1])
abline(a=0,b=1)
plot(rowMeans(out$gamma0)[1:n][x2==cat.x2 & x3==cat.x3 & g==0],gamma_0[x2==cat.x2 & x3==cat.x3 & g==0])
abline(a=0,b=1)
# par(mfrow=c(1,2),mar=c(5, 4, 2, 2) + 0.1,bty="n")
# hist(colMeans(bta1[g==1 & S==1,]) - colMeans(b0[g==0 & S==1,]),xlab="Posterior draws",main=bquote(X %in% S[1]))
# abline(v=mean((tau_1+alpha_1)[g==1 & S==1]),col="firebrick",lwd=2)
# hist(colMeans(bta1[g==1 & S==2,]) - colMeans(b0[g==0 & S==2,]),xlab="Posterior draws",main=bquote(X %in% S[2]))
# abline(v=mean((tau_1+alpha_1)[g==1 & S==2]),col="firebrick",lwd=2)


bta1 <- out$beta1_tau1_alpha1[1:n,]
b0 <- out$beta0[1:n,]
g1 <- out$gamma1[1:n,]
g0 <- out$gamma0[1:n,]
eps <- 0.01
bart_catt_region <- lapply(1:draws, function(i) searchPTA::searchPTA(df,g1[,i],g0[,i],bta1[,i],b0[,i],epsilon=eps,minsplit=1,minbucket=1,cp=0))
post.cart.posterior <- matrix(NA,n,draws)
for (j in 1:n)
{
  for (i in 1:(draws))
  {
    region.pj <- bart_catt_region[[i]]$regions[j,]
    region.pj <- which(region.pj==TRUE)
    if (sum(region.pj)>0) post.cart.posterior[j,i] <- bart_catt_region[[i]]$catt[region.pj]
  }
}

num.points <- 5 ## number of points to choose
prob_id <- rowMeans(!is.na(post.cart.posterior)) ## posterior probability of being in PTA region
#### Region S1
points.s1 <- order(prob_id)
points.s1 <- points.s1[(g==1 & S==1)[points.s1]]
points.s1 <- points.s1[seq(1,length(points.s1),length.out=num.points)]
#### Region S2
points.s2 <- order(prob_id)
points.s2 <- points.s2[(g==1 & S==2)[points.s2]]
points.s2 <- points.s2[seq(1,length(points.s2),length.out=num.points)]
#### Unidentified region
points.snull <- order(prob_id)
points.snull <- points.snull[(g==1 & S==3)[points.snull]]
points.snull <- points.snull[seq(1,length(points.snull),length.out=num.points)]
### Plot
library(ggridges)
library(ggplot2)
library(gridExtra)

plotlist <- list(3)
for (j in 1:3){
  ## Select region
  if (j == 1){points <- points.s1}
  if (j == 2){points <- points.s2}
  if (j == 3){points <- points.snull}
  ## Create data frame for plot
  post.mat <- post.cart.posterior[points,]
  m <- ncol(post.mat)
  k <- nrow(post.mat)
  tempdf <- data.frame(catt = c(post.mat),obs = factor(rep(1:k,each=1)),prob_id=prob_id[points])
  if (j!=3)
  {
    ## Plots for PTA regions with vertical line indicating true CATT in region
    plotlist[[j]] <- ggplot(tempdf, aes(x = `catt`, y = `obs`, fill = prob_id)) +
      xlim(1.5, 3.85) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, na.rm = TRUE) +
      geom_vline(xintercept = mean((tau_1+alpha_1)[g==1 & S==j]),color="firebrick") +
      labs(title = 'Posterior Densities')
  } else
  {
    ## Plots for non-PTA regions not indicating true CATT in region
    plotlist[[j]] <- ggplot(tempdf, aes(x = `catt`, y = `obs`, fill = prob_id)) +
      xlim(1.5, 3.85) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, na.rm = TRUE) +
      labs(title = 'Posterior Densities')
  }
}
### Print plot
grid.arrange(plotlist[[1]],plotlist[[2]],plotlist[[3]],nrow = 3)
