## Setup-----------------------------------------------------------------------
seed <- 007
set.seed(seed)
n <- 1000
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
  return(a+b*x1)
}
beta_fun <- function(x1,x2,x3,g)
{
  ## Intercept
  a0 <- ifelse(x2=="large",0.4,0.1)
  a1 <- ifelse(x3=="a",0.6, ifelse(x3=="b",0.45,0.15))
  a2 <- ifelse(g==1,0.1,0.3)
  a <- a0+a1+a2
  ## x1 coefficient
  b <- ifelse(g==1,0.4,0.1)
  ## Return beta
  return(a+b*x1)
}
gamma_fun <- function(x1,x2,x3,g) beta_fun(x1,x2,x3,g)
tau_fun <- function(x2,x3)
{
  a0 <- ifelse(x2=="large",0.1,0.8)
  a1 <- ifelse(x3=="a",0.3, 0.1)
  return(a0+a1)
}
alpha_fun <- function(x2)
{
  a0 <- ifelse(x2=="large",0,0.4)
  # a1 <- ifelse(x3=="c",0.5,0)
  return(a0)
}
### Calculate function values
mu <- mu_fun(x1,x2,g)
beta <- beta_fun(x1,x2,x3,g)
gamma <- gamma_fun(x1,x2,x3,g)
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
beta <- beta - (S==1 & g==1)*(beta_avg_1_1-beta_avg_1_0) - (S==2 & g==1)*(beta_avg_2_1-beta_avg_2_0) - (S==3 & g==1)*(beta_avg_3_1-beta_avg_3_0 - 3)
gamma_avg_1_1 <- mean(gamma[S==1 & g==1])
gamma_avg_1_0 <- mean(gamma[S==1 & g==0])
gamma_avg_2_1 <- mean(gamma[S==2 & g==1])
gamma_avg_2_0 <- mean(gamma[S==2 & g==0])
gamma_avg_3_1 <- mean(gamma[S==3 & g==1])
gamma_avg_3_0 <- mean(gamma[S==3 & g==0])
gamma <- gamma - (S==1 & g==1)*(gamma_avg_1_1-gamma_avg_1_0) - (S==2 & g==1)*(gamma_avg_2_1-gamma_avg_2_0) - (S==3 & g==1)*(gamma_avg_3_1-gamma_avg_3_0 - 3)
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
library(stochtree)
library(foreach)
library(doParallel)
### First stage: fit BART to obtain (\gamma_1, \gamma_0), (\beta_1, \beta_0) predictions
#### BART setup
num_chains <- 20
draws <- 50
burnin <- 10
keep.draws <- (burnin+1):(draws+burnin)
ncores <- 10
num_trees <- 100
max_depth <- 20
bart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
bart.mean.parmlist <- list(num_trees=num_trees, min_samples_leaf=20, alpha=0.95, beta=2,
                           max_depth=max_depth, sample_sigma2_leaf=FALSE,
                           sigma2_leaf_init = diag(rep(0.1/num_trees,6)))
#### Organize data
model.df <- data.frame(y,t,g=c(g,g,g),x1=c(x1,x1,x1),x2=c(x2,x2,x2),x3=c(x3,x3,x3))
model.df$x2 <- factor(model.df$x2,ordered=TRUE)
model.df$x3 <- factor(model.df$x3,ordered=FALSE)
df.x <- subset(model.df,select=4:ncol(model.df))
model.df$pre <- -as.numeric(model.df$t==-1)
model.df$post <- as.numeric(model.df$t==1)
#### Set basis vector for leaf regressions
Psi <- cbind(1-model.df$g,model.df$g,(1-model.df$g)*model.df$pre,model.df$g*model.df$pre,(1-model.df$g)*model.df$post,model.df$g*model.df$post)
#### Run model
cl <- makeCluster(ncores)
registerDoParallel(cl)
bart_model_outputs <- foreach (i = 1:num_chains) %dopar% {
  out <- list(gamma1=NULL,gamma0=NULL,beta1=NULL,beta0=NULL)
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
  out
}
stopCluster(cl)
#### Join results from each run
bart_catt_i <- do.call("cbind",lapply(bart_model_outputs, function(i) (i$beta1_tau1_alpha1-i$beta0)[1:n,]))
bta1 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$beta1_tau1_alpha1[1:n,]))
b0 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$beta0[1:n,]))
g1 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$gamma1[1:n,]))
g0 <- do.call("cbind",lapply(bart_model_outputs, function(i) i$gamma0[1:n,]))
#### Check the quality of the BART estimates
par(mfrow=c(2,2),bty="n")
plot(rowMeans(bta1)[1:n][g==1],(beta+tau+alpha)[g==1])
abline(a=0,b=1)
plot(rowMeans(b0)[1:n][g==0],(beta)[g==0])
abline(a=0,b=1)
plot(rowMeans(g1)[1:n][g==1],gamma[g==1])
abline(a=0,b=1)
plot(rowMeans(g0)[1:n][g==0],gamma[g==0])
abline(a=0,b=1)
### Second stage: apply our CART search procedure to the posterior draws of BART
eps <- 0.03
df <- data.frame(X1=x1,X2=factor(x2,ordered=TRUE),X3=factor(x3,ordered=FALSE),g=g)
# df <- data.frame(X=S,g=g)
cl <- makeCluster(ncores)
registerDoParallel(cl)
bart_catt_region <- foreach(i = 1:(draws*num_chains)) %dopar%
{
  searchPTA::searchPTA(df,g1[,i],g0[,i],bta1[,i],b0[,i],saveCART = FALSE,epsilon=eps,minsplit=1,minbucket=1,cp=0,maxdepth=30)
}
stopCluster(cl)
#### Calculate CATT posterior
post.cart.posterior <- sapply(bart_catt_region, function(i) i$catt)
#### Check distribution of posterior probability of identification for points within each region
prob_id <- rowMeans(!is.na(post.cart.posterior)) ## posterior probability of being in PTA region
par(mfrow=c(1,1),bty="n",cex.lab=0.7)
boxplot(prob_id~S,ylab="Posterior prob. of point being in identified region")
### Ridge plot for posterior distribution of a handful of points
library(ggridges)
library(ggplot2)
library(gridExtra)
num.points <- 10 ## number of points to choose
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
#### Plot
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
      xlim(0, 2) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, na.rm = TRUE) +
      geom_vline(xintercept = mean((tau+alpha)[g==1 & S==j]),color="firebrick") +
      labs(title = 'Posterior Densities')
  } else
  {
    ## Plots for non-PTA regions not indicating true CATT in region
    plotlist[[j]] <- ggplot(tempdf, aes(x = `catt`, y = `obs`, fill = prob_id)) +
      xlim(0, 2) +
      geom_density_ridges(scale = 3, rel_min_height = 0.01, na.rm = TRUE) +
      labs(title = 'Posterior Densities')
  }
}
##### Print plot
grid.arrange(plotlist[[1]],plotlist[[2]],plotlist[[3]],nrow = 3)
