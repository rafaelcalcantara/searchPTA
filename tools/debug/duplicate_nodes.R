###############################################################################
## Problem: in some situations, when X is continuous, the pta.nodes function
## assigns two nodes to the same point; it is supposed to determine the most
## disaggregated PTA regions as possible, so each point should be assigned to
## *at most* one node
###############################################################################
devtools::install_github("rafaelcalcantara/searchPTA")
library(searchPTA)
epsilon <- 0.05
### Function to fit S-Learner BART with DiD regression in the leaves
fit.bart.did <- function(df_train,df_test,features,num_trees=50,max_depth=20)
{
  out <- list(b1=NULL,b0=NULL)
  df.x.train <- subset(df_train,select=features)
  df.x.test <- subset(df_test,select=features)
  ## Set G to (0,1)
  df_train$g <- df_train$g - min(df_train$g)
  df_test$g <- df_test$g - min(df_test$g)
  ## Lists of parameters for the Stochtree BART function
  bart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
  bart.mean.parmlist <- list(num_trees=num_trees, min_samples_leaf=20, alpha=0.95, beta=2,
                             max_depth=max_depth, sample_sigma2_leaf=FALSE, sigma2_leaf_init = diag(rep(0.1/150,4)))
  ## Set basis vector for leaf regressions
  Psi <- cbind(1-df_train$g,df_train$g,(1-df_train$g)*df_train$t,df_train$g*df_train$t)
  Psi_test <- cbind(1-df_test$g,df_test$g,(1-df_test$g)*df_test$t,df_test$g*df_test$t)
  ## Model fit
  bart.fit = stochtree::bart(X_train=df.x.train, y_train=df_train$y,
                             leaf_basis_train = Psi, mean_forest_params=bart.mean.parmlist,
                             general_params=bart.global.parmlist,
                             num_mcmc=500,num_gfr=50)
  ## Extract beta for g=1 and g=0
  ### We use the 'predict_raw' function from the bart.fit object
  ### 1) Pre-process test data to the format expected by the predict_raw function
  X_test <- stochtree::preprocessPredictionData(df.x.test,bart.fit$train_set_metadata)
  ### 2) Create 'ForestDataset' structure to be read by the predict_raw function
  X_test <- stochtree::createForestDataset(X_test,basis=Psi_test)
  ### 3) Use predict_raw to extract betas for each group. Because the basis is (1,g,t,g*t), we need the 3rd and 4th coefficients
  betas <- bart.fit$mean_forests$predict_raw(X_test)
  ### 4) Save results; the predict_raw output is a (n,p,m) array, where n=num_obs, p is the length of the basis vector, m is the number of MCMC samples
  #### Extracting betas[,3,] gives us a n*m matrix of individual draws for the coefficient on t, and similarly betas[,4,] for g*t
  #### The raw predictions need to be scaled back to match the actual data (i.e. the raw predictions are based on (Y-Y_bar)/sd(Y))
  #### Note that, because those two coefficients represent differences, the Y_bar term cancels out and we only need to multiply by sd(Y):
  #### E.g. the coefficient on t is equal to E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=1]-E[(Y-Y_bar)/sd(Y)|X=x,G=0,T=0] = (E[Y|X=x,G=0,T=1]-E[Y|X=x,G=0,T=0])/sd(Y)
  out$b1 <- betas[,4,]*sd(df_train$y)
  out$b0 <- betas[,3,]*sd(df_train$y)
  return(out)
}
# 1) First, we check the functions with discrete data--------------------------
## Generate data
### Define parameters for each group
mu_1 <- c(1,3,1,3)
mu_0 <- c(1,2,3,4)
mu <- cbind(mu_0,mu_1)
#########
beta_1 <- c(1,1,1,1)
beta_0 <- c(7,1,2,0)
beta <- cbind(beta_0,beta_1)
#########
alpha_1 <- seq(0,1,length.out = 4)
alpha_0 <- seq(1,0,length.out = 4)
alpha <- cbind(alpha_0, alpha_1)
#########
tau_1 <- c(1,2,3,4)
tau_0 <- c(2,1,4,3)
tau <- cbind(tau_0, tau_1)
### Generate features
n <- 1500
t <- c(rep(-1,n/3),rep(0,n/3),rep(1,n/3))
g <- rep(c(0,1),n/2)
z <- g*(t==1)
x <- sample(1:4,n/3,replace=TRUE)
x <- c(x,x,x)
### Generate outcome
idx <- cbind(x,g+1)
Ey <- mu[idx] + beta[idx]*t + tau[idx]*z + alpha[idx]*z*t
y <- Ey + 0.5*sd(Ey)*rnorm(n)
df <- data.frame(y,t,g,x=x)
df_main <- subset(df,t>=0)
df_placebo <- subset(df,t<=0)
### Checking the functions for true beta_1 and beta_0
b1 <- beta_1[x]
b0 <- beta_0[x]
#### placebo.cart
cart <- placebo.cart(data.frame(x=as.factor(x[t<=0])),b1[t<=0],b0[t<=0],epsilon)
rpart.plot::rpart.plot(cart)
### pta.nodes
nodes <- pta.nodes(cart,epsilon)
table(x[t<=0][nodes[,1]])
table(x[t<=0][nodes[,2]])
table(x[t<=0])
### catt.per.region
catt <- catt.per.region(b1[t>=0],b0[t>=0],nodes)
catt[x[t>=0] == 2]
catt[x[t>=0] %in% c(3,4)]
