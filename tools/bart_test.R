## Setup-----------------------------------------------------------------------
# devtools::load_all()
# library(searchPTA)
# seed <- 007
n <- 10000
# set.seed(seed)
## Toggles for which case to run
continuous <- FALSE
discrete <- TRUE
## Continuous X----------------------------------------------------------------
# For this example, we generate two continuous features X_1,X_2 \in [0,1]
# We define the PTA regions as:
# S_1: X_1 < 0.5 & X_2 < 0.5
# S_2: X_2 > 0.5
# We make the points in region S_1 more dispersed while the points in region S_2
# are more concentrated near the center of that region
# This is to highlight that, with continuous X, having many 'border' points can
# cause trouble: if beta_1 and beta_0 are continuous, points near the PTA regions
# will have beta_1-beta_0 \approx 0, so the risk of misclassification is higher,
# implying that our predictions might borrow too much from unidentified regions
if (continuous)
{
  ### Generate data--------------------------------------------------------------
  #### Fixed features
  t <- c(rep(-1,n),rep(0,n),rep(1,n))
  n1 <- n %/% 2
  n0 <- n-n1
  g <- c(rep(1,n1),rep(0,n0))
  z <- c(g,g,g)*(t==1)
  #### Generate X
  sig.x <- 0.02
  s1.g1 <- cbind(rnorm(n,0.2,sig.x),rnorm(n,0.2,sig.x))
  s2.g1 <- cbind(rnorm(n,0.45,0.05),rnorm(n,0.7,sig.x))
  snull.g1 <- cbind(rnorm(n,0.7,sig.x),rnorm(n,0.2,sig.x))
  s1.g0 <- cbind(rnorm(n,0.3,sig.x),rnorm(n,0.3,sig.x))
  s2.g0 <- cbind(rnorm(n,0.55,0.05),rnorm(n,0.8,sig.x))
  snull.g0 <- cbind(rnorm(n,0.8,sig.x),rnorm(n,0.3,sig.x))
  regions <- sample(1:3,n,replace=TRUE)
  x1 <- (regions==1 & g==1)*punif(s1.g1[,1]) + (regions==2 & g==1)*punif(s2.g1[,1]) + (regions==3 & g==1)*punif(snull.g1[,1]) + (regions==1 & g==0)*punif(s1.g0[,1]) + (regions==2 & g==0)*punif(s2.g0[,1]) + (regions==3 & g==0)*punif(snull.g0[,1])
  x2 <- (regions==1 & g==1)*punif(s1.g1[,2]) + (regions==2 & g==1)*punif(s2.g1[,2]) + (regions==3 & g==1)*punif(snull.g1[,2]) + (regions==1 & g==0)*punif(s1.g0[,2]) + (regions==2 & g==0)*punif(s2.g0[,2]) + (regions==3 & g==0)*punif(snull.g0[,2])
  x <- rep(2,n)
  for (i in 1:n)
  {
    if (x1[i]>0.5 & x2[i]<0.5) x[i] <- 1
    if (x1[i]>0.5 & x2[i]>0.5) x[i] <- 3
    if (x1[i]<0.5 & x2[i]>0.5) x[i] <- 4
  }
  #### Define functions for each group
  ##### 1) mu
  mu_1 <- 0.7*sin(1.5*pi*x1) + exp(x2)
  mu_0 <- 0.8*log(0.5+x1) + pnorm(x2,0.5,0.25)
  ##### 2) beta
  beta_1 <- cos(0.5*pi*x1) + 3*log(x2+0.5)
  beta_0 <- (x2+0.5)^3 - (x1+0.5)^2
  ##### Demeaning beta in the right regions
  beta_1_avg_2 <- mean(beta_1[g==1 & x==2])
  beta_0_avg_2 <- mean(beta_0[g==0 & x==2])
  beta_1_avg_34 <- mean(beta_1[g==1 & x %in% 3:4])
  beta_0_avg_34 <- mean(beta_0[g==0 & x %in% 3:4])
  beta_1 <- beta_1 - (x==2 & g==1)*(beta_1_avg_2-beta_0_avg_2) - (x %in% 3:4)*(beta_1_avg_34-beta_0_avg_34)
  ##### Checking that betas average to 0 in the right regions
  mean(beta_1[g==1 & x==2])-mean(beta_0[g==0 & x==2])
  mean(beta_1[g==1 & x %in% 3:4])-mean(beta_0[g==0 & x %in% 3:4])
  ##### 3) gamma
  rho <- 2
  gamma_1 <- beta_1*rho
  gamma_0 <- beta_0*rho
  ##### Checking that gammas average to 0 in the right regions
  mean(gamma_1[g==1 & x==2])-mean(gamma_0[g==0 & x==2])
  mean(gamma_1[g==1 & x %in% 3:4])-mean(gamma_0[g==0 & x %in% 3:4])
  ##### 4) alpha
  alpha_1 <- sqrt(x1) + x2
  alpha_0 <- 2/(5+0.2*exp(5*x1)) + 0.7*(x2+0.3)^2
  ##### 5) tau
  tau_1 <- pnorm(x1,0.5,0.2) + sin(pi*x2)
  tau_0 <- 0.4*x1^2 + 0.25*cos(0.75*pi*x2)
  ##### Scaling tau and alpha
  tau_1 <- tau_1/max(tau_1)*mean(mu_1)*0.5
  tau_0 <- tau_0/max(tau_0)*mean(mu_0)*0.5
  alpha_1 <- alpha_1/max(alpha_1)*mean(mu_1)*0.5
  alpha_0 <- alpha_0/max(alpha_0)*mean(mu_0)*0.5
  ### Perform CART search with true functions------------------------------------
  bta1 <- beta_1+tau_1+alpha_1
  ## Set epsilon value; for this, we look at quantiles of the gamma_1-gamma_0 distribution
  ### We use the values of gamma here because beta_1 is not identified
  epsilon <- quantile(abs(gamma_1-gamma_0),0)
  ## Perform CART search
  X <- data.frame(X1=x1,X2=x2,g=g)
  # X <- data.frame(X=factor(x,ordered=TRUE),g=g)
  catt_region <- searchPTA::searchPTA(x=X,delta1_aux=gamma_1,delta0_aux=gamma_0,delta1_main=bta1,delta0_main=beta_0,epsilon=epsilon,minsplit=1,minbucket=1,cp=0)
  ### Plot regions
  col <- rep(NA,nrow(catt_region$regions))
  for (i in 1:length(col))
  {
    temp <- which(catt_region$regions[i,])
    if (length(temp)==0) temp <- 1
    else temp <- temp+1
    col[i] <- temp
  }
  plot(x1,x2,type="n",cex.axis=0.7,xlab=bquote(X[1]),ylab=bquote(X[2]),xlim=c(0,1),ylim=c(0,1))
  abline(v=0.5,h=0.5,lwd=0.5)
  points(x1,x2,bg=col,cex=0.7,cex.axis=0.7,pch=21,xlab=bquote(X[1]),ylab=bquote(X[2]))
}
## Discrete X------------------------------------------------------------------
if (discrete)
{
  ### Generate data--------------------------------------------------------------
  #### Fixed features
  t <- c(rep(-1,n),rep(0,n),rep(1,n))
  n1 <- n %/% 2
  n0 <- n-n1
  g <- c(rep(1,n1),rep(0,n0))
  z <- c(g,g,g)*(t==1)
  #### Generate X
  sig.x <- 0.1
  s1.g1 <- cbind(rnorm(n,0.2,sig.x),rnorm(n,0.2,sig.x))
  s2.g1 <- cbind(rnorm(n,0.45,0.05),rnorm(n,0.7,sig.x))
  snull.g1 <- cbind(rnorm(n,0.7,sig.x),rnorm(n,0.2,sig.x))
  s1.g0 <- cbind(rnorm(n,0.3,sig.x),rnorm(n,0.3,sig.x))
  s2.g0 <- cbind(rnorm(n,0.55,0.05),rnorm(n,0.8,sig.x))
  snull.g0 <- cbind(rnorm(n,0.8,sig.x),rnorm(n,0.3,sig.x))
  regions <- sample(1:3,n,replace=TRUE)
  x1 <- (regions==1 & g==1)*punif(s1.g1[,1]) + (regions==2 & g==1)*punif(s2.g1[,1]) + (regions==3 & g==1)*punif(snull.g1[,1]) + (regions==1 & g==0)*punif(s1.g0[,1]) + (regions==2 & g==0)*punif(s2.g0[,1]) + (regions==3 & g==0)*punif(snull.g0[,1])
  x2 <- (regions==1 & g==1)*punif(s1.g1[,2]) + (regions==2 & g==1)*punif(s2.g1[,2]) + (regions==3 & g==1)*punif(snull.g1[,2]) + (regions==1 & g==0)*punif(s1.g0[,2]) + (regions==2 & g==0)*punif(s2.g0[,2]) + (regions==3 & g==0)*punif(snull.g0[,2])
  x1.raw <- x1
  x2.raw <- x2
  ncat <- 3
  x1 <- (floor(x1*ncat))/ncat
  x2 <- (floor(x2*ncat))/ncat
  # x1.round <- round(x1,1)
  # x2.round <- round(x2,1)
  # x1 <- floor(x1*10)/10
  # x2 <- floor(x2*10)/10
  x <- rep(2,n)
  for (i in 1:n)
  {
    # if (x1[i]>0.5 & x2[i]<0.5) x[i] <- 1
    # if (x1[i]>0.5 & x2[i]>=0.5) x[i] <- 3
    # if (x1[i]<=0.5 & x2[i]>=0.5) x[i] <- 4
    #####
    # if (x1[i]==0.5 & x2[i]==0) x[i] <- 1
    # if (x1[i]==0.5 & x2[i]==0.5) x[i] <- 3
    # if (x1[i]==0 & x2[i]==0.5) x[i] <- 4
    #####
    if (x1[i]>0.5 & x2[i]<0.5) x[i] <- 1
    if (x1[i]<0.5 & x2[i]>=0.5) x[i] <- 3
    if (x1[i]>=0.5 & x2[i]>=0.5) x[i] <- 4
  }
  # x1 <- floor(x1/0.25)
  # x2 <- floor(x2/0.25)
  #### Define functions for each group
  ##### 1) mu
  mu_1 <- 0.7*sin(1.5*pi*x1) + exp(x2)
  mu_0 <- 0.8*log(0.5+x1) + pnorm(x2,0.5,0.25)
  ##### 2) beta
  beta_1 <- cos(0.5*pi*x1) + 3*log(x2+0.5)
  beta_0 <- (x2+0.5)^3 - (x1+0.5)^2
  ##### Demeaning beta in the right regions
  beta_1_avg_2 <- mean(beta_1[g==1 & x==2])
  beta_0_avg_2 <- mean(beta_0[g==0 & x==2])
  beta_1_avg_34 <- mean(beta_1[g==1 & x %in% 3:4])
  beta_0_avg_34 <- mean(beta_0[g==0 & x %in% 3:4])
  beta_1 <- beta_1 - (x==2 & g==1)*(beta_1_avg_2-beta_0_avg_2) - (x %in% 3:4)*(beta_1_avg_34-beta_0_avg_34)
  ##### Checking that betas average to 0 in the right regions
  mean(beta_1[g==1 & x==2])-mean(beta_0[g==0 & x==2])
  mean(beta_1[g==1 & x %in% 3:4])-mean(beta_0[g==0 & x %in% 3:4])
  ##### 3) gamma
  rho <- 2
  gamma_1 <- beta_1*rho
  gamma_0 <- beta_0*rho
  ##### Checking that gammas average to 0 in the right regions
  mean(gamma_1[g==1 & x==2])-mean(gamma_0[g==0 & x==2])
  mean(gamma_1[g==1 & x %in% 3:4])-mean(gamma_0[g==0 & x %in% 3:4])
  ##### 4) alpha
  alpha_1 <- sqrt(x1) + x2
  alpha_0 <- 2/(5+0.2*exp(5*x1)) + 0.7*(x2+0.3)^2
  ##### 5) tau
  tau_1 <- pnorm(x1,0.5,0.2) + sin(pi*x2)
  tau_0 <- 0.4*x1^2 + 0.25*cos(0.75*pi*x2)
  ##### Scaling tau and alpha
  tau_1 <- tau_1/max(tau_1)*mean(mu_1)*0.5
  tau_0 <- tau_0/max(tau_0)*mean(mu_0)*0.5
  alpha_1 <- alpha_1/max(alpha_1)*mean(mu_1)*0.5
  alpha_0 <- alpha_0/max(alpha_0)*mean(mu_0)*0.5
  ### Perform CART search with true functions------------------------------------
  bta1 <- beta_1+tau_1+alpha_1
  ## Set epsilon value; for this, we look at quantiles of the gamma_1-gamma_0 distribution
  ### We use the values of gamma here because beta_1 is not identified
  epsilon <- quantile(abs(gamma_1-gamma_0),0)
  # epsilon <- 0.05
  ## Perform CART search
  par(mfrow=c(2,2))
  for (i in 1:4)
  {
    # i <- ifelse(j%%2==0,2,4)
    if (i==1) {X <- data.frame(X1=factor(x1,ordered=FALSE),X2=factor(x2,ordered=FALSE),g=g); title <- "rounded X, unordered factor"}
    if (i==2) {X <- data.frame(X1=factor(x1,ordered=TRUE),X2=factor(x2,ordered=TRUE),g=g); title <- "rounded X, ordered factor"}
    if (i==3) {X <- data.frame(X1=x1.raw,X2=x2.raw,g=g); title <- "raw continuous X"}
    if (i==4) {X <- data.frame(X1=x1,X2=x2,g=g); title <- "rounded X, numeric"}
    if (i==5) {X <- data.frame(X=factor(x,ordered=FALSE),g=g); title <- "true regions, unordered factor"}
    if (i==6) {X <- data.frame(X=factor(x,ordered=TRUE),g=g); title <- "true regions, ordered factor"}
    catt_region <- searchPTA::searchPTA(x=X,delta1_aux=gamma_1,delta0_aux=gamma_0,delta1_main=bta1,delta0_main=beta_0,epsilon=epsilon,minsplit=1,minbucket=1,cp=0)
    if (i==1) CART <- catt_region$cart
    ### Plot regions
    col <- rep(NA,nrow(catt_region$regions))
    for (i in 1:length(col))
    {
      temp <- which(catt_region$regions[i,])
      if (length(temp)==0) temp <- 1
      else temp <- temp+1
      col[i] <- temp
    }
    shapes <- rep(21,n)
    for (i in 1:n)
    {
      if (x[i]==2) shapes[i] <- 22
      if (x[i] %in% 3:4) shapes[i] <- 24
    }
    plot(x1,x2,type="n",cex.axis=0.7,xlab=bquote(X[1]),ylab=bquote(X[2]),xlim=c(-0.05,0.05)+quantile(x1,c(0,1)),ylim=c(-0.05,0.05)+quantile(x2,c(0,1)),main=title,cex.main=0.7)
    points(x1,x2,bg=col,cex=0.7,cex.axis=0.7,pch=shapes,xlab=bquote(X[1]),ylab=bquote(X[2]))
  }
  # par(mfrow=c(1,1))
  # rpart.plot::rpart.plot(CART)
}
