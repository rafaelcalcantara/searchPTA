###############################################################################
## Benchmarking functions time and memory-wise
###############################################################################
n_vec <- seq(1000,10000,length.out=10)
time.list <- vector("list",length=length(n_vec))
for (N in 1:length(n_vec))
{
  n <- n_vec[N]
  n_regions <- 2
  if (n_regions==2) lims <- c(-1,5)
  if (n_regions==3) lims <- c(-0.5,3.2)
  ## Generate data---------------------------------------------------------------
  ### Fixed features
  t <- c(rep(-1,n),rep(0,n),rep(1,n))
  n1 <- n %/% 2
  n0 <- n-n1
  g <- c(rep(1,n1),rep(0,n0))
  z <- c(g,g,g)*(t==1)
  x1 <- g*rgamma(n,3,3) + (1-g)*rgamma(n,6,3)
  x2 <- ifelse(g==1,sample(c("small","large"),n,replace=TRUE,prob=c(0.6,0.4)),sample(c("small","large"),n,replace=TRUE,prob=c(0.4,0.6)))
  x3 <- ifelse(g==1,sample(c("a","b","c"),n,replace=TRUE,prob=c(0.4,0.25,0.35)),sample(c("a","b","c"),n,replace=TRUE,prob=rep(1/3,3)))
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
    out <- rep(NA,length(g))
    for (i in 1:length(g))
    {
      if (x2[i] == "large" & g[i] == 1)
      {
        a <- 1; b <- 0.2
      } else if (x2[i] == "large" & g[i] == 0)
      {
        a <- 0.5; b <- 0.6
      } else if (x2[i] == "small" & x3[i] %in% c("a","b") & g[i] == 1)
      {
        a <- 3.5; b <- 0.3
      } else if (x2[i] == "small" & x3[i] %in% c("a","b") & g[i] == 0)
      {
        a <- 0.7; b <- 0.1
      } else if (x2[i] == "small" & x3[i] == "c" & g[i] == 1)
      {
        a <- 1.5; b <- 0.8
      } else if (x2[i] == "small" & x3[i] == "c" & g[i] == 0)
      {
        a <- 0.5; b <- 0.5
      }
      out[i] <- a + b*x1[i]
    }
    return(out)
  }
  gamma_fun <- function(x1,x2,x3,g) beta_fun(x1,x2,x3,g)
  tau_fun <- function(x2,x3)
  {
    a <- ifelse(x2=="large",0.05,ifelse(x3=="c",0.3,0.7))
    return(a)
  }
  alpha_fun <- function(x2)
  {
    a0 <- ifelse(x2=="large" & x3=="a",0.25,0.4)
    return(a0)
  }
  ### Calculate function values
  mu <- mu_fun(x1,x2,g)
  beta <- beta_fun(x1,x2,x3,g)
  gamma <- gamma_fun(x1,x2,x3,g)
  tau <- tau_fun(x2,x3)
  alpha <- alpha_fun(x2)
  catt <- (tau+alpha)[g==1]
  ### Define PTA regions
  S <- rep(3,n)
  for (i in 1:n)
  {
    if (x2[i] == "large") S[i] <- 1
    if (x2[i] == "small" & x3[i] %in% c("a","b")) S[i] <- 2
  }
  ### Establish PTA in regions
  beta_avg_1_1 <- mean(beta[S==1 & g==1])
  beta_avg_1_0 <- mean(beta[S==1 & g==0])
  beta_avg_2_1 <- mean(beta[S==2 & g==1])
  beta_avg_2_0 <- mean(beta[S==2 & g==0])
  beta_avg_3_1 <- mean(beta[S==3 & g==1])
  beta_avg_3_0 <- mean(beta[S==3 & g==0])
  gamma_avg_1_1 <- mean(gamma[S==1 & g==1])
  gamma_avg_1_0 <- mean(gamma[S==1 & g==0])
  gamma_avg_2_1 <- mean(gamma[S==2 & g==1])
  gamma_avg_2_0 <- mean(gamma[S==2 & g==0])
  gamma_avg_3_1 <- mean(gamma[S==3 & g==1])
  gamma_avg_3_0 <- mean(gamma[S==3 & g==0])
  sig.beta <- 0
  if (n_regions==2)
  {
    beta <- beta - (S==1 & g==1)*(beta_avg_1_1-beta_avg_1_0) - (S==2 & g==1)*(beta_avg_2_1-beta_avg_2_0) - (S==3 & g==1)*(beta_avg_3_1-beta_avg_3_0 - 3)
    gamma <- gamma - (S==1 & g==1)*(gamma_avg_1_1-gamma_avg_1_0) - (S==2 & g==1)*(gamma_avg_2_1-gamma_avg_2_0) - (S==3 & g==1)*(gamma_avg_3_1-gamma_avg_3_0 - 3)
  } else if (n_regions==3)
  {
    beta <- beta - (S==1 & g==1)*(beta_avg_1_1-beta_avg_1_0+rnorm(n,0,sig.beta)) - (S==2 & g==1)*(beta_avg_2_1-beta_avg_2_0) - (S==3 & g==1)*(beta_avg_3_1-beta_avg_3_0)
    gamma <- gamma - (S==1 & g==1)*(gamma_avg_1_1-gamma_avg_1_0+rnorm(n,0,sig.beta)) - (S==2 & g==1)*(gamma_avg_2_1-gamma_avg_2_0) - (S==3 & g==1)*(gamma_avg_3_1-gamma_avg_3_0)
  }
  ### Generate outcome
  Ey <- rep(mu,3) - rep(gamma,3)*(t==-1) + rep(beta,3)*(t==1) +  rep(tau,3)*z + rep(alpha,3)*z*(t==1)
  error.sd <- 0.2*sd(Ey)
  y <- Ey + rnorm(n,0,error.sd)
  df <- data.frame(X1=x1,X2=factor(x2,ordered=TRUE),X3=factor(x3,ordered=FALSE),g=g)
  ######
  eps <- 0.05
  m <- 3
  time <- rep(0,m)
  for (i in 1:m)
  {
    t0 <- Sys.time()
    test <- searchPTA::searchPTA(df,gamma,gamma,beta+tau+alpha,beta,epsilon=eps,minsplit=1,minbucket=1,cp=0,maxdepth=30,saveCART = FALSE)
    t1 <- Sys.time()
    time[i] <- as.numeric(t1-t0)
  }
  time.list[[N]] <- time
}
## Check times
plot(n_vec/1000,sapply(time.list,mean),type="b",xlab="N (in 1k)",ylab="Time (sec)")
abline(a=0,b=1,lty=2,col="gray")
