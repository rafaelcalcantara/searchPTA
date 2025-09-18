###############################################################################
## Illustrations of our procedure
###############################################################################
# devtools::load_all()
epsilon <- 0.01
cp <- 0.01
minsplit <- 20
n <- 15000
dependent <- FALSE ## If X depends on G
## Generate data
### Define parameters for each group
mu_1 <- c(1,3,1,3)
mu_0 <- c(1,2,3,4)
mu <- cbind(mu_0,mu_1)
#########
alpha_1 <- seq(0,1,length.out = 4)
alpha_0 <- seq(1,0,length.out = 4)
alpha <- cbind(alpha_0, alpha_1)
#########
tau_1 <- c(1,2,3,4)
tau_0 <- c(2,1,4,3)
tau <- cbind(tau_0, tau_1)
### Generate features
t <- c(rep(-1,n/3),rep(0,n/3),rep(1,n/3))
g <- rep(c(0,1),n/2)
z <- g*(t==1)
#########
if (dependent)
{
  x1 <- runif(n/3)*g[1:(n/3)] + rbeta(n/3,3,1)*(1-g[1:(n/3)])
  x2 <- runif(n/3)*g[1:(n/3)] + rbeta(n/3,1,3)*(1-g[1:(n/3)])
  ## Quadrant probabilities
  p1.g1 <- p2.g1 <- p3.g1 <- p4.g1 <- 0.25
  p1.g0 <- pbeta(0.5,3,1)*pbeta(0.5,1,3)
  p2.g0 <- (1-pbeta(0.5,3,1))*pbeta(0.5,1,3)
  p3.g0 <- (1-pbeta(0.5,3,1))*(1-pbeta(0.5,1,3))
  p4.g0 <- pbeta(0.5,3,1)*(1-pbeta(0.5,1,3))
  #########
  beta_1 <- c(1/p1.g1,1,1*p3.g1/(p3.g1+p4.g1),1.01*p4.g1/(p3.g1+p4.g1))
  beta_0 <- c(7/p1.g0,1,2*p3.g0/(p3.g0+p4.g0),-5.96*p4.g0/(p3.g0+p4.g0))
  beta <- cbind(beta_0,beta_1)
} else
{
  ## Generate data
  x1 <- runif(n/3)
  x2 <- runif(n/3)
  #########
  beta_1 <- c(1,1,1,1)*0.1
  beta_0 <- c(7,1,2,0)*0.1
  beta <- cbind(beta_0,beta_1)
}
####
x <- rep(1,n/3)
for (i in 1:(n/3))
{
  if (x1[i]>0.5 & x2[i]<0.5) x[i] <- 2
  if (x1[i]>0.5 & x2[i]>0.5) x[i] <- 3
  if (x1[i]<0.5 & x2[i]>0.5) x[i] <- 4
}
x <- c(x,x,x)
idx <- cbind(x,g+1)
### Generate outcome
Ey <- mu[idx] + beta[idx]*t + tau[idx]*z + alpha[idx]*z*t
y <- Ey + 0.5*sd(Ey)*rnorm(n)
df <- data.frame(y,t,g,x1=c(x1,x1,x1),x2=c(x2,x2,x2))
# df <- data.frame(y,t,g,x1=x)
df_main <- subset(df,t>=0)
df_placebo <- subset(df,t<=0)
### Checking the functions for true beta_1 and beta_0
b1 <- beta_1[x][t==0]
b0 <- beta_0[x][t==0]
alpha1 <- alpha_1[x][t==0]
tau1 <- tau_1[x][t==0]
test <- searchPTA::searchPTA(x=subset(df_main,t==1,select=c("x1","x2","g")),delta1_aux=b1,delta0_aux=b0,delta1_main=b1+alpha1+tau1,delta0_main=b0,epsilon=epsilon,saveCART = TRUE,cp=cp,minsplit=minsplit)
######
print("Feature summary per region")
lapply(test$x.per.regions,summary)
print("True average trend difference per region")
apply(test$regions,2, function(i) mean(b1[i & g[t==0]==1])-mean(b0[i & g[t==0]==0]))
print("True CATT per region")
apply(test$regions,2, function(i) mean(tau1[i & g[t==0]==1]+alpha1[i & g[t==0]==1]))
print("Estimated CATT per region")
test$catt
####
n_regions <- ncol(test$regions)
col <- seq(2,length.out=n_regions)
bg <- scales::alpha(col,0.5)
col.labels <- lapply(1:n_regions,
                     function(i) bquote(X[1] %in%~"("~.(round(min(x1[test$regions[1:(n/3),i]]),3))~","~.(round(max(x1[test$regions[1:(n/3),i]]),3))~")"~";"~X[2] %in% ~"("~.(round(min(x2[test$regions[1:(n/3),i]]),3))~","~.(round(max(x2[test$regions[1:(n/3),i]]),3))~")"))
col.labels[[length(col.labels)+1]] <- "NA"
####
layout(matrix(c(1,2),ncol=1),heights = c(3,1))
par(mar=c(3,4,1,1),mgp=c(2, 0.5, 0))
plot(x1,x2,type="n",bty="n",xlab=bquote(X[1]),ylab=bquote(X[2]),cex.axis=0.7,cex.lab=0.7)
polygon(x=c(0.5,0.5,1,1),y=c(0,0.5,0.5,0))
polygon(x=c(0,0,1,1),y=c(0.5,1,1,0.5))
for (i in 1:n_regions) points(x1[test$regions[1:(n/3),i]],x2[test$regions[1:(n/3),i]],col=col[i],bg=bg[i],pch=21,cex=0.7)
points(x1[apply(test$regions,1,function(i) sum(i)==0)],
       x2[apply(test$regions,1,function(i) sum(i)==0)],
       col="black",bg=scales::alpha("black",0.25),pch=21,cex=0.7)
plot.new()
legend("center",legend=col.labels,col=c(col,1),pt.bg=c(bg,scales::alpha("black",0.25)),
       pt.cex=0.7,pch=21,ncol=1,cex=0.64,bty="n")
