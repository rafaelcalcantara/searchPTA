gam1 <- function(xvals){

  value <-  (xvals[,1]-0.5)^2 - xvals[,2]

  return(value)

}

gam0 <- function(xvals){

  value <-  xvals[,1] + xvals[,2]

  return(value)

}


f1 <- function(nsamp){

  x1 <- rbeta(nsamp,3,1)
  x2 <- rbeta(nsamp,1,3)

  x <- cbind(x1,x2)
  return(x)

}

f0 <- function(nsamp){

  x1 <- rbeta(nsamp,1,1)
  x2 <- rbeta(nsamp,5,2)

  x <- cbind(x1,x2)
  return(x)

}

# region <- function(xvals,parms){
#
#   lower1 <- parms[1]
#   upper1 <- parms[1] + exp(parms[2])
#   lower2 <- parms[3]
#   upper2 <- parms[3] + exp(parms[4])
#
#   indicies <- xvals[,1] > lower1 & xvals[,1] < upper1 & xvals[,2] > lower2 & xvals[,2] < upper2
#
#   return(indicies)
#
# }




n <- 50000

x1 <- f1(n)
x0 <- f0(n)

X <- data.frame(rbind(x1,x0),g = c(rep(1,n),rep(0,n)))

raf <- searchPTA::searchPTA(X,gam1(rbind(x1,x0)),gam0(rbind(x1,x0)),gam1(rbind(x1,x0)),gam0(rbind(x1,x0)),epsilon = 5)

# biascomp <- function(params){
#
#   idx1 <- region(x1,params)
#   idx0 <- region(x0,params)
#
# biasval <-  abs(mean(gam1(x1[idx1,])) - mean(gam0(x0[idx0,])))
#
# return(biasval)
#
# }
#
#
# bias <- biascomp(c(0,1,0,1))
#
# print(bias)
#
#
# foropt <- function(pars){
#
#   temp <- biascomp(pars) + 0.1*exp(pars[2]) + 0.1*exp(pars[3])
#   return(temp)
#
# }
#
# optrect <- optim(c(0.3,0,0.3,0),foropt)$par
#
#
# print(biascomp(optrect))
#
#
# temp <- cbind(optrect[1],optrect[1] + exp(optrect[2]), optrect[3], optrect[3] + exp(optrect[4]))
#
# temp <- pmin(temp,c(1,1,1,1))
#
# print(temp)
#


biascalc <- function(idx1,idx0){



  biasval <-  abs(mean(gam1(x1[idx1,])) - mean(gam0(x0[idx0,])))

  return(biasval)

}

x0 <- data.frame(x0)
names(x0) <- c("x.x1","x.x2")

x1 <- data.frame(x1)
names(x1) <- c("x.x1","x.x2")

preds1 <- predict(raf$CART,newdata = x1)
preds0 <- predict(raf$CART,newdata = x0)


for (j in 1:length(unique(preds1))){
  idx1 <- preds1==sort(unique(preds1)[j])
  idx0 <- preds0==sort(unique(preds0)[j])

  print(biascalc(idx1,idx0))}



