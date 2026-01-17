#' cart_split
#'
#' This code is adapted from the vignette in https://github.com/cran/rpart/blob/master/tests/usersplits.R
#' The only thing we change from that example is the splitting function (stemp), which is the focus of our approach
#' For more details on the other functions and rpart more generally, check the vignette and package documentation
#' @param y Outcome variable
#' @param x Explanatory variable
cart_split <- function(y,x)
{
  ### 1) init function
  itemp <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=parms, numresp=1, numy=1,
         summary= function(yval, dev, wt, ylevel, digits ) {
           paste("  mean=", format(signif(yval, digits)),
                 ", MSE=" , format(signif(dev/wt, digits)),
                 sep='')
         })
  }
  ### 2) 'evaluation' function
  etemp <- function(y, wt, parms) {
    wmean <- sum(y*wt)/sum(wt)
    rss <- sum(wt*(y-wmean)^2)
    w <- as.numeric(wt==1)
    w1 <- w/sum(w) + (1-w)/sum(1-w)
    yval <- sum(y*w1)
    list(label = yval, deviance = rss)
  }
  ### 3) splitting function
  stemp <- function(y, wt, x, parms, continuous) {
    epsilon <- parms
    # We perform the calculations per G=g
    w <- as.numeric(wt==1)
    n <- length(y)
    w1 <- w/sum(w) + (1-w)/sum(1-w)

    temp1 <- cumsum(y*w)/cumsum(w)
    temp0 <- cumsum(y*(1-w))/cumsum(1-w)
    temp3 <- (-cumsum(y*w))/(sum(w)-cumsum(w))
    temp2 <- (-cumsum(y*(1-w)))/(sum(1-w)-cumsum(1-w))
    lmean <- temp1+temp0
    lmean <- ifelse(is.na(lmean)|!is.finite(lmean),NA,lmean)
    rmean <- temp3+temp2
    rmean <- ifelse(is.na(rmean)|!is.finite(rmean),NA,rmean)
    ####
    ym <- sum(y*w1)
    ####
    delta.diff <- sapply(1:(n-1), function(i) max(abs(rmean[i]),abs(lmean[i])))
    delta.diff <- ifelse(is.na(delta.diff),Inf,delta.diff)
    if (is.na(ym))
    {
      goodness <- rep(0,n-1)
    } else if (any(delta.diff < abs(ym)))
    {
      goodness <- 1/delta.diff
    } else
    {
      # delta.diff <- sapply(1:(n-1), function(i) max(abs(rmean[i]),abs(lmean[i])))
      # delta.diff <- ifelse(is.na(delta.diff),Inf,delta.diff)
      # goodness <- 1/delta.diff
      goodness <- runif(n-1,0,1)/10
    }
    # goodness <- goodness[-n]
    ## Store results
    list(goodness=goodness, direction=sign(lmean[-n]))
  }
  ### 4) List to be passed on to the rpart function with our custom splitting criteria
  ulist <- list(eval = etemp, split = stemp, init = itemp)
  ### Return list for rpart function
  return(ulist)
}
#' placebo.cart
#'
#' Fit CART to placebo (pre-trend) predictions
#' @param x Main data feature set
#' @param gamma1 gamma_1 prediction
#' @param gamma0 gamma_0 prediction
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @param wt weight vector used for calculating the g=1 and g=0 means
#' @param ... Additional arguments for rpart fit (see rpart.control)
placebo.cart <- function(x,gamma1,gamma0,epsilon,wt,...)
{
  ## Define output
  y <- c(gamma1,-gamma0)
  ## Adjust X format
  ### We keep numeric variables and ordered factor variables and one-hot encode unordered factor variables
  ### This implies that all variables will be evaluated using the continuous split criteria. This is the rpart
  #### default for numeric and ordered factors, but not for unordered factors, for which there is a different routine
  #### based on ordering the factors by the mean of y per level of the factor
  temp <- rbind(subset(x,g==1,select=which(names(x)!="g")),subset(x,g==0,select=which(names(x)!="g")))
  unord.factors <- which(sapply(temp, function(i) is.factor(i) & !is.ordered(i)))
  not.unord.factors <- which(sapply(temp, function(i) !is.factor(i) | is.ordered(i)))
  if (length(unord.factors)>0)
  {
    ## If there are unordered categorical features
    dummy.mat <- lapply(subset(temp,select=unord.factors),contrasts,contrasts=FALSE)
    xm <- model.matrix(~.-1, data = subset(temp,select=unord.factors), contrasts.arg=dummy.mat)
    xm <- cbind(subset(temp,select=not.unord.factors),xm)
  } else
  {
    xm <- temp
  }
  xm <- data.frame(y,x=xm)
  ## Fit CART tree
  out <- rpart::rpart(y~., data = xm,method = cart_split(y,x), weights = wt, parms=epsilon, ...)
  return(out)
}
#' results
#'
#' Obtain CATT and delta gamma for each point based on the deepest PTA node they reach
#' @param x feature set
#' @param gamma1 gamma_1 prediction
#' @param gamma0 gamma_0 prediction
#' @param bta1 beta_1 + tau_1 + alpha_1 prediction
#' @param beta0 beta_0 prediction
#' @param placebo_cart rpart object
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @param saveCART boolean: should rpart object be stored in the output
results <- function(x,gamma1,gamma0,bta1,beta0,placebo_cart,epsilon,saveCART)
{
  g <- x$g
  N <- length(g)
  ## Determine which nodes have b1-b0 close enough to zero
  pta.nodes <- which(abs(placebo_cart$frame$yval) < epsilon)
  if (length(pta.nodes)==0)
  {
    if (saveCART) return(list(beta.diff=rep(NA,N),catt=rep(NA,N),CART=placebo_cart))
    return(list(beta.diff=rep(NA,N),catt=rep(NA,N)))
  }
  pta.nodes <- rownames(placebo_cart$frame)[pta.nodes]
  ## Get the path for every node (internal or terminal)
  # path.per.node <- get_node_paths(placebo_cart)
  path.per.node <- get_node_paths(placebo_cart)
  ## Get the path for each point in sample (placebo_cart$where returns the leaf node each point falls into)
  path.per.point <- path.per.node[placebo_cart$where]
  ## Obtain n x k matrix, where each column is one of 1:k PTA regions, values are TRUE if point passes through a region, FALSE otherwise
  # t0 <- Sys.time()
  # points.in.pta.nodes <- sapply(pta.nodes, function(i) sapply(path.per.point, function(j) i %in% j))
  # t1 <- Sys.time()
  points.in.pta.nodes <- checkPTA(as.numeric(pta.nodes),path.per.point)
  # t2 <- Sys.time()
  ## Obtain CATT per PTA region
  catt.per.pta.node <- apply(points.in.pta.nodes, 2, function(i) mean(bta1[g==1 & i])-mean(beta0[g==0 & i]))
  ## Obtain delta gamma per PTA region
  dgamma.per.pta.node <- apply(points.in.pta.nodes, 2, function(i) mean(gamma1[g==1 & i])-mean(gamma0[g==0 & i]))
  ## Get deepest PTA node for each point (NA if never crosses a PTA node)
  max.pta.node <- sapply(path.per.point, function(i) ifelse(length(pta.nodes[pta.nodes %in% i])==0,NA,max(pta.nodes[pta.nodes %in% i])))
  ## Create output list
  if (saveCART) return(list(beta.diff=dgamma.per.pta.node[max.pta.node],catt=catt.per.pta.node[max.pta.node],n_regions=length(unique(max.pta.node[!is.na(max.pta.node)])),CART=placebo_cart))
  return(list(beta.diff=dgamma.per.pta.node[max.pta.node],catt=catt.per.pta.node[max.pta.node],n_regions=length(unique(max.pta.node[!is.na(max.pta.node)]))))
}

#' searchPTA
#'
#' Wrapper function to perform all operations
#' @param x Feature set from main data
#' @param gamma1 gamma_1 prediction from auxiliary regression
#' @param gamma0 gamma_0 prediction from auxiliary regression
#' @param bta1 beta_1 + tau_1 + alpha_1 prediction from main regression
#' @param beta0 beta_0 prediction from main regression
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @param saveCART boolean: should rpart object be stored in the output
#' @param ... Additional arguments for rpart fit (see rpart.control)
#' @export
searchPTA <- function(x,gamma1,gamma0,bta1,beta0,epsilon,saveCART=TRUE,...)
{
  gamma1_temp <- gamma1[x$g==1]
  gamma0_temp <- gamma0[x$g==0]
  bta1_temp <- bta1[x$g==1]
  beta0_temp <- beta0[x$g==0]
  wt <- c(rep(1,length(gamma1_temp)),rep(2,length(gamma0_temp)))
  placebo_cart <- placebo.cart(x=x,gamma1=gamma1_temp,gamma0=gamma0_temp,epsilon=epsilon,wt=wt,...)
  out <- results(x,gamma1=gamma1,gamma0=gamma0,bta1=bta1,beta0=beta0,placebo_cart=placebo_cart,epsilon=epsilon,saveCART=saveCART)
  return(out)
}
