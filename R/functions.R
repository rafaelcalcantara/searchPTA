#' get_node_paths
#'
#' Auxiliary function to determine the paths to all nodes -- internal or leaf -- in CART
#' We make use of the fact that the nodes are labeled as follows:
#' parent node i splits into left node 2i  and right node 2i+1
#' @param rpart_obj Output from the rpart call
get_node_paths <- function(rpart_obj) {
  node_ids <- as.numeric(rownames(rpart_obj$frame))

  paths <- lapply(node_ids, function(n) {
    path <- c()
    while (n >= 1) {
      path <- c(n, path)
      if (n == 1) break
      n <- n %/% 2
    }
    path
  })

  names(paths) <- node_ids
  return(paths)
}
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

    if (continuous) {
      temp1 <- cumsum(y*w)/cumsum(w)
      temp0 <- cumsum(y*(1-w))/cumsum(1-w)
      temp3 <- (-cumsum(y*w))/(sum(w)-cumsum(w))
      temp2 <- (-cumsum(y*(1-w)))/(sum(1-w)-cumsum(1-w))
      temp0 <- ifelse(is.na(temp0),0,temp0)
      temp1 <- ifelse(is.na(temp1),0,temp1)
      temp2 <- ifelse(is.na(temp2),0,temp2)
      temp3 <- ifelse(is.na(temp3),0,temp3)
      lmean <- temp1+temp0
      rmean <- temp3+temp2
      ####
      ym <- sum(y*w1)
      goodness <- ((lmean-ym)^2 + (rmean-ym)^2)/(sum((y-ym)^2))
      ####
      delta.diff <- n*var(y)*((abs(rmean) < epsilon) | (abs(lmean) < epsilon))
      goodness <- goodness[-n] + delta.diff[-n]
      # goodness <- delta.diff[-n]
      # goodness <- rep(-1,length(delta.diff[-n]))
      goodness <- ifelse(is.na(goodness),0,goodness)
      goodness <- ifelse(!is.finite(goodness),0,goodness)
      ## Store results
      list(goodness=goodness, direction=sign(lmean[-n]))
    }
    else {
      # Categorical X variable
      ux <- sort(unique(x))
      wsum <- tapply(w1, x, sum)
      ysum  <- tapply(y*w1, x, sum)
      means <- ysum/wsum

      # For anova splits, we can order the categories by their means
      #  then use the same code as for a non-categorical
      ord <- order(means)
      n <- length(ord)
      temp1 <- cumsum(ysum[ord]*wsum[ord])/cumsum(wsum[ord])
      temp0 <- cumsum(ysum[ord]*(1-wsum[ord]))/cumsum(1-wsum[ord])
      temp3 <- (sum(ysum[ord]*wsum[ord])-cumsum(ysum[ord]*wsum[ord]))/(sum(wsum[ord])-cumsum(wsum[ord]))
      temp2 <- (sum(ysum[ord]*(1-wsum[ord]))-cumsum(ysum[ord]*(1-wsum[ord])))/(sum(1-wsum[ord])-cumsum(1-wsum[ord]))
      lmean <- temp1+temp0
      rmean <- temp3+temp2
      ####
      ym <- sum(y*w1)
      goodness <- ((lmean-ym)^2 + (rmean-ym)^2)/(sum((y-ym)^2))
      ####
      delta.diff <- n*var(y)*((abs(rmean) < epsilon) | (abs(lmean) < epsilon))
      goodness <- goodness[-n] + delta.diff[-n]
      goodness <- ifelse(is.na(goodness),0,goodness)
      ## Store results
      list(goodness=goodness, direction=ux[ord])
    }
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
  pta.nodes <- which(abs(placebo_cart$frame$yval)<epsilon)
  if (length(which(abs(placebo_cart$frame$yval)<epsilon))==0)
  {
    if (saveCART) return(list(beta.diff=rep(NA,N),catt=rep(NA,N),CART=placebo_cart))
    return(list(beta.diff=rep(NA,N),catt=rep(NA,N)))
  }
  pta.nodes <- rownames(placebo_cart$frame)[pta.nodes]
  ## Get the path for every node (internal or terminal)
  path.per.node <- get_node_paths(placebo_cart)
  ## Get the path for each point in sample (placebo_cart$where returns the leaf node each point falls into)
  path.per.point <- path.per.node[placebo_cart$where]
  ## Obtain n x k matrix, where each column is one of 1:k PTA regions, values are TRUE if point passes through a region, FALSE otherwise
  points.in.pta.nodes <- sapply(pta.nodes, function(i) sapply(path.per.point, function(j) i %in% j))
  ## Obtain CATT per PTA region
  catt.per.pta.node <- apply(points.in.pta.nodes, 2, function(i) mean(bta1[g==1 & i])-mean(beta0[g==0 & i]))
  ## Obtain delta gamma per PTA region
  dgamma.per.pta.node <- apply(points.in.pta.nodes, 2, function(i) mean(gamma1[g==1 & i])-mean(gamma0[g==0 & i]))
  ## Get deepest PTA node for each point (NA if never crosses a PTA node)
  max.pta.node <- sapply(path.per.point, function(i) ifelse(length(pta.nodes[pta.nodes %in% i])==0,NA,max(pta.nodes[pta.nodes %in% i])))
  ## Create output list
  if (saveCART) return(list(beta.diff=dgamma.per.pta.node[max.pta.node],catt=catt.per.pta.node[max.pta.node],CART=placebo_cart))
  return(list(beta.diff=dgamma.per.pta.node[max.pta.node],catt=catt.per.pta.node[max.pta.node]))
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
  g1.rows <- as.numeric(rownames(subset(x,g==1)))
  g0.rows <- as.numeric(rownames(subset(x,g==0)))
  wt <- c(rep(1,length(gamma1_temp)),rep(2,length(gamma0_temp)))
  placebo_cart <- placebo.cart(x=x,gamma1=gamma1_temp,gamma0=gamma0_temp,epsilon=epsilon,wt=wt,...)
  return(results(x,gamma1=gamma1,gamma0=gamma0,bta1=bta1,beta0=beta0,placebo_cart=placebo_cart,epsilon=epsilon,saveCART=saveCART))
}
