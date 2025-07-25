#' custom_cart_split
#'
#' This code is adapted from the vignette in https://github.com/cran/rpart/blob/master/tests/usersplits.R
#' The only thing we change from that example is the splitting function (stemp), which is the focus of our approach
#' For more details on the other functions and rpart more generally, check the vignette and package documentation
#' @param y Outcome variable
#' @param x Explanatory variable
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' For example, we will consider that PTA holds in regions where beta_1 - beta_0 <= epsilon
cart_split <- function(y,x,epsilon)
{
  ### 1) init function
  itemp <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1, numy=1,
         summary= function(yval, dev, wt, ylevel, digits ) {
           paste("  mean=", format(signif(yval, digits)),
                 ", MSE=" , format(signif(dev/wt, digits)),
                 sep='')
         })
  }
  ### 2) 'evaluation' function
  etemp <-  function(y, wt, parms) {
    wmean <- sum(y*wt)/sum(wt)
    rss <- sum(wt*(y-wmean)^2)
    list(label= wmean, deviance=rss)
  }
  ### 3) splitting function
  stemp <- function(y, wt, x, parms, continuous) {
    # Center y
    n <- length(y)
    ym <-  sum(y*wt)/sum(wt)
    y <- y-ym

    if (continuous) {
      # continuous x variable
      temp <- cumsum(y*wt)[-n]

      left.wt  <- cumsum(wt)[-n]
      right.wt <- sum(wt) - left.wt
      lmean <- temp/left.wt
      rmean <- -temp/right.wt
      goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2) + n*var(y)*((abs(rmean+ym) < epsilon) | (abs(lmean+ym) < epsilon))
      list(goodness= goodness, direction=sign(lmean))
    }
    else {
      # Categorical X variable
      ux <- sort(unique(x))
      wtsum <- tapply(wt, x, sum)
      ysum  <- tapply(y*wt, x, sum)
      means <- ysum/wtsum

      # For anova splits, we can order the categories by their means
      #  then use the same code as for a non-categorical
      ord <- order(means)
      n <- length(ord)
      temp <- cumsum(ysum[ord])[-n]
      left.wt  <- cumsum(wtsum[ord])[-n]
      right.wt <- sum(wt) - left.wt
      lmean <- temp/left.wt
      rmean <- -temp/right.wt
      list(goodness= (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)  + n*var(y)*((abs(rmean+ym) < epsilon) | (abs(lmean+ym) < epsilon))  ,
           direction = ux[ord])
    }
  }
  ### 4) List to be passed on to the rpart function with our custom splitting criteria
  ulist <- list(eval = etemp, split = stemp, init = itemp)
  ### Return list for rpart function
  return(ulist)
}
#' placebo.cart
#'
#' Fit CART to predictions for main data from placebo regression
#' @param x Main data feature set
#' @param beta1 beta_1 prediction for main data from placebo regression
#' @param beta0 beta_0 prediction for main data from placebo regression
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
placebo.cart <- function(x,beta1,beta0,epsilon)
{
  ## Define output
  y <- c(beta1,-beta0)
  ## Adjust X format (e.g. one-hot encode levels of factors etc.)
  xm <- model.matrix(~.-1, data = x)
  ## data.frame with y,x (x has to be duplicated since we use both beta_1(x) and -beta_0(x) as outputs)
  xm <- data.frame(y,x = rbind(xm,xm))
  ## Fit CART tree
  out <- rpart::rpart(y~., data = xm,method = cart_split(y,x,epsilon) , cp = 0)
  return(out)
}
#' pta.nodes
#'
#' This function extracts the finest partitions of the data which CART flags as PTA regions
#' and stores which points belong in each region
#' @param placebo_cart CART tree, output from placebo.cart function
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
pta.nodes <- function(placebo_cart,epsilon)
{
  N <- length(placebo_cart$y)/2
  ## Determine which nodes have b1-b0 close enough to zero
  ind <- which(abs(placebo_cart$frame$yval)<epsilon)
  ind <- rownames(placebo_cart$frame)[ind]
  ## If there are no regions identified by CART, we simply return a nx1 matrix with all elements equal to FALSE
  if (length(ind)==0) return(matrix(FALSE,N,1))
  ## The path.rpart function gets the splitting rule for a specified node or set of nodes
  ## We will use it to determine the regions where PTA is likely to hold
  ## First, we get the path for every node (internal or terminal)
  path <- rpart::path.rpart(placebo_cart, node = as.numeric(rownames(placebo_cart$frame)),print.it = FALSE)
  ## The output of path.rpart is a character vector of the form c("root","split.rule.1","split.rule.2") etc.
  ## Now, we map these vectors to node names as defined by rpart
  ## The object final.node.in.path obtains the last node which appears in a certain path
  ### For example, node 1 is associated with "root". Suppose node 2 is defined by root -> x2 < 0.5.
  ### We want node 2 to be associated with "x2 < 0.5", so the path "root -> x2 < 0.5" becomes "1 -> 2".
  ### Similarly, suppose node 4 is defined by root -> x2 < 0.5 -> x4 < 0.5.
  ### We want node 4 to be associated with "x4 < 0.5", so the path "root -> x2 < 0.5 -> x4 < 0.5" becomes "1 -> 2 -> 4".
  final.node.in.path <- do.call("rbind",lapply(path,tail,1))
  final.node.in.path <- cbind(node=rownames(final.node.in.path),split=final.node.in.path)
  path.nodes <- path
  for (i in 1:length(path))
  {
    ## This step merges the split rules in a given path to the nodes associated with each rule
    temp <- merge(final.node.in.path,path[[i]],by.x=2,by.y=1,sort=FALSE)
    path.nodes[[i]] <- temp[,2]
  }
  ## Keeping only nodes that are flagged as PTA regions by CART
  nodes.pta <- path.nodes[ind]
  ## Now, we filter to keep only the finest partition possible. E.g. if {2}, {3,4} and {2,3,4} are all flagged,
  ### we want to keep only {2} and {3,4}
  ## To do that, we remove all nodes that are on the path to other nodes. E.g. if nodes 3,6,7 are flagged,
  ### node 3 is on the path to nodes 6 and 7, so it is removed
  if (length(ind)>1) # If there is only one node flagged, this step is innocuous
  {
    index <- NULL
    for (i in 1:(length(ind)-1))
    {
      for (j in (i+1):length(ind))
      {
        if (ind[i] %in% nodes.pta[[ind[j]]])
        {
          index <- append(index,i)
          break
        }
      }
    }
    ind <- ind[1:length(ind) %in% index == FALSE]
  }
  ## Now we get the full path for each point i \in {1,...,N}
  ## Element "where" in the rpart object gives:
  ### "the row number of frame corresponding to the leaf node that each observation falls into"
  ### rownames of the "frame" element of the rpart object gives the names of the nodes
  ## The object leaf.vec stores the name of the leaf node each i falls into
  leaf.vec <- rownames(placebo_cart$frame)[placebo_cart$where]
  ## Now, we make a list of size N which stores the path that reaches the leaf node of each i
  leaf.per.point <- vector("list",N)
  for (i in 1:N)
  {
    leaf.per.point[[i]] <- path.nodes[[leaf.vec[i]]]
  }
  ## Finally, we create a N x length(ind) matrix. Each column of this matrix stores a boolean vector
  ### which equals TRUE when i crosses a given PTA region. E.g. if the columns refer to regions x = 2 and x \in {3,4}
  ### the matrix has 2 columns, one which tracks points with x=2, and one for points with x \in {3,4}.
  regions <- matrix(FALSE,N,length(ind))
  for (i in 1:N)
  {
    for (j in 1:ncol(regions))
    {
      if (ind[j] %in% leaf.per.point[[i]]) regions[i,j] <- TRUE
    }
  }
  return(regions)
}
#' catt.per.region
#' Obtain the CATT predictions for the main dataset for the points in each PTA region
#' @param beta1 beta_1 prediction from main regression
#' @param beta0 beta_0 prediction from main regression
#' @param regions PTA regions, output from pta.nodes function
catt.per.region <- function(beta1,beta0,regions)
{
  tau <- beta1 - beta0
  out <- rep(NA,length(tau))
  for (i in 1:nrow(regions))
  {
    part <- which(regions[i,])
    out[i] <- mean(tau[regions[,part]])
  }
  return(out)
}
#' searchPTA
#'
#' Wrapper function to perform all operations
#' @param x Feature set from main data
#' @param beta1_placebo beta_1 prediction for main data from placebo regression
#' @param beta0_placebo beta_0 prediction for main data from placebo regression
#' @param beta1_main beta_1 prediction from main regression
#' @param beta0_main beta_0 prediction from main regression
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @return List with 1 - CART tree exploring PTA regions; 2 - matrix of PTA regions; 3 - CATT predictions for x_i given its PTA region
#' @export
searchPTA <- function(x,beta1_placebo,beta0_placebo,beta1_main,beta0_main,epsilon)
{
  placebo_cart <- placebo.cart(x,beta1_placebo,beta0_placebo,epsilon)
  regions <- pta.nodes(placebo_cart,epsilon)
  catt <- catt.per.region(beta1_main,beta0_main,regions)
  return(list(cart=placebo_cart,regions=regions,catt=catt))
}
# ## ----pta.catt-----------------------------------------------------------------
# ### Function to perform operations on all nodes
# pta.catt <- function(df,x.ind,main_posterior,placebo_posterior,eps)
# {
#   x <- subset(df, select=x.ind)
#   ## Sample one placebo draw for each main draw
#   placebo_draws <- sample(1:ncol(placebo_posterior$b1),ncol(main_posterior$b1))
#   ## Perform operation on every pair (main_draw,placebo_draw)
#   temp <- lapply(1:ncol(main_posterior$b1), function(i) node.wrapper(x,main_posterior,i,placebo_posterior,placebo_draws[i],eps))
#   ## Organize results
#   out <- list(cart=vector("list",ncol(main_posterior$b1)),
#               regions=vector("list",ncol(main_posterior$b1)),
#               catt=matrix(0,nrow(main_posterior$b1),ncol(main_posterior$b1)))
#   for (i in 1:ncol(main_posterior$b1))
#   {
#     out$cart[[i]] <- temp[[i]]$cart
#     out$regions[[i]] <- temp[[i]]$regions
#     out$catt[,i] <- temp[[i]]$catt
#   }
#   return(out)
# }
