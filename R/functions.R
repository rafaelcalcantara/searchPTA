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
#' fine_part
#'
#' Auxiliary function to obtain the finest partition for which delta1-delta0 ~ 0 in the CART tree
#' @param list List that contains all regions where PTA holds, from the most aggregated to most disaggregated ones
fine_part <- function(list)
{
  out <- sapply(1:length(list), function(i) {
    any(sapply((1:length(list))[-i], function(j) {
      all(list[[i]] %in% list[[j]])
    }))
  })
  return(!out)
}
#' points.per.region
#'
#' Auxiliary function to create a matrix with N rows and ncol = number of PTA regions
#' Each column is a boolean vector of size N equal to TRUE if a point falls into the region
#' associated to that column
#' @param listA List of finest partitions for which PTA holds
#' @param listB List which stores the path that reaches the leaf node of each i
points.per.region <- function(listA,listB)
{
  out <- matrix(FALSE,nrow=length(listB),ncol=length(listA))
  for (i in 1:nrow(out))
  {
    ind <- which(as.numeric(names(listA)) %in% listB[[i]])
    if (length(ind)>0) out[i,ind] <- TRUE
  }
  colnames(out) <- names(listA)
  rownames(out) <- names(listB)
  return(out)
}
#' custom_cart_split
#'
#' This code is adapted from the vignette in https://github.com/cran/rpart/blob/master/tests/usersplits.R
#' The only thing we change from that example is the splitting function (stemp), which is the focus of our approach
#' For more details on the other functions and rpart more generally, check the vignette and package documentation
#' @param y Outcome variable
#' @param x Explanatory variable
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' i.e. we will consider that PTA holds in regions where delta_1 - delta_0 <= epsilon
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
    # We perform the calculations per G=g
    w <- as.numeric(wt==1)
    n <- length(y)
    w1 <- w/sum(w) + (1-w)/sum(1-w)

    if (continuous) {
      # continuous x variable
      # Calculate sum-of-squares
      ## G=1
      # temp <- cumsum(y1*w)[-n]
      # left.wt  <- cumsum(w)[-n]
      # right.wt <- n - left.wt
      # lmean1 <- temp/left.wt
      # rmean1 <- -temp/right.wt
      # ## If weight was 0, we get Inf for mean; it should be 0, as up until that point in the cumsum,
      # ## the unit still had g==0 (conv. g==1). Set it to 0
      # lmean1[is.nan(lmean1)] <- 0
      # rmean1[is.nan(rmean1)] <- 0
      # goodness <- (left.wt*lmean1^2 + right.wt*rmean1^2)/(sum((y1*w)^2)/sum(w) + sum((y1*(1-w))^2)/sum(1-w))
      # ## G=0
      # temp <- cumsum(y0*w)[-n]
      # left.wt  <- cumsum(1-w)[-n]
      # right.wt <- n - left.wt
      # lmean0 <- temp/left.wt
      # rmean0 <- -temp/right.wt
      # ## If weight was 0, we get Inf for mean; it should be 0, as up until that point in the cumsum,
      # ## the unit still had g==0 (conv. g==1). Set it to 0
      # lmean0[is.nan(lmean0)] <- 0
      # rmean0[is.nan(rmean0)] <- 0
      # goodness <- goodness + (left.wt*lmean0^2 + right.wt*rmean0^2)/(sum((y1*w)^2)/sum(w) + sum((y1*(1-w))^2)/sum(1-w))
      ## Add term for difference in delta means
      temp1 <- cumsum(y*w)/cumsum(w)
      temp0 <- cumsum(y*(1-w))/cumsum(1-w)
      temp3 <- (sum(y*w)-cumsum(y*w))/(sum(w)-cumsum(w))
      temp2 <- (sum(y*(1-w))-cumsum(y*(1-w)))/(sum(1-w)-cumsum(1-w))
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
      list(goodness= goodness, direction=rep(-1,n-1))
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
#' @param delta1 delta_1 prediction for main data from placebo regression
#' @param delta0 delta_0 prediction for main data from placebo regression
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @param cp Complexity parameter from the rpart function
#' @param ... Additional arguments for rpart.control --- control CART tree growth
placebo.cart <- function(x,delta1,delta0,epsilon,wt,...)
{
  ## Define output
  y <- c(delta1,-delta0)
  ## Adjust X format (e.g. one-hot encode levels of factors etc.)
  xm <- model.matrix(~.-1, data = rbind(subset(x,g==1,select=which(names(x)!="g")),subset(x,g==0,select=which(names(x)!="g"))))
  ## data.frame with y,x
  xm <- data.frame(y,x = xm)
  ## Fit CART tree
  out <- rpart::rpart(y~., data = xm,method = cart_split(y,x,epsilon), weights = wt, ...)
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
  N <- length(placebo_cart$y)
  ## Determine which nodes have b1-b0 close enough to zero
  ind <- which(abs(placebo_cart$frame$yval)<epsilon)
  ind <- rownames(placebo_cart$frame)[ind]
  ## If there are no regions identified by CART, we simply return a nx1 matrix with all elements equal to FALSE
  if (length(ind)==0) return(matrix(FALSE,N,1))
  ## First, we get the path for every node (internal or terminal)
  path.nodes <- get_node_paths(placebo_cart)
  ## Keeping only nodes that are flagged as PTA regions by CART
  nodes.pta <- path.nodes[ind]
  ## Now, we filter to keep only the finest partition possible. E.g. if {2}, {3,4} and {2,3,4} are all flagged,
  ### we want to keep only {2} and {3,4}
  ## To do that, we remove all nodes that are on the path to other nodes. E.g. if nodes 3,6,7 are flagged,
  ### node 3 is on the path to nodes 6 and 7, so it is removed
  finest.partitions <- fine_part(nodes.pta)
  nodes.pta <- nodes.pta[finest.partitions]
  ## Now we get the full path for each point i \in {1,...,N}
  ## Element "where" in the rpart object gives:
  ### "the row number of frame corresponding to the leaf node that each observation falls into"
  ### rownames of the "frame" element of the rpart object gives the names of the nodes
  ## The object leaf.vec stores the name of the leaf node each i falls into
  leaf.vec <- rownames(placebo_cart$frame)[placebo_cart$where]
  ## Now, we make a list of size N which stores the path that reaches the leaf node of each i
  leaf.per.point <- vector("list",N)
  for (i in 1:N) leaf.per.point[[i]] <- path.nodes[[leaf.vec[i]]]
  ## Finally, we create a N x length(ind) matrix. Each column of this matrix stores a boolean vector
  ### which equals TRUE when i crosses a given PTA region. E.g. if the columns refer to regions x = 2 and x \in {3,4}
  ### the matrix has 2 columns, one which tracks points with x=2, and one for points with x \in {3,4}.
  regions <- points.per.region(nodes.pta,leaf.per.point)
  return(regions)
}
#' results.per.region
#' Output results
#' @param x feature set
#' @param delta1_aux delta_1 prediction from auxiliary regression
#' @param delta0_aux delta_0 prediction from auxiliary regression
#' @param delta1_main delta_1 prediction from main regression
#' @param delta0_main delta_0 prediction from main regression
#' @param regions PTA regions, output from pta.nodes function
#' @param CART rpart object, either NULL, if saveCART==FALSE, or the object, if saveCART==TRUE
#' @return list with: 1- Difference in betas per region; 2- CATT estimates per region; 3- matrix with regions; 4- lists with x per region and x with no flagged region; 6- CART tree
results.per.region <- function(x,delta1_aux,delta0_aux,delta1_main,delta0_main,regions,CART=NULL)
{
  g <- x$g
  dbetas <- apply(regions,2, function(i) mean(delta1_aux[g==1 & i])-mean(delta0_aux[g==0 & i]))
  catt <- apply(regions,2, function(i) mean(delta1_main[g==1 & i])-mean(delta0_main[g==0 & i]))
  x.in.reg <- lapply(1:ncol(regions), function(i) x[regions[,i],])
  x.not.in.reg <- x[rowSums(regions)==0,]
  if (is.null(CART)) out <- list(beta.diff=dbetas,catt=catt,regions=regions,x.per.regions=x.in.reg,x.not.in.regions=x.not.in.reg)
  else out <- list(beta.diff=dbetas,catt=catt,regions=regions,x.per.regions=x.in.reg,x.not.in.regions=x.not.in.reg,cart=CART)
  return(out)
}
#' searchPTA
#'
#' Wrapper function to perform all operations
#' @param x Feature set from main data
#' @param delta1_aux delta_1 prediction from auxiliary regression
#' @param delta0_aux delta_0 prediction from auxiliary regression
#' @param delta1_main delta_1 prediction from main regression
#' @param delta0_main delta_0 prediction from main regression
#' @param epsilon Parameter that dictates how close the trends between two groups needs to be for PTA to be considered valid
#' @param saveCART Whether or not to save the CART tree fit in the first step
#' @param ... Additional arguments for rpart.control --- control CART tree growth
#' @export
searchPTA <- function(x,delta1_aux,delta0_aux,delta1_main,delta0_main,epsilon,saveCART=TRUE,...)
{
  delta1_aux <- delta1_aux[x$g==1]
  delta0_aux <- delta0_aux[x$g==0]
  g1.rows <- as.numeric(rownames(subset(x,g==1)))
  g0.rows <- as.numeric(rownames(subset(x,g==0)))
  wt <- c(rep(1,length(delta1_aux)),rep(2,length(delta0_aux)))
  placebo_cart <- placebo.cart(x=x,delta1=delta1_aux,delta0=delta0_aux,epsilon=epsilon,wt=wt,...)
  regions <- pta.nodes(placebo_cart=placebo_cart,epsilon=epsilon)
  regions <- as.matrix(regions[order(c(g1.rows,g0.rows)),])
  if (saveCART)
  {
    return(results.per.region(x,delta1_aux,delta0_aux,delta1_main,delta0_main,regions,CART=placebo_cart))
  } else
  {
    return(results.per.region(x,delta1_aux,delta0_aux,delta1_main,delta0_main,regions,CART=NULL))
  }
}
