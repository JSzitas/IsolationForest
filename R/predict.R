#' Path Length of a tree node (lm splitting)
#'
#' @description Calculate the path length of an observation through a given tree
#'
#' @param X The observations to use.
#' @param Tree A given tree to take the path through.
#' @param current_depth The current depth of the search
#' (how many nodes we have passed in total).
#' @param node_index The current node.
#' @details Calculates how deep in a tree an observation has to travel before it either
#'          * reaches a terminal node
#'          * we reach max tree depth
#'          using lm splitting.
#' @return Maximal depth + a factor calculated using **c_factor**.
#'
#'


path_length_lm <- function(X, Tree, current_depth = 0, node_index = 1)
{
  if (Tree[[1]][node_index,"Type"] == -1 ){
    return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  }

  ifelse(
    (X[, Tree[[2]][current_depth+1, 1]] - Tree[[2]][current_depth+1,2]) < 0,
    path_length_lm(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
    path_length_lm(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
}

#' Path Length of a tree node (lm splitting)
#'
#' @description Calculate the path length of an observation through a given tree
#'
#' @param X The observations to use.
#' @param Tree A given tree to take the path through.
#' @param current_depth The current depth of the search
#' (how many nodes we have passed in total).
#' @param node_index The current node.
#' @details Calculates how deep in a tree an observation has to travel before it either
#'          * reaches a terminal node
#'          * we reach max tree depth
#'          using lm splitting.
#' @return Maximal depth + a factor calculated using **c_factor**.
#'
#'


path_length_residual <- function(X, Tree, current_depth = 0, node_index = 1, residual_degree)
{
  if (Tree[[1]][node_index,"Type"] == -1 ){
    return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  }


  fitted <- rowSums( X[, Tree[[2]][current_depth+1, 2]] %*% t(Tree[[2]][current_depth+1,
                                                            (4:(4+residual_degree-1))]))

  ifelse(
    (abs(X[, Tree[[2]][current_depth+1, 1]] - fitted)) < Tree[[2]][current_depth+1, 3],
    path_length_residual( X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"],
                          residual_degree),
    path_length_residual( X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"],
                          residual_degree))
}



#' Path Length of a tree node (extended splitting)
#'
#' @description Calculate the path length of an observation through a given tree
#'
#' @param X The observations to use.
#' @param Tree A given tree to take the path through.
#' @param current_depth The current depth of the search
#' (how many nodes we have passed in total).
#' @param node_index The current node.
#' @details Calculates how deep in a tree an observation has to travel before it either
#'          * reaches a terminal node
#'          * we reach max tree depth
#'          using extended splitting.
#' @return Maximal depth + a factor calculated using **c_factor**.
#'
#'


path_length <- function(X, Tree, current_depth = 0, node_index = 1)
{
  if (Tree[[1]][node_index,"Type"] == -1 ){
    return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  }

  ifelse(
    ((X - Tree[[2]][current_depth+1,1:ncol(X)]) %*%
          Tree[[2]][current_depth+1,(ncol(X)+1):(2*ncol(X))]) < 0,
    path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
    path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
}

#' Path Length of a tree node (vanilla splitting)
#'
#' @description Calculate the path length of an observation through a given tree
#'
#' @param X The observations to use.
#' @param Tree A given tree to take the path through.
#' @param current_depth The current depth of the search
#' (how many nodes we have passed in total).
#' @param node_index The current node.
#' @details Calculates how deep in a tree an observation has to travel before it either
#'          * reaches a terminal node
#'          * we reach max tree depth
#'          using vanilla splitting.
#' @return Maximal depth + a factor calculated using **c_factor**.
#'
#'

path_length_vanilla <- function(X, Tree, current_depth = 0, node_index = 1)
{
  if (Tree[[1]][node_index,"Type"] == -1 ){
    return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  }
cat(node_index, current_depth)
  ifelse(
    ( X[,Tree[[2]][current_depth+1,1]] - Tree[[2]][current_depth+1,2] ) < 0,
    path_length_vanilla(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
    path_length_vanilla(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
}


#' Predict from an Isolation Forest
#'
#' @description Calculate the Anomaly Scores from a fitted Isolation Forest
#'
#' @param object A fitted object of class "Isolation Forest"
#' @param newdata Data to use for the prediction.
#' @param knn_smoothed Whether to use k-nearest neighbor (knn) smoothing
#' on the final predictions, defaults to **FALSE**.
#' @param knn_k Number of clusters to use for knn smoothing.
#' @param knn_method Method to use for knn smoothing - "average"/"median"
#' @param knn_distance Distance to use for knn smoothing "euclidean"/"manhattan".
#'
#' @return A vector of anomaly scores for **newdata** fitted from the trees in **object**.
#' @export
#'
#'
#'
#'

predict.isolationForest <- function( object,
                                     newdata,
                                     knn_smoothed = FALSE,
                                     knn_k = 5,
                                     knn_method = "average",
                                     knn_distance = "euclidean" )
{
  # check is object is of class Isolation Forest
  if(class(object) != "Isolation Forest"){
    stop("Model is not of class Isolation Forest!")
  }

  # check if ncol(newdata) == ncol(training_data)
  if(ncol(newdata) != object$n_variables ){
    stop(paste("Newdata has", ncol(newdata),
               "columns, but original training data had",
               object$n_variables, "columns" ))
  }
  if(sum(unlist(is.na(newdata))) != 0){
    newdata[is.na(newdata)] <- sample(c(-1e9,1e9),1)
  }
  # if(object$parallel == TRUE){
  #   future::plan( object$future_plan )
  #   on.exit(future::plan("default"), add = TRUE)
  # }

  # parallel tree prediction
  # if( object$vanilla ){
    # vanilla
    paths <- #future.apply::future_
    sapply(object$forest, function(i){
      path_length_vanilla(as.matrix(newdata), i)
    })
  # }
  # else if(object$lm){
  #   paths <- #future.apply::future_sapply
  #   sapply(object$forest, function(i){
  #     path_length_lm(as.matrix(newdata), i)
  #   })
  # }
  # else if(object$residual){
  #   paths <- #future.apply::future_sapply
  #     sapply(object$forest, function(i){
  #       path_length_residual(as.matrix(newdata), i,
  #                            residual_degree = object$residual_degree)
  #     })
  # }
  # else{
  #   # parallel tree prediction
  #   # standard, non-vanilla
  #   paths <- #future.apply::future_sapply
  #     sapply(object$forest, function(i){
  #     path_length(as.matrix(newdata), i)
  #   })
  # }
  res <- 2^(-rowMeans(paths)/c_factor(object$Phi))
  if( knn_smoothed ){
    res <- Rfast::knn( xnew = as.matrix(newdata),
                       y = res,
                       x = as.matrix(newdata),
                       dist.type = knn_distance,
                       type = "R",
                       k = knn_k,
                       method = knn_method )
  }
  return(res)
}


#' Calculates the 'path' factor
#'
#' @description Calculate the path factor of n
#'
#' @param n A positive integer.
#' @return A positive numeric of the path factor
#'

c_factor <- function(n) {
  if (n == 2){
    res <- 1
  }
  else if (n < 2){
    res <- 0
  }
  else {
    H <- log(n - 1) + 0.5772156649
    res <- 2 * H - (2*(n - 1)/n)
  }
  return(res)
}


