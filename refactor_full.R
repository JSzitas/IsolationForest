# double split tree building

parse_degree <- function( function_degree, passed ){
  X_smpl <- passed
  degree_vec <- paste( "X_smpl^", 1:function_degree, sep = "")
  result <- Reduce("cbind",lapply(degree_vec, FUN = function(i){
    eval(parse(text = i))
  }))
  result <- as.matrix(result, ncol = function_degree)
  return(result)
}

# residual_split <- function( X,
#                             residual_degree,
#                             lambda = 1)
# {
#     reg_cols <- sample( 1:ncol(X), 2)
#     Y_smpl <- X[,reg_cols[1]]
#     X_smpl <- X[,reg_cols[2]]
#
#     res <- parse_degree( function_degree = residual_degree, passed = X_smpl )
#
#     coef <- c( (solve(t(res) %*% res) + diag(lambda,
#                                            nrow = ncol(res),
#                                            ncol = ncol(res)) ) %*% t(res) %*% Y_smpl)
#
#     res <- abs( suppressWarnings( Y_smpl - coef * X_smpl))
#     res_smpl <- runif(1, min( res), max(res))
#
#     res <- res - res_smpl
#
#
#     # different polynomials, kernel? (ie gaussian/poly kernel)
#
#
#     return(list( filter = which( res < 0),
#                  fit_values =  #Y_var =
#                    c( reg_cols[1], # X_var =
#                       reg_cols[2], # comparison =
#                       res_smpl, # coef =
#                       coef )))
# }

residual_split <- function( X,
                            residual_degree,
                            lambda = 1)
{
  reg_cols <- sample( 1:ncol(X), residual_degree + 1)
  Y_smpl <- X[,reg_cols[1]]
  X_smpl <- X[,reg_cols[-1]]

  X_smpl <- as.matrix( cbind(rep(1, nrow(X)), X_smpl))
  # res <- parse_degree( function_degree = residual_degree, passed = X_smpl )
  # MASS::ginv
  coef <- c( (MASS::ginv(t(X_smpl) %*% X_smpl) + diag(rnorm(1),
                                                      nrow = ncol(X_smpl),
                                                      ncol = ncol(X_smpl)) ) %*% t(X_smpl) %*% Y_smpl)

  res <- abs( suppressWarnings(
    Reduce("+", lapply( 1:ncol(X_smpl), FUN = function(i){
      Y_smpl - coef * X_smpl[,i]}))))
  res_smpl <- runif(1, min( res), max(res))

  res <- res - res_smpl
  # different polynomials, kernel? (ie gaussian/poly kernel)

  return(list( filter = which( res < 0),
               fit_values =  #Y_var =
                 c( reg_cols, # sampled
                    res_smpl, # coef =
                    coef )))
}



vanilla_split <- function( X ){

  var_sample <- sample(1:ncol(X), 1)
  min_split <- min( unlist(X[,var_sample]))
  max_split <- max( unlist(X[,var_sample]))

  vanilla_comparison <- runif(1, min_split, max_split)

  res <- unlist(X[,var_sample]) - vanilla_comparison

  return(list( filter = which( res < 0),
               fit_values = c( var_sample,
                               vanilla_comparison )
  ))
}


extended_split <- function( X,
                            extension_level ){
  # split performs splits on a variable
  # this find the lower and upper bounds for the values in columns of X
  mins <- unlist( lapply(1:ncol(X), function(i){
    min( unlist(X[,i]) )
  }))
  maxes <- unlist( lapply(1:ncol(X), function(i){
    max( unlist(X[,i]) )
  }))

  index <- sample(1:ncol(X), (ncol(X) - extension_level - 1), replace = FALSE)
  # Pick the indices for the normal vector elements
  normal_vector <- rnorm(ncol(X), mean = 0, sd = 1)
  normal_vector[index] <- 0
  # use indexes to pick the dimensions on which to use this
  intercept_points <- unlist(lapply(1:ncol(X),function(i){
    runif(1, min = mins[i], maxes[i])
  }))

  res <- as.matrix(X - intercept_points) %*% normal_vector

  return(list( filter = which( res < 0 ),
               fit_values = c( normal_vector,
                               intercept_points)
  ))
}

# double_split <- function( X ){
#   var_sample <- sample(1:ncol(X), 1)
#   min_split <- min( unlist(X[,var_sample]))
#   max_split <- max( unlist(X[,var_sample]))
#
#   comparison <- runif(2, min_split, max_split)
#
#   res <- union( which((unlist(X[,var_sample]) < min(comparison))),
#                 which((unlist(X[,var_sample]) > max(comparison))) )
#
#   return(list( filter = res,
#                fit_values = c( var_sample,
#                                min(comparison),
#                                max(comparison))
#                ))
#
# }

# kernel_split <- function( X ){
#
#   kernel <-  exp(-1/2*(1^2) *(
#     as.matrix( dist( sapply( X, FUN = function(i){
#       (i - mean(i))/sd(i) }
#       )))^2))
#
#   var_sample <- sample(1:ncol(kernel), 1)
#
#   norm()
#
#   return(list( filter = res,
#                fit_values = c( var_sample,
#                                min(comparison),
#                                max(comparison))
#   ))
# }


# difference_split <- function( X ){
#
#   var_sample <- sample(1:ncol(X), 1)
#   vec <- unlist( X[, var_sample] )
#
#   vec <- data.frame( vec, id = 1:length(vec))
#   vec <- vec[ do.call( order, vec), ]
#   return(vec)
#
#   vec_differences <- vec[2:nrow(vec)] - vec[1:(nrow(vec) - 1), 1]
#   res <- runif(1, min = min( vec_differences ), max = max( vec_differences ) )
#
#   res <- vec_differences < res
#
#
#
#   return( list( filter =  ))
#
# }



# quantgrid_split <- function( X ){
#
#
#   quantile_grids
#
#
# }


variance_split <- function( X ){

  var_sample <- sample(1:ncol(X), 1)
  med <- median(unlist(X[, var_sample]))

  res <- abs( med - unlist(X[, var_sample]) )
  res_comparison <- runif(1, min = min(res), max = max(res))

  return( list( filter = which( res < res_comparison),
                fit_values = c( var_sample,
                                med,
                                res_comparison ))
  )

}


split <- function( X,
                   type = c( "vanilla",
                             "extended",
                             "residual",
                             "variance"),
                   # "double_split"),
                   # "kernel"),
                   # "difference_split"),
                   residual_degree = 1,
                   lambda = 1,
                   extension_level = 1 )
{
  split_type <- list()
  split_type[["vanilla"]] <- quote( vanilla_split(X) )
  split_type[["extended"]] <- quote( extended_split(X, extension_level))
  split_type[["residual"]] <- quote( residual_split(X, residual_degree, lambda))
  split_type[["variance"]] <- quote( variance_split(X))
  # split_type[["kernel_split"]] <- quote( kernel_split(X))

  result <- eval(split_type[[type]])

  return(result)
}


recurse <- function( index_data,
                     current_depth,
                     max_depth,
                     node_index = 0,
                     envir,
                     type_split,
                     residual_degree,
                     lambda,
                     extension_level )
{
  ## don't sample columns with all duplicates
  duplicates <- sapply( envir$X[ index_data, , drop = FALSE ],
                        function(x){
                          all(duplicated(x)[-1L])
                        })
  # Termination - either we are over the max depth limit, or we
  # we have come to the a split where there is only 1 observation
  # or all of the data in the split is the same - ie all same category,
  # or same value.
  if (current_depth >= max_depth || length(index_data) <= 1 || all(duplicates) ){
    envir$mat[ node_index,
               c("Type", "Size")] <- c(-1, length(index_data))
    # the empty return is here so that the function actually has a return somewhere.
    return()
  }

  # random MIA recoding >>
  # this randomly chooses the split to which the NA values will go
  if(sum(unlist(is.na(envir$X[index_data,]))) != 0){
    MIA_direction <- sample(c(-1e9,1e9),1)
    MIA_data <- envir$X[index_data,]
    MIA_data[is.na(envir$X[index_data,])] <- MIA_direction
  }
  else{
    MIA_data <- envir$X[index_data,]
  }

  # perform splitting, using the current indexing of the data, and the correct
  # extension level
  res <- split( MIA_data,
                type = type_split,
                residual_degree,
                lambda,
                extension_level )
  # modify matrix in place
  envir$mat[node_index, c("Left")] <- nodeLeft <- 2 * node_index
  envir$mat[node_index, c("Right")] <- nodeRight <- 2 * node_index + 1
  envir$mat[node_index, c( "Type")] <- c( 1)

  envir$pred_mat[ current_depth + 1,] <- c( res$fit_values )

  # recurse to the left and to the right until termination is reached -
  # thus the function iteratively calls first its left nodes and then
  # its right nodes, until we are done.
  recurse( index_data[res$filter ], current_depth + 1,
           max_depth, nodeLeft,  envir, type_split,
           residual_degree,
           lambda,
           extension_level )
  recurse( index_data[-res$filter], current_depth + 1,
           max_depth, nodeRight, envir, type_split,
           residual_degree,
           lambda,
           extension_level )
}

iTree <- function( X,
                   max_tree_depth,
                   split_type,
                   residual_degree,
                   lambda,
                   extension_level )
{
  env <- new.env() # pass everything in this environment to avoid copies
  env$mat <- matrix( 0,
                     nrow = max_nodes(max_tree_depth),
                     ncol = 5,
                     dimnames =
                       list(NULL, c( "TerminalID", "Type","Size","Left",
                                     "Right"))
  )

  res_table <- list()
  res_table[["vanilla"]] <- 2
  res_table[["extended"]] <- 2*ncol(X)
  res_table[["residual"]] <- 3 + residual_degree + 1
  res_table[["variance"]] <- 3
  # res_table[["double_split"]] <- 3

  ncol_pred_mat <- unlist( res_table[[split_type]] )

  env$pred_mat <- matrix( 0,
                          nrow = max_tree_depth,
                          ncol = ncol_pred_mat )

  env$X <- X

  recurse( index_data = 1:nrow(X),
           current_depth = 0,
           max_depth = max_tree_depth,
           node_index = 1,
           envir = env,
           type_split = split_type,
           residual_degree,
           lambda,
           extension_level)

  return( list(env$mat,
               env$pred_mat))
}


max_nodes <- function( max_tree_depth )
{
  return( 2*( 2^max_tree_depth ) - 1)
}


isolationForest <- function( X,
                             n_trees = 100,
                             max_depth = NULL,
                             Phi = 256,
                             seed = 1071,
                             type = c( "vanilla",
                                       "extended",
                                       "residual",
                                       "variance" ),
                             # "double_split",
                             # "random_sphere" ),
                             residual_degree = 1,
                             lambda = 1,
                             extension_level = 1,
                             encode_categories = FALSE,
                             categorical_variables = NULL,
                             category_encodings_methods = NULL,
                             parallel = TRUE)#,
# future_plan = "multiprocess" )
{
  set.seed(seed)

  if(is.null(max_depth)){
    max_depth <- ceiling( log2( Phi ) )
  }
  if(encode_categories){
    X <- data.frame(categoryEncodings::encode_categories( X,
                                                          fact = categorical_variables,
                                                          method = category_encodings_methods ))
  }
  # if(parallel == TRUE){
  #   future::plan(future_plan)
  #   on.exit(future::plan("default"), add = TRUE)
  # }

  forest = vector("list", n_trees)

  # forest <- future.apply::future_lapply(1:n_trees, function(i){
  #   iTree(X[sample(1:nrow(X),Phi),], max_depth, extension_level, vanilla, lm )
  # }, future.seed = TRUE )

  # if(cpp){
  #   for(i in 1:n_trees){
  #     forest[[i]] <- iTree_cpp( as.matrix( X[sample(1:nrow(X),Phi),] ),
  #                                           ext = extension_level,
  #                                           type = type,
  #                                           max_depth = max_depth)
  #   }
  #   forest <- lapply(forest, FUN = function(i){ colnames(i[[1]]) <- c( "TerminalID",
  #                                                                      "Type",
  #                                                                      "Size",
  #                                                                      "Left",
  #                                                                      "Right")
  #   return( i ) })
  #   if(type == "vanilla"){
  #     forest <- lapply(forest, FUN = function(i){ i[[2]][,1] <- i[[2]][,1] + 1
  #     return( i )
  #     })
  #   }
  # }
  # else{
  if(parallel){
    forest <-  replicate(n = n_trees, .f = function(i){iTree(
      X[sample(1:nrow(X),Phi),],
      max_depth,
      split_type = type,
      residual_degree,
      lambda,
      extension_level)}, parallel = TRUE)
  }
  else{
    for(i in 1:n_trees){
      forest[[i]] <- iTree(
        X[sample(1:nrow(X),Phi),],
        max_depth,
        split_type = type,
        residual_degree,
        lambda,
        extension_level)
    }
  }



  isolation_forest_object <- list( forest = forest,
                                   Phi    = Phi,
                                   max_depth = max_depth,
                                   extension_level = extension_level,
                                   n_trees = n_trees,
                                   n_variables  = ncol(X),
                                   type = type,
                                   residual_degree = residual_degree,
                                   lambda = lambda,
                                   extension_level = extension_level,
                                   parallel = parallel )
  # future_plan = future_plan )
  class(isolation_forest_object) <- "Isolation Forest"
  return( isolation_forest_object )
}




# parallel replicate
replicate <- function( n, .f, parallel = FALSE )
{
  if(parallel == FALSE)
  {
    result <- list()
    while(n > 1)
    {
      result[[n]] <- do.call(.f, args = list() )
      n <- n - 1
    }
  }
  else
  {
    ncore <- parallel::detectCores()
    clust <- parallel::makeCluster(ncore)
    parallel::clusterExport( cl = clust,
                             varlist = ls(envir = .GlobalEnv),
                             envir = .GlobalEnv  )
    result <- parallel::parLapply( clust, X = 1:n, fun = .f )
    parallel::stopCluster(clust)
  }
  return(result)
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



comparator_vanilla <- function( Tree, X, current_depth ){
  return((X[, Tree[[2]][current_depth+1, 1]] - Tree[[2]][current_depth+1,2]) < 0)
}

comparator_extended <- function(Tree, X, current_depth ){
  return( ((X - Tree[[2]][current_depth+1,1:ncol(X)]) %*%
             Tree[[2]][current_depth+1,(ncol(X)+1):(2*ncol(X))]) < 0)
}

comparator_residual <- function( Tree, X, current_depth ){

  fitted <- rowSums( X[, Tree[[2]][current_depth+1, 2]] %*%
                       t(Tree[[2]][ current_depth+1,
                                    (4:ncol(Tree[[2]]))]))

  return(
    (abs(X[, Tree[[2]][current_depth+1, 1]] - fitted)) < Tree[[2]][current_depth+1, 3])

}


comparator_variance <- function( Tree, X, current_depth ){
  return(
    ( abs(Tree[[2]][current_depth+1, 2] - X[, Tree[[2]][current_depth+1, 1]] )) <
      Tree[[2]][current_depth+1, 3]
  )
}


# comparator_double <- function( Tree, X, current_depth ){
#   return(
#     union( which( unlist(X[, Tree[[2]][current_depth+1, 1]] ) <
#                   Tree[[2]][current_depth+1, 2]),
#            which( unlist(X[, Tree[[2]][current_depth+1, 1]] ) >
#                   Tree[[2]][current_depth+1, 2]))
#   )
# }
#
#
# comparator_rbf_kernel <- function( Tree, X, current_depth ){
#   return(
#
#   )
# }



tree_type_selector <- function( type ){
  split_type <- list()
  split_type[["vanilla"]] <- comparator_vanilla
  split_type[["extended"]] <- comparator_extended
  split_type[["residual"]] <- comparator_residual
  split_type[["variance"]] <- comparator_variance
  # split_type[["rbf_kernel"]] <- comparator_rbf_kernel
  return(split_type[[type]])
}


path_length <- function(X, Tree, current_depth = 0, node_index = 1, path_len_fun ){

  if (Tree[[1]][node_index,"Type"] == -1 ){
    return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  }
  #   if (current_depth >7 ){ #|| current_depth + 1 >= ncol(Tree[[2]])  ){
  #     return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
  #   }
  # cat(current_depth, node_index)
  ifelse(
    do.call( path_len_fun, list( Tree, X, current_depth )),
    path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"], path_len_fun),
    path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"], path_len_fun))

}



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
  paths <- #future.apply::future_sapply
    sapply(object$forest, function(i){
      path_length( as.matrix(newdata), i,
                   path_len_fun = tree_type_selector( object$type))
    })

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



#'
#' path_length_lm <- function(X, Tree, current_depth = 0, node_index = 1)
#' {
#'   if (Tree[[1]][node_index,"Type"] == -1 ){
#'     return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
#'   }
#'
#'   ifelse(
#'     (X[, Tree[[2]][current_depth+1, 1]] - Tree[[2]][current_depth+1,2]) < 0,
#'     path_length_lm(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
#'     path_length_lm(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
#' }
#'
#' #' Path Length of a tree node (lm splitting)
#' #'
#' #' @description Calculate the path length of an observation through a given tree
#' #'
#' #' @param X The observations to use.
#' #' @param Tree A given tree to take the path through.
#' #' @param current_depth The current depth of the search
#' #' (how many nodes we have passed in total).
#' #' @param node_index The current node.
#' #' @details Calculates how deep in a tree an observation has to travel before it either
#' #'          * reaches a terminal node
#' #'          * we reach max tree depth
#' #'          using lm splitting.
#' #' @return Maximal depth + a factor calculated using **c_factor**.
#' #'
#' #'
#'
#'
#' path_length_residual <- function(X, Tree, current_depth = 0, node_index = 1, residual_degree)
#' {
#'   if (Tree[[1]][node_index,"Type"] == -1 ){
#'     return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
#'   }
#'
#'
#'   fitted <- rowSums( X[, Tree[[2]][current_depth+1, 2]] %*% t(Tree[[2]][current_depth+1,
#'                                                                         (4:(4+residual_degree-1))]))
#'
#'   ifelse(
#'     (abs(X[, Tree[[2]][current_depth+1, 1]] - fitted)) < Tree[[2]][current_depth+1, 3],
#'     path_length_residual( X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"],
#'                           residual_degree),
#'     path_length_residual( X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"],
#'                           residual_degree))
#' }
#'
#'
#'
#' #' Path Length of a tree node (extended splitting)
#' #'
#' #' @description Calculate the path length of an observation through a given tree
#' #'
#' #' @param X The observations to use.
#' #' @param Tree A given tree to take the path through.
#' #' @param current_depth The current depth of the search
#' #' (how many nodes we have passed in total).
#' #' @param node_index The current node.
#' #' @details Calculates how deep in a tree an observation has to travel before it either
#' #'          * reaches a terminal node
#' #'          * we reach max tree depth
#' #'          using extended splitting.
#' #' @return Maximal depth + a factor calculated using **c_factor**.
#' #'
#' #'
#'
#'
#' path_length <- function(X, Tree, current_depth = 0, node_index = 1)
#' {
#'   if (Tree[[1]][node_index,"Type"] == -1 ){
#'     return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
#'   }
#'
#'   ifelse(
#'     ((X - Tree[[2]][current_depth+1,1:ncol(X)]) %*%
#'        Tree[[2]][current_depth+1,(ncol(X)+1):(2*ncol(X))]) < 0,
#'     path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
#'     path_length(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
#' }
#'
#' #' Path Length of a tree node (vanilla splitting)
#' #'
#' #' @description Calculate the path length of an observation through a given tree
#' #'
#' #' @param X The observations to use.
#' #' @param Tree A given tree to take the path through.
#' #' @param current_depth The current depth of the search
#' #' (how many nodes we have passed in total).
#' #' @param node_index The current node.
#' #' @details Calculates how deep in a tree an observation has to travel before it either
#' #'          * reaches a terminal node
#' #'          * we reach max tree depth
#' #'          using vanilla splitting.
#' #' @return Maximal depth + a factor calculated using **c_factor**.
#' #'
#' #'
#'
#' path_length_vanilla <- function(X, Tree, current_depth = 0, node_index = 1)
#' {
#'   if (Tree[[1]][node_index,"Type"] == -1 ){
#'     return(current_depth + c_factor(Tree[[1]][node_index,"Size"]))
#'   }
#'
#'   ifelse(
#'     ( X[,Tree[[2]][current_depth+1,1]] - Tree[[2]][current_depth+1,2] ) < 0,
#'     path_length_vanilla(X, Tree, current_depth + 1, Tree[[1]][node_index,"Left"]),
#'     path_length_vanilla(X, Tree, current_depth + 1, Tree[[1]][node_index,"Right"]))
#' }
#'

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






