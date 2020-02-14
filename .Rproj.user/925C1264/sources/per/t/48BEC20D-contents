split <- function( X,
                   ext_level,
                   vanilla = FALSE )
{
  if(vanilla){
    normal_vector <- rep(0, ncol(X))
    intercept_points <- rep(0, ncol(X))
    var_sample <- sample(1:ncol(X), 1)

    min_split <- min(unlist(X[,var_sample]))
    max_split <- max(unlist(X[,var_sample]))

    vanilla_comparison <- runif(1, min_split, max_split)

    vanilla_sample <- c( var_sample, vanilla_comparison )

    res <- unlist(X[,var_sample]) - vanilla_comparison
    return(list( #value = res,
      filter = which( res < 0 ),
      vanilla_fitted = vanilla_sample ))
  }
  else{
    # split performs splits on a variable
    # this find the lower and upper bounds for the values in columns of X
    mins <- unlist( lapply(1:ncol(X), function(i){
      min(unlist(X[,i]))
    }))
    maxes <- unlist( lapply(1:ncol(X), function(i){
      max(unlist(X[,i]))
    }))

    index <- sample(1:ncol(X), (ncol(X) - ext_level - 1), replace = FALSE)
    # Pick the indices for which the normal vector elements
    normal_vector <- rnorm(ncol(X), mean = 0, sd = 1)
    normal_vector[index] <- 0
    # use indexes to pick the dimensions on which to use this
    intercept_points <- unlist(lapply(1:ncol(X),function(i){
      runif(1, min = mins[i], maxes[i])
    }))

    res <- (X - intercept_points) %*% normal_vector

    return(list( #value = res,
      filter = which( res < 0 ),
      normal = normal_vector,
      intercept = intercept_points))
  }
}


# X = data, e = current depth, l = max depth, ni = node index
recurse <- function( index_data,
                     current_depth,
                     max_depth,
                     node_index = 0,
                     envir,
                     ext_level,
                     vanilla )
{
  ## don't sample columns with all duplicates
  duplicates <- sapply( envir$X[ index_data, , drop = F ],
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

  # perform splitting, using the current indexing of the data, and the correct
  # extension level
  res <- split( as.matrix( envir$X[index_data,] ),
                ext_level, vanilla )

  # modify matrix in place
  envir$mat[node_index, c("Left")] <- nodeLeft <- 2 * node_index
  envir$mat[node_index, c("Right")] <- nodeRight <- 2 * node_index + 1
  envir$mat[node_index, c( "Type")] <- c( 1)

  if(vanilla == TRUE){
    envir$normal_intercept_mat[ current_depth + 1,] <- c( res$vanilla_fitted )
  }
  else{
    envir$normal_intercept_mat[ current_depth + 1,] <- c( res$normal,
                                                          res$intercept )
  }

  # recurse to the left and to the right until termination is reached -
  # thus the function iteratively calls first its left nodes and then
  # its right nodes, until we are done.
  recurse( index_data[res$filter ], current_depth + 1,
           max_depth, nodeLeft,  envir, ext_level, vanilla )
  recurse( index_data[-res$filter], current_depth + 1,
           max_depth, nodeRight, envir, ext_level, vanilla )
}


iTree <- function( X,
                   max_tree_depth,
                   extension_level,
                   vanilla )
{
  env <- new.env() # pass everything in this environment to avoid copies
  env$mat <- matrix( 0,
                     nrow = max_nodes(max_tree_depth),
                     ncol = 5,
                     dimnames =
                       list(NULL, c( "TerminalID", "Type","Size","Left",
                                     "Right"))
  )
if(vanilla == TRUE){
  env$normal_intercept_mat <- matrix( 0,
                                      nrow = max_tree_depth,
                                      ncol = 2,
                                      dimnames =
                                        list( NULL,
                                              c( "Split_Variable_Vanilla",
                                                 "Split_Value_Vanilla")))
  }
else{
  env$normal_intercept_mat <- matrix( 0,
                                      nrow = max_tree_depth,
                                      ncol = 2*ncol(X),
                                      dimnames =
                                        list( NULL,
                                              c(paste("Normals ", 1:ncol(X),
                                                      col = ""),
                                                paste("Intercepts", 1:ncol(X),
                                                      col = "" ))))
  }
  env$X <- X

  recurse( index_data = 1:nrow(X),
           current_depth = 0,
           max_depth = max_tree_depth,
           node_index = 1,
           envir = env,
           ext_level = extension_level,
           vanilla )
  ## lots of allocated matrix space goes unused because max_nodes is an upper bound
  ## on how large the tree can be. This function subsets to only the rows needed
  ## for prediction
  # compress_matrix <- function(m) {
  #   m = cbind(seq.int(nrow(m)), m)[m[,"Type"] != 0,,drop=FALSE]
  #   m[,"Left"] = match(m[,"Left"], m[,1], nomatch = 0)
  #   m[,"Right"] = match(m[,"Right"], m[,1], nomatch = 0)
  #   m[m[,"Type"] == -1,"TerminalID"] = seq.int(sum(m[,"Type"] == -1))
  #   m[,-1,drop=FALSE]
  # }
  return( list(env$mat,
               env$normal_intercept_mat))
}


isolationForest <- function( X,
                             n_trees = 100,
                             max_depth = NULL,
                             extension_level = 0,
                             Phi = 256,
                             seed = 1071,
                             vanilla = FALSE )
{
  set.seed(seed)

  if(is.null(max_depth)){
    max_depth <- ceiling( log2( Phi ) )
  }

  forest = vector("list", n_trees)

  forest <- future.apply::future_lapply(1:n_trees, function(i){
    iTree(X[sample(1:nrow(X),Phi),], max_depth, extension_level, vanilla )
  }, future.seed = TRUE )


  isolation_forest_object <- list( forest = forest,
                                   Phi    = Phi,
                                   max_depth = max_depth,
                                   extension_level = extension_level,
                                   n_trees = n_trees,
                                   n_variables  = ncol(X),
                                   vanilla = vanilla )
  class(isolation_forest_object) <- "Isolation Forest"
  return( isolation_forest_object )
}

path_length <- function(X, Tree, e = 0, ni = 1)
{
  if (Tree[[1]][ni,"Type"] == -1 ){
    return(e + c_factor(Tree[[1]][ni,"Size"]))
  }

  ifelse(
    ((X - Tree[[2]][e+1,1:ncol(X)]) %*% Tree[[2]][e+1,(ncol(X)+1):(2*ncol(X))]) < 0,
    path_length(X, Tree, e + 1, Tree[[1]][ni,"Left"]),
    path_length(X, Tree, e + 1, Tree[[1]][ni,"Right"]))
}


path_length_vanilla <- function(X, Tree, e = 0, ni = 1)
{
  if (Tree[[1]][ni,"Type"] == -1 ){
    return(e + c_factor(Tree[[1]][ni,"Size"]))
  }

  ifelse(
    ( X[,Tree[[2]][e+1,1]] - Tree[[2]][e+1,2] ) < 0,
    path_length_vanilla(X, Tree, e + 1, Tree[[1]][ni,"Left"]),
    path_length_vanilla(X, Tree, e + 1, Tree[[1]][ni,"Right"]))
}

predict.isolationForest <- function(object, newdata, ...)
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
# parallel tree prediction
  if( object$vanilla ){
    # vanilla
    paths <- future.apply::future_sapply(object$forest, function(i){
      path_length_vanilla(as.matrix(newdata), i)
    })
  }
  else{
    # parallel tree prediction
    # standard, non-vanilla
    paths <- future.apply::future_sapply(object$forest, function(i){
      path_length(as.matrix(newdata), i)
    })
  }
  res <- 2^(-rowMeans(paths)/cn(object$Phi))
  return(res)
}


# calculates the 'path' factor for a given number
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




