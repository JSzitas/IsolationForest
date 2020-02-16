#


split <- function( X,
                   ext_level,
                   vanilla = FALSE )
{
  if(vanilla){
    normal_vector <- rep(0, ncol(X))
    intercept_points <- rep(0, ncol(X))
    var_sample <- sample(1:ncol(X), 1)
    min_split <- min( unlist(X[,var_sample]))
    max_split <- max( unlist(X[,var_sample]))

    vanilla_comparison <- runif(1, min_split, max_split)

    vanilla_sample <- c( var_sample, vanilla_comparison )

    res <- unlist(X[,var_sample]) - vanilla_comparison

    return(list( filter = which( res < 0),
                 vanilla_fitted = vanilla_sample ))
  }
  else{
    # split performs splits on a variable
    # this find the lower and upper bounds for the values in columns of X
    mins <- unlist( lapply(1:ncol(X), function(i){
      min( unlist(X[,i]) )
    }))
    maxes <- unlist( lapply(1:ncol(X), function(i){
      max( unlist(X[,i]) )
    }))

    index <- sample(1:ncol(X), (ncol(X) - ext_level - 1), replace = FALSE)
    # Pick the indices for the normal vector elements
    normal_vector <- rnorm(ncol(X), mean = 0, sd = 1)
    normal_vector[index] <- 0
    # use indexes to pick the dimensions on which to use this
    intercept_points <- unlist(lapply(1:ncol(X),function(i){
      runif(1, min = mins[i], maxes[i])
    }))

    res <- (X - intercept_points) %*% normal_vector

    return(list( filter = which( res < 0 ),
                 normal = normal_vector,
                 intercept = intercept_points ))
  }
}


recurse <- function( index_data,
                     current_depth,
                     max_depth,
                     node_index = 0,
                     envir,
                     ext_level,
                     vanilla )
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
  res <- split( as.matrix( MIA_data),
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

  return( list(env$mat,
               env$normal_intercept_mat))
}




