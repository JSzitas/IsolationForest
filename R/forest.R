
isolationForest <- function( X,
                             n_trees = 100,
                             max_depth = NULL,
                             extension_level = 0,
                             Phi = 256,
                             seed = 1071,
                             vanilla = FALSE,
                             encode_categories = FALSE,
                             categorical_variables = NULL,
                             category_encodings_methods = NULL,
                             parallel = TRUE,
                             future_plan = "multiprocess" )
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
  if(parallel == TRUE){
    future::plan(future_plan)
    on.exit(future::plan("default"), add = TRUE)
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
                                   vanilla = vanilla,
                                   parallel = parallel,
                                   future_plan = future_plan )
  class(isolation_forest_object) <- "Isolation Forest"
  return( isolation_forest_object )
}
