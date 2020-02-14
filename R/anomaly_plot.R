
anomaly_plot <- function( x,
                          y,
                          scores = NULL,
                          data = NULL,
                          forest = NULL,
                          x.min= NULL,
                          x.max = NULL,
                          y.min = NULL,
                          y.max = NULL,
                          point_plot = TRUE,
                          contamination = 0.05 )
{
  if(point_plot){
    if(is.null(data)){
      data <- cbind(x,y)
    }
    if(is.null(scores)){
      scores <- predict.isolationForest(forest, data)
    }

      Anomalies <- as.factor(scores > quantile(scores,(1-contamination)))
      levels(Anomalies) <- c("Normal","Anomaly")


      data_fm <- data.frame( cbind( data,
                                    Anomalies))

      ggplot2::ggplot(data_fm, ggplot2::aes( x, y, colour = Anomalies, shape = Anomalies)) +
        ggplot2::geom_point(size = 1.9) +
        ggplot2::scale_colour_manual(name = "Scoring", values = c("#2554C7","#E42217")) +
        ggplot2::scale_shape_manual(name = "Scoring", values = c(15,17)) +
        ggplot2::xlab("X") +
        ggplot2::ylab("Y") +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text( angle = 0, vjust = 0.5, size = 12),
          axis.title.x = ggplot2::element_text( size = 12) )



  }
  else{
  # contour
      if(is.character(x)){
        x_index <- grep( pattern = x, x = colnames(data) )
        x <- data$x
      }
      else{
        x_index <- x
        x <- data[, x]
      }
      if(is.character(y)){
        y_index <- grep( pattern = y, x = colnames(data) )
        y <- data$y
      }
      else{
        y_index <- y
        y <- data[, y]
      }

      if(is.null(x.min)){
        x.min <- min(unlist(x))
      }
      if(is.null(y.min)){
        y.min <- min(unlist(y))
      }
      if(is.null(x.max)){
        x.max <- max(unlist(x))
      }
      if(is.null(y.max)){
        y.max <- max(unlist(y))
      }

      # get column means for later
      means <- t(colMeans(data[,-c( x_index, y_index )]))
      var_names <- colnames( data )

      # define the grid on which the contour will be plotted
      x_cont <- seq(x.min, x.max, length.out = 100 )
      y_cont <- seq(y.min, y.max, length.out = 100 )
      new_grid <- data.frame( expand.grid(x_cont, y_cont ) )

      data <- data.frame( new_grid, means)
      colnames(data) <- var_names
      # predict
      scored <- predict.isolationForest(forest, data)

      data_fm <- data.frame( cbind( data, scored))
      colnames(data_fm) <- c("x","y","scores")

      ggplot2::ggplot(data_fm, ggplot2::aes( x = x, y = y, z = scores)) +
        ggplot2::geom_raster(ggplot2::aes(fill = scores), interpolate = TRUE) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Anomaly Score")) +
        ggplot2::scale_fill_gradientn(colours= c( "#F5793A","#A95AA1","#85C0F9","#0F2080")) +
        ggplot2::xlab("X") +
        ggplot2::ylab("Y") +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text( angle = 0, vjust = 0.5, size = 12),
          axis.title.x = ggplot2::element_text( size = 12) )
      #  c("#56B4E9","#0072B2","#E69F00", "#D55E00")

  }
}


