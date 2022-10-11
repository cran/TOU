
utils::globalVariables(c("Time", "Adsorption", "Replicate"))
#' @import ggplot2
#' @title A plot with the experimental data and the Transformed
#' Ornstein-Uhlenbeck (TOU) model.
#' @description Plots the experimental data and the fitted TOU model.
#' @param x a tou object.
#' @param time.unit an optional label indicating the unit name for the time, the
#' allowed labels are "seconds", "minutes" and "hours".
#' @param adsorption.unit a label indicating for unit name for the adsorption,
#' the allowed labels are "mg/g" and "mg/mmol".
#' @return a plot with the experimental data and the fitted TOU model.
#' @examples
#'
#' # an example with one trajectory of experimental adsorption
#' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
#'                    250, 300)
#' observed.values <- c(0, 1.684, 2.341, 2.581, 2.842, 2.890, 2.959, 3.042,
#'                      3.083, 3.043, 3.017, 2.954, 2.996, 2.886, 2.844)
#' observed.process <- cbind(observed.time, observed.values)
#' # fitting the model without any fixed parameters
#' result <- fit.model.tou(w=observed.process)
#' # default units time in minutes and default units adsorption in mg/g
#' draw.model.tou(result)
#' # changing units time to seconds and units adsorption to mg/mmol
#' draw.model.tou(result, time.unit="seconds", adsorption.unit="mg/mmol")
#'
#' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
#'                    250, 300)
#' observed.values.1 <- c(0, 1.684, 2.341, 2.581, 2.842, 2.890, 2.959, 3.042,
#'                        3.083, 3.043, 3.017, 2.954, 2.996, 2.886, 2.844)
#' observed.values.2 <- c(0, 1.618, 2.217, 2.571, 2.763, 2.841, 2.866, 2.898,
#'                        2.935, 2.973, 2.906, 2.919, 2.910, 3.012, 3.071)
#' observed.values.3 <- c(0, 1.596, 2.333, 2.611, 2.750, 2.731, 2.829, 2.838,
#'                        2.864, 2.884, 2.886, 2.911, 2.896, 2.877, 2.969)
#' observed.processes <- cbind(observed.time, observed.values.1,
#'                             observed.values.2, observed.values.3)
#' # fitting the model without any fixed parameters
#' result <- fit.model.tou(w=observed.processes)
#' draw.model.tou(result)
#'
#' @export
draw.model.tou <- function(x, time.unit=NULL, adsorption.unit=NULL) {
  n.plot <- nrow(x$observed.data)
  m.plot <- ncol(x$observed.data)
  df.ggplot <- data.frame(Time=rep(0, n.plot*(m.plot-1)), Adsorption=0,
                          Replicate=0)
  df.ggplot$Time <- rep(x$observed.data[,1], m.plot-1)
  for(i in 2:m.plot) {
    a <- (i-2)*n.plot + 1
    b <- (i-1)*n.plot
    df.ggplot$Adsorption[a:b] <- x$observed.data[,i]
    df.ggplot$Replicate[a:b] <- rep(paste("Replicate", i-1), n.plot)
  }
  x.range <- range(df.ggplot$Time)
  x.01 <- x.range[1] + 0.50*(x.range[2]-x.range[1])
  x.02 <- x.range[1] + 0.55*(x.range[2]-x.range[1])
  x.03 <- x.range[1] + 0.46*(x.range[2]-x.range[1])
  x.04 <- x.range[1] + 0.52*(x.range[2]-x.range[1])
  y.range <- range(df.ggplot$Adsorption)
  y.01 <- y.range[1] + 0.10*(y.range[2]-y.range[1])
  y.02 <- y.range[1] + 0.05*(y.range[2]-y.range[1])
  p <- ggplot2::ggplot(data=df.ggplot, ggplot2::aes(x=Time, y=Adsorption, colour=Replicate)) +
       ggplot2::geom_point() +
       ggplot2::theme_bw() +
       ggplot2::theme(legend.position = c(0.8, (35/300)+(1/30)*m.plot)) +
       ggplot2::labs(subtitle="Fitted model in continuous line, experimental data in points",
                     color="Experiment") +
       ggplot2::stat_function(fun=mean.tou,
                              args=list(parameters=as.numeric(c(x$qe,x$lambda,x$a))),
                              colour="black")
  if(is.null(time.unit)) {
    p <- p + ggplot2::xlab("Time")
  } else if(time.unit=="seconds") {
    p <- p + ggplot2::xlab("Time (seconds)")
  } else if(time.unit=="minutes") {
    p <- p + ggplot2::xlab("Time (minutes)")
  } else if(time.unit=="hours") {
    p <- p + ggplot2::xlab("Time (hours)")
  }
  if(is.null(adsorption.unit)) {
    p <- p + ggplot2::ylab( expression(q[t]) )
  } else if(adsorption.unit=="mg/g") {
    p <- p + ggplot2::ylab( expression( paste( q[t], " (mg/g)" ) ) )
  } else if(adsorption.unit=="mg/mmol") {
    p <- p + ggplot2::ylab( expression( paste( q[t], " (mg/mmol)" ) ) )
  }
  plot(p)
}

