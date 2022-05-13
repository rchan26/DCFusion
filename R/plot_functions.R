#' Acceptance rate comparison between results
#'
#' Creates plots to compare the acceptance probabilities and time taken to run
#'
#' @param hier1 hierarchical/progressive Monte Carlo fusion result 
#'              (will be in blue on plots by default)
#' @param hier2 hierarchical/progressive Monte Carlo fusion result 
#'              (will be in green on plots by default)
#' @param colours vector of 2 colours, where colours[1] is for hier1 and 
#'                colours[2] is for hier2
#' @param time time that you are comparing, if not appropriate, can use NA
#' @param hierarchical logical value to determine if it was a result from 
#'                     hierarchical fusion (TRUE) or progressive fusion (FALSE)
#'
#' return None
#'
#' @export
acceptance_rate_plots <- function(hier1, 
                                  hier2,
                                  colours = c('blue', 'green'), 
                                  time = NA, 
                                  hierarchical) {
  par(mfrow=c(2,2), oma = c(0, 0, 4, 0))
  if (hierarchical) {
    #rho
    plot(1:length(hier1$overall_rho), hier1$overall_rho, ylim = c(0,1), col = colours[1],
         ylab = 'Rho Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$overall_rho), hier1$overall_rho, col = colours[1])
    points(1:length(hier2$overall_rho), hier2$overall_rho, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$overall_rho), hier2$overall_rho, col = colours[2])
    # Q
    plot(1:length(hier1$overall_Q), hier1$overall_Q, ylim = c(0,1), col = colours[1],
         ylab = 'Q Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$overall_Q), hier1$overall_Q, col = colours[1])
    points(1:length(hier2$overall_Q), hier2$overall_Q, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$overall_Q), hier2$overall_Q, col = colours[2])
    # rho*Q (overall)
    plot(1:length(hier1$overall_rhoQ), hier1$overall_rhoQ, ylim = c(0,1), col = colours[1],
         ylab = 'Overall Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$overall_rhoQ), hier1$overall_rhoQ, col = colours[1])
    points(1:length(hier2$overall_rhoQ), hier2$overall_rhoQ, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$overall_rhoQ), hier2$overall_rhoQ, col = colours[2])
    # time
    plot(1:length(hier1$overall_time), hier1$overall_time, col = colours[1],
         ylab = 'Time Taken (sec)', xlab = 'Level',
         ylim = c(min(c(hier1$overall_time, hier2$overall_time)), 
                  max(c(hier1$overall_time, hier2$overall_time))), 
         pch = 4, cex = 0.5)
    lines(1:length(hier1$overall_time), hier1$overall_time, col = colours[1])
    points(1:length(hier2$overall_time), hier2$overall_time, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$overall_time), hier2$overall_time, col = colours[2])
  } else {
    #rho
    plot(1:length(hier1$rho_acc), hier1$rho_acc, ylim = c(0,1), col = colours[1],
         ylab = 'Rho Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$rho_acc), hier1$rho_acc, col = colours[1])
    points(1:length(hier2$rho_acc), hier2$rho_acc, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$rho_acc), hier2$rho_acc, col = colours[2])
    # Q
    plot(1:length(hier1$Q_acc), hier1$Q_acc, ylim = c(0,1), col = colours[1],
         ylab = 'Q Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$Q_acc), hier1$Q_acc, col = colours[1])
    points(1:length(hier2$Q_acc), hier2$Q_acc, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$Q_acc), hier2$Q_acc, col = colours[2])
    # rho*Q (overall)
    plot(1:length(hier1$rhoQ_acc), hier1$rhoQ_acc, ylim = c(0,1), col = colours[1],
         ylab = 'Overall Acceptance Rate', xlab = 'Level', pch = 4, cex = 0.5)
    lines(1:length(hier1$rhoQ_acc), hier1$rhoQ_acc, col = colours[1])
    points(1:length(hier2$rhoQ_acc), hier2$rhoQ_acc, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$rhoQ_acc), hier2$rhoQ_acc, col = colours[2])
    # time
    plot(1:length(hier1$time), hier1$time, col = colours[1],
         ylab = 'Time Taken (sec)', xlab = 'Level',
         ylim = c(min(c(hier1$time, hier2$time)), max(c(hier1$time, hier2$time))), 
         pch = 4, cex = 0.5)
    lines(1:length(hier1$time), hier1$time, col = colours[1])
    points(1:length(hier2$time), hier2$time, col = colours[2], pch = 4, cex = 0.5)
    lines(1:length(hier2$time), hier2$time, col = colours[2])
  }
  title(main = paste('Time =', time), outer = T)
  par(mfrow=c(1,1))
}

#' @export
plot_fusion_results <- function(plot_rows, 
                                plot_columns, 
                                full_post, 
                                fusion_post, 
                                ylimit = NULL, 
                                dimensions = FALSE,
                                bw = NULL) {
  original_par <- par()
  if (is.null(bw)) {
    bw <- rep("nrd0", ncol(full_post))
  } else if (length(bw)!=ncol(full_post)) {
    stop('plot_fusion_results: bw must be a vector of length ncol(full_post)')
  }
  if (dimensions) {
    # first plot each of the dimensions 
    par(mfrow=c(1, 2))
    for (i in 1:ncol(full_post)) {
      plot(density(full_post[,i], bw = bw[i]), xlab = paste(expression(beta),i-1), main = 'full')
      abline(v=mean(full_post[,i]), lty = 4)
      plot(density(fusion_post[,i], bw = bw[i]), xlab = paste(expression(beta),i-1), col = 'red', lty = 2, 'fusion')
      abline(v=mean(full_post[,i]), lty = 4)
    }
  }
  par(mfrow=c(plot_rows, plot_columns))
  for (i in 1:ncol(full_post)) {
    if (is.null(ylimit)) {
      plot(density(full_post[,i], bw = bw[i]), xlab = paste(expression(beta),i-1), main = '')
      lines(density(fusion_post[,i], bw = bw[i]), col = 'red', lty = 2)
    } else {
      if (length(ylimit)!=ncol(full_post)) {
        stop('plot_fusion_results: ylimit must be a vector of length ncol(full_post)')
      }
      plot(density(full_post[,i], bw = bw[i]), xlab = paste(expression(beta),i-1), main = '', ylim = c(0, ylimit[i]))
      lines(density(fusion_post[,i], bw = bw[i]), col = 'red', lty = 2)
    }
  }
  # reset plots
  tryCatch(par(original_par), warning = function(w) {})
}

#' @export
plot_fusion_matrix <- function(full_post, 
                               fusion_post, 
                               common_limit, 
                               title = NULL,
                               bw = NULL) {
  original_par <- par()
  if (ncol(full_post)!=ncol(fusion_post)) {
    stop('plot_fusion_matrix: full_post and fusion_post must be matrices of same dimension')
  }
  # create dimensions of plot
  dimensions <- ncol(full_post)
  # if title is provided, we need to have more space on the top to accommodate
  if (!is.null(title)) {
    par(mfrow=c(dimensions, dimensions), 
        mai = c(0.25, 0.25, 0.25, 0.25),
        oma = c(0.75, 1, 2.5, 1))
  } else {
    par(mfrow=c(dimensions, dimensions), 
        mai = c(0.25, 0.25, 0.25, 0.25),
        oma = c(0.75, 1, 0.75, 1))
  }
  if (is.null(bw)) {
    bw <- rep("nrd0", dimensions)
  } else if (length(bw)!=dimensions) {
    stop('plot_fusion_matrix: bw must be a vector of length ncol(full_post)')
  }
  for (i in 1:dimensions) {
    for (j in 1:dimensions) {
      # plot the univariate samples
      if (i==j) {
        plot(density(full_post[,i], bw = bw[i]), xlab = '', ylab = '', main = '')
        lines(density(fusion_post[,i], bw = bw[i]), col = 'red', lty = 2)
      } else if (i < j) {
        # get the limits of x and y that just should just fit in the kdes
        xlim_down <- min(min(full_post[,i]), min(fusion_post[,i]))
        xlim_up <- max(max(full_post[,i]), max(fusion_post[,i]))
        ylim_down <- min(min(full_post[,j]), min(fusion_post[,j]))
        ylim_up <- max(max(full_post[,j]), max(fusion_post[,j]))
        # plot kde for full posterior
        plot(ks::kde(full_post[,c(i,j)]), xlab = '', ylab = '', main = '', 
             xlim = c(xlim_down, xlim_up), ylim = c(ylim_down, ylim_up))
        # plot kde for fusion posterior
        plot(ks::kde(fusion_post[,c(i,j)]), col = 'red', add = T)
      } else {
        # plot kde for full posterior
        plot(ks::kde(full_post[,c(i,j)]), xlab = '', ylab = '', main = '', 
             xlim = common_limit, ylim = common_limit)
        # plot kde for fusion posterior
        plot(ks::kde(fusion_post[,c(i,j)]), col = 'red', add = T)
      }
    }
  }
  # plot title if provided
  if (!is.null(title)) mtext(title, outer = TRUE, cex = 1.5)
  # reset plots
  tryCatch(par(original_par), warning = function(w) {})
}

#' @export
compare_samples_univariate <- function(plot_rows, 
                                       plot_columns, 
                                       posteriors, 
                                       colours, 
                                       ylimit = NULL, 
                                       bw = NULL) {
  original_par <- par()
  if (length(posteriors)!=length(colours)) {
    stop('compare_samples_univariate: the length of posteriors and colours must be the same')
  } 
  original_par <- par()
  if (is.null(bw)) {
    bw <- rep("nrd0", ncol(posteriors[[1]]))
  } else if (length(bw)!=ncol(posteriors[[1]])) {
    stop('compare_samples_univariate: bw must be a vector of length ncol(posteriors[[1]])')
  }
  par(mfrow=c(plot_rows, plot_columns))
  for (i in 1:ncol(posteriors[[1]])) {
    if (is.null(ylimit)) {
      plot(density(posteriors[[1]][,i], bw = bw[i]), xlab = paste(expression(beta),i-1), main = '', col = colours[1])
      for (index in 2:length(posteriors)) {
        lines(density(posteriors[[index]][,i], bw = bw[i]), col = colours[index], lty = 2)
      }
    } else {
      if (length(ylimit)!=ncol(posteriors[[1]])) {
        stop('compare_samples_univariate: ylimit must be a vector of length ncol(posteriors[[1]])')
      }
      plot(density(posteriors[[1]][,i], bw = bw[i]), xlab = paste(expression(beta),i-1), main = '', 
           ylim = c(0, ylimit[i]), col = colours[1])
      for (index in 2:length(posteriors)) {
        lines(density(posteriors[[index]][,i], bw = bw[i]), col = colours[index], lty = 2)
      }
    }
  }
  # reset plots
  tryCatch(par(original_par), warning = function(w) {})
}

#' @export
compare_samples_bivariate <- function(posteriors, 
                                      colours, 
                                      common_limit, 
                                      title = NULL,
                                      bw = NULL) {
  original_par <- par()
  if (length(posteriors)!=length(colours)) {
    stop('compare_samples_bivariate: the length of posteriors and colours must be the same')
  } 
  # create dimensions of plot
  dimensions <- ncol(posteriors[[1]])
  # if title is provided, we need to have more space on the top to accommodate
  if (!is.null(title)) {
    par(mfrow=c(dimensions, dimensions), 
        mai = c(0.25, 0.25, 0.25, 0.25),
        oma = c(0.75, 1, 2.5, 1))
  } else {
    par(mfrow=c(dimensions, dimensions), 
        mai = c(0.25, 0.25, 0.25, 0.25),
        oma = c(0.75, 1, 0.75, 1))
  }
  if (is.null(bw)) {
    bw <- rep("nrd0", dimensions)
  } else if (length(bw)!=dimensions) {
    stop('compare_samples_bivariate: bw must be a vector of length ncol(full_post)')
  }
  for (i in 1:dimensions) {
    for (j in 1:dimensions) {
      # plot the univariate samples
      if (i==j) {
        plot(density(posteriors[[1]][,i], bw = bw[i]), xlab = '', ylab = '', main = '', col = colours[1])
        for(index in 2:length(posteriors)) {
          lines(density(posteriors[[index]][,i], bw = bw[i]), col = colours[index], lty = index)
        }
      } else if (i < j) {
        # get the limits of x and y that just should just fit in the kdes
        xlim_down <- min(sapply(posteriors, function(post) min(post[,i])))
        xlim_up <- max(sapply(posteriors, function(post) max(post[,i])))
        ylim_down <- min(sapply(posteriors, function(post) min(post[,j])))
        ylim_up <- max(sapply(posteriors, function(post) max(post[,j])))
        # plot kde for full posterior
        plot(ks::kde(posteriors[[1]][,c(i,j)]), xlab = '', ylab = '', main = '', 
             xlim = c(xlim_down, xlim_up), ylim = c(ylim_down, ylim_up), col = colours[1])
        # plot kde for fusion posterior
        for (index in 2:length(posteriors)) {
          plot(ks::kde(posteriors[[index]][,c(i,j)]), col = colours[index], add = T)
        }
      } else {
        # plot kde for full posterior
        plot(ks::kde(posteriors[[1]][,c(i,j)]), xlab = '', ylab = '', main = '', 
             xlim = common_limit, ylim = common_limit, col = colours[1])
        # plot kde for fusion posterior
        for (index in 2:length(posteriors)) {
          plot(ks::kde(posteriors[[index]][,c(i,j)]), col = colours[index], add = T)
        }
      }
    }
  }
  # plot title if provided
  if (!is.null(title)) mtext(title, outer = TRUE, cex = 1.5)
  # reset plots
  tryCatch(par(original_par), warning = function(w) {})
}

#' @export
integrated_abs_distance <- function(full_post, fusion_post, bw = NULL, print_res = TRUE) {
  if (is.vector(full_post) & is.vector(fusion_post)) {
    full_post <- matrix(full_post)
    fusion_post <- matrix(fusion_post)
  }
  if (ncol(full_post)!=ncol(fusion_post)) {
    stop('integrated_abs_distance: full_post and fusion_post must be matrices of same dimension')
  } 
  dimensions <- ncol(full_post)
  if (is.null(bw)) {
    bw <- rep("nrd0", dimensions)
  } else if (length(bw)!=dimensions) {
    stop('integrated_abs_distance: bw must be a vector of length ncol(full_post)')
  }
  int_abs_dist <- rep(0, dimensions)
  bandwidth_chosen <- rep(0, dimensions)
  for (i in 1:dimensions) {
    min_value <- min(c(full_post[,i], fusion_post[,i]))
    max_value <- max(c(full_post[,i], fusion_post[,i]))
    # calculate kde for ith dimension in posterior samples
    baseline_kde <- density(full_post[,i], bw = bw[i], from = min_value, to = max_value, n = 10000)
    bandwidth_chosen[i] <- baseline_kde$bw
    # calculate kde for ith dimension in fusion samples
    fusion_kde <- density(fusion_post[,i], bw = bandwidth_chosen[i], from = min_value, to = max_value, n = 10000)
    if (!identical(baseline_kde$x, fusion_kde$x)) {
      stop('integrated_abs_distance: the grid for x values are not the same')
    }
    # obtain f(x) value from kdes
    baseline_y <- baseline_kde$y
    fusion_y <- fusion_kde$y
    # calculate differences between baseline_y and fusion_y
    diff <- abs(baseline_y - fusion_y)
    # calculating the total variation
    int_abs_dist[i] <- sfsmisc:::integrate.xy(x = baseline_kde$x, fx = diff)
  }
  IAD <- sum(int_abs_dist)/(2*dimensions)
  if (print_res) {
    print(paste('bandwidths chosen:', bandwidth_chosen))
    print(paste('IAD:', IAD))
  }
  return(IAD)
}

#' @export
integrated_abs_distance_exp_4 <- function(fusion_post, mean = 0, beta = 1, bw = NULL, print_res = TRUE) {
  if (!is.vector(fusion_post) & (length(fusion_post)>= 2)) {
    stop("integrated_abs_distance_exp_4: fusion_post must be a vector with length greater than or equal to 2")
  }
  if (is.null(bw)) {
    bw <- "nrd0"
  }
  # calculate kde for posterior samples
  min_value <- min(c(fusion_post, mean-2/sqrt(beta)))
  max_value <- max(c(fusion_post, mean+2/sqrt(beta)))
  fusion_kde <- density(fusion_post, bw = bw, from = min_value, to = max_value, n = 10000)
  bandwidth_chosen <- fusion_kde$bw
  # obtain f(x) value from kde and compute the target density at the same values
  fusion_y <- fusion_kde$y
  baseline_y <- exp_4_density(x = fusion_kde$x, mean = mean, beta = beta)
  # calculate differences between baseline_y and fusion_y
  diff <- abs(baseline_y - fusion_y)
  # calculating the total variation
  IAD <- 0.5*sfsmisc:::integrate.xy(x = fusion_kde$x, fx = diff)
  if (print_res) {
    print(paste('bandwidths chosen:', bandwidth_chosen))
    print(paste('IAD:', IAD)) 
  }
  return(IAD)
}

#' @export
integrated_abs_distance_uniGaussian <- function(fusion_post, mean = 0, sd = 1, beta, bw = NULL, print_res = TRUE) {
  if (!is.vector(fusion_post) & (length(fusion_post)>= 2)) {
    stop("integrated_abs_distance_uniGaussian: fusion_post must be a vector with length greater than or equal to 2")
  }
  if (is.null(bw)) {
    bw <- "nrd0"
  }
  # calculate kde for posterior samples
  min_value <- min(c(fusion_post, mean-4*sd/sqrt(beta)))
  max_value <- max(c(fusion_post, mean+4*sd/sqrt(beta)))
  fusion_kde <- density(fusion_post, bw = bw, from = min_value, to = max_value, n = 10000)
  bandwidth_chosen <- fusion_kde$bw
  # obtain f(x) value from kde and compute the target density at the same values
  fusion_y <- fusion_kde$y
  baseline_y <- dnorm_tempered(x = fusion_kde$x, mean = mean, sd = sd, beta = beta)
  # calculate differences between baseline_y and fusion_y
  diff <- abs(baseline_y - fusion_y)
  # calculating the total variation
  IAD <- 0.5*sfsmisc:::integrate.xy(x = fusion_kde$x, fx = diff)
  if (print_res) {
    print(paste('bandwidths chosen:', bandwidth_chosen))
    print(paste('IAD:', IAD)) 
  }
  return(IAD)
}

#' @export
integrated_abs_distance_biGaussian <- function(fusion_post, marg_means, marg_sds, bw) {
  dimensions <- ncol(fusion_post)
  if (is.null(bw)) {
    bw <- rep("nrd0", 2)
  } else if (length(bw)!=2) {
    stop('integrated_abs_distance_biGaussian: bw must be a vector of length 2')
  }
  d1_IAD <- integrated_abs_distance_uniGaussian(fusion_post = fusion_post[,1],
                                                mean = marg_means[1],
                                                sd = marg_sds[1],
                                                beta = 1,
                                                bw = bw[1],
                                                print_res = FALSE)
  d2_IAD <- integrated_abs_distance_uniGaussian(fusion_post = fusion_post[,2],
                                                mean = marg_means[2],
                                                sd = marg_sds[2],
                                                beta = 1,
                                                bw = bw[2],
                                                print_res = FALSE)
  IAD <- mean(c(d1_IAD, d2_IAD))
  print(paste('IAD:', IAD))
  return(IAD)
}

#' @export
integrated_abs_distance_multiGaussian <- function(fusion_post, marg_means, marg_sds, bw) {
  dimensions <- ncol(fusion_post)
  if (is.null(bw)) {
    bw <- rep("nrd0", dimensions)
  } else if (length(bw)!=dimensions) {
    stop('integrated_abs_distance_biGaussian: bw must be a vector of length dimensions')
  }
  if (length(marg_means)!=dimensions) {
    stop('integrated_abs_distance_biGaussian: marg_means must be a vector of length dimensions')
  } else if (length(marg_sds)!=dimensions) {
    stop('integrated_abs_distance_biGaussian: marg_sds must be a vector of length dimensions')
  }
  IAD <- mean(sapply(1:dimensions, function(d) {
    integrated_abs_distance_uniGaussian(fusion_post = fusion_post[,d],
                                        mean = marg_means[d],
                                        sd = marg_sds[d],
                                        beta = 1,
                                        bw = bw[d],
                                        print_res = FALSE)
  }))
  print(paste('IAD:', IAD))
  return(IAD)
}
