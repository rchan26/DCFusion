d2product_norm <- function(x, mu1, mu2, sd1, sd2) {
   mu <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
   std_dv <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
   return(dnorm(x, mean = mu, sd = std_dv))
}

d3product_norm <- function(x, mu1, mu2, mu3, sd1, sd2, sd3) {
   mu12 <- (sd1*sd1*mu2 + sd2*sd2*mu1) / (sd1*sd1 + sd2*sd2)
   std_dv12 <- sqrt(1 / ((1/(sd1*sd1)) + (1/(sd2*sd2))))
   mu123 <- (std_dv12*std_dv12*mu3 + sd3*sd3*mu12) / (std_dv12*std_dv12 + sd3*sd3)
   std_dv123 <- sqrt(1 / ((1/(std_dv12*std_dv12)) + (1/(sd3*sd3))))
   return(dnorm(x, mean = mu123, sd = std_dv123))
}

# setting parameter values for each of the sub-posteriors
mean_values <- c(-1, 1, 8)
sd_values <- sqrt(c(1, 0.5, 4))

##################################################

# sub-posteriors
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', 
      n = 10000, col = '#003262', lwd = 3, lty = 2)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T,
      n = 10000, col = '#a41034', lwd = 3, lty = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, 
      n = 10000, col = '#fcb515', lwd = 3, lty = 1)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]), 
                                   expression(f[3])),
       col = c('#003262', '#a41034', '#fcb515'), 
       lwd = 3, 
       lty = c(2, 3, 1),
       cex = 1.1,
       bty = 'n')

# Method 1
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf',
      n = 10000, lty = 3, lwd = 3)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T,
      n = 10000, lty = 3, lwd = 3)
curve(d2product_norm(x, mean_values[1], mean_values[2], sd_values[1], sd_values[2]), add = T,
      n = 10000, col = '#01653a', lty = 1, lwd = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, 
      n = 10000, col = '#9e292f', lty = 2, lwd = 3)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[2])),
       col = c('black', 'black', '#9e292f', '#01653a'), 
       lwd = 3, 
       lty = c(3, 3, 2, 1),
       cex = 1.1,
       bty = 'n')

# Method 2
curve(dnorm(x, mean_values[1], sd_values[1]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf', 
      n = 10000, lty = 3, lwd = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T,
      n = 10000, lty = 3, lwd = 3)
curve(d2product_norm(x, mean_values[1], mean_values[3], sd_values[1], sd_values[3]), add = T, 
      n = 10000, col = '#01653a', lty = 1, lwd = 3)
curve(dnorm(x, mean_values[2], sd_values[2]), add = T, 
      n = 10000, col = '#9e292f', lty = 2, lwd = 3)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[1]*f[3])),
       col = c('black', '#9e292f', 'black', '#01653a'),
       lwd = 3,
       lty = c(3, 2, 3, 1),
       cex = 1.1,
       bty = 'n')

# Method 3
curve(dnorm(x, mean_values[2], sd_values[2]), -5, 15, ylim = c(0, 0.8), ylab = 'pdf',
      n = 10000, lty = 3, lwd = 3)
curve(dnorm(x, mean_values[3], sd_values[3]), add = T, 
      n = 10000, lty = 3, lwd = 3)
curve(d2product_norm(x, mean_values[2], mean_values[3], sd_values[2], sd_values[3]), add = T, 
      n = 10000, col = '#01653a', lty = 1, lwd = 3)
curve(dnorm(x, mean_values[1], sd_values[1]), add = T, 
      n = 10000, col = '#9e292f', lty = 2, lwd = 3)
legend(x = -5, y = 0.8, legend = c(expression(f[1]), 
                                   expression(f[2]),
                                   expression(f[3]),
                                   expression(f[2]*f[3])),
       col = c('#9e292f', 'black', 'black', '#01653a'), 
       lwd = 2, 
       lty = c(2, 3, 3, 1),
       cex = 1.1,
       bty = 'n')
