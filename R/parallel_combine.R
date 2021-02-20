######################### conensus Monte Carlo (Scott et al. 2016) #########################

#' @export
consensus_scott <- function(S, 
                            samples_to_combine, 
                            indep = FALSE,
                            shuff = FALSE,
                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  pcm <- proc.time()
  # check all matrices in samples_to_combine are of the same length
  dimension <- ncol(samples_to_combine[[1]])
  num_samples <- nrow(samples_to_combine[[1]])
  if (!(is.list(samples_to_combine)) && length(samples_to_combine)!=S) {
    stop('consensus_Scott: samples_to_combine must be a list of length S')
  } else {
    for (i in 2:S) {
      if (any(dim(samples_to_combine[[i]]) != c(num_samples, dimension))) {
        stop('consensus_Scott: samples_to_combine must contain matrices of the same size')
      }
    }
  }
  
  # define array to store the samples to use the parallelMCMCcombine package
  batch_samples <- array(NA, dim = c(dimension, num_samples, S))
  for (s in 1:S) {
    samples_to_combine[[s]]
    batch_samples[,,s] <- t(samples_to_combine[[s]])
  }
  
  # use appropriate method depending if parameters are independent or not
  if (!indep) {
    consensus_samples <- parallelMCMCcombine::consensusMCcov(subchain = batch_samples, 
                                                             shuff = shuff)
  } else {
    consensus_samples <- parallelMCMCcombine::consensusMCindep(subchain = batch_samples,
                                                               shuff = shuff)
  }
  
  final <- proc.time()-pcm
  # return samples as an n x d matrix
  return(list('samples' = t(consensus_samples), 
              'time' = final['elapsed']))
}

######################### semi-parametric method (KDEMC) (Neiswanger et al. 2013) #########################

#' @export
neiswanger <- function(S, 
                       samples_to_combine,
                       bw = NULL, 
                       anneal = TRUE,
                       shuff = FALSE,
                       seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  pcm <- proc.time()
  # check all matrices in samples_to_combine are of the same length
  dimension <- ncol(samples_to_combine[[1]])
  num_samples <- nrow(samples_to_combine[[1]])
  if (!(is.list(samples_to_combine)) && length(samples_to_combine)!=S) {
    stop('neiswanger: samples_to_combine must be a list of length S')
  } else {
    for (i in 2:S) {
      if (any(dim(samples_to_combine[[i]]) != c(num_samples, dimension))) {
        stop('neiswanger: samples_to_combine must contain matrices of the same size')
      }
    }
  }
  
  if (is.null(bw)) {
    bw <- rep(1, dimension)
  } else if (length(bw)!=ncol(samples_to_combine[[1]])) {
    stop('neiswanger: bw must be a vector of length ncol(samples_to_combine[[1]])')
  }
  
  # define array to store the samples to use the parallelMCMCcombine package
  batch_samples <- array(NA, dim = c(dimension, num_samples, S))
  for (s in 1:S) {
    batch_samples[,,s] <- t(samples_to_combine[[s]])
  }
  
  # use appropriate method depending if parameters are independent or not
  neiswanger_samples <- parallelMCMCcombine::semiparamDPE(subchain = batch_samples,
                                                          bandw = bw,
                                                          anneal = anneal, 
                                                          shuff = shuff)
  
  final <- proc.time()-pcm
  # return samples as an n x d matrix
  return(list('samples' = t(neiswanger_samples), 
              'time' = final['elapsed']))
}

######################### weierstrass sampler (Wang and Dunson 2013) #########################

# taken from: https://github.com/david-dunson/weierstrass/blob/master/R/weierstrass.r
repmat <- function(mat, nrow = 1, ncol = 1) {
  if(is.null(dim(mat))){
    cat('it is not a matrix')
    break
  }
  
  r = dim(mat)[1]
  c = dim(mat)[2]
  return(matrix(rep(mat, nrow*ncol), nrow = nrow*r, ncol = ncol*c))
}

# taken from: https://github.com/david-dunson/weierstrass/blob/master/R/weierstrass.r
#' @export
weierstrass <- function(Samples, 
                        num.sets = NULL, 
                        para.dim = NULL,
                        method = 'reject',
                        kernel = 'uniform', 
                        kernel.corr = TRUE,
                        accept = 0.1, 
                        weight = TRUE,
                        average = TRUE, 
                        matching = TRUE, 
                        resampling = TRUE, 
                        resampling.corr = FALSE, 
                        display = TRUE,
                        level = 1, 
                        sets = 1,
                        shuff = FALSE,
                        seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  pcm <- proc.time()
  m = length(Samples)
  if(!is.null(num.sets) && m!=num.sets){
    cat('Number of sets not consistent')
    break
  }
  if(is.null(dim(Samples[[1]])))
    p = 1
  else
    p = dim(Samples[[1]])[2]
  
  # if we want to shuffle around the sub-posterior samples
  if (shuff) {
    if (p==1) {
      for (i in 1:length(Samples)) {
        Samples[[i]] <- Samples[[i]][sample(length(Samples[[i]]))]
      }
    } else {
      for (i in 1:length(Samples)) {
        Samples[[i]] <- Samples[[i]][sample(nrow(Samples[[i]])),]
      }
    }
  }
  
  if(!is.null(para.dim) && p!=para.dim){
    cat('Number of parameters not consistent')
    break
  }
  if(p == 1)
    kernel.corr = FALSE
  if(m == 1)
    return(Samples[[1]])
  else if(m == 2){
    X1 = matrix(Samples[[1]],ncol = p)
    X2 = matrix(Samples[[2]],ncol = p)
    N = min(dim(X1)[1], dim(X2)[1])
    
    X1 = matrix(X1[1:N,],N,p)
    X2 = matrix(X2[1:N,],N,p)
    if(matching == TRUE) matching = N
    
    output = matrix(0, N, p)
    X = NULL
    w1 = w = matrix(1/2, p, 2)
    sigmah = matrix(0, p, 1)
    for(j in 1:p){
      w1[j,] = c(1/var(X1[,j]), 1/var(X2[,j]))
      sigmah[j,1] =1/(1/var(X1[,j])+1/var(X2[,j]))
      w1[j,] = w1[j,]/sum(w1[j,])
    }
    if(weight == TRUE) w = w1
    else sigmah = matrix(rep(mean(sigmah),p),p,1)
    
    if(kernel.corr == TRUE){
      #print(dim(X1))
      invC1 = solve(cov(X1))
      invC2 = solve(cov(X2))
      newC = solve(invC1 + invC2)
      if(weight == TRUE){
        W1 = invC1%*%newC
        W2 = invC2%*%newC
      }else{
        W1 = W2 = diag(rep(1/2,p), p,p)
      }
    }
    
    if(method == 'reject'){
      if(matching != FALSE)
        iter = ceiling(matching/N/accept)
      else
        iter = 1
      for(k in 1:iter){
        permut = sample(1:N, N)
        X2 = matrix(X2[permut,],N,p)
        if(average == TRUE){
          if(kernel.corr == TRUE){
            #print(dim(X1), dim(X2), dim(W1), dim(W2))
            output = X1%*%W1 + X2%*%W2
          }
          else
            output = X1*repmat(t(w[,1]),N,1) + X2*repmat(t(w[,2]),N,1)
        }
        else{
          judge = matrix(runif(N*p), N, p)<repmat(t(w[,1]),N,1)
          output = X2
          output[judge] = X1[judge]
        }
        
        if(kernel == 'uniform'){
          if(kernel.corr == TRUE)
            diff = abs(X1 - X2)%*%(invC1+invC2)
          else
            diff = abs(X1 - X2)/repmat(t(sigmah),N,1)
          max.diff = apply(diff, 1, max)
          cutoff = quantile(max.diff, prob = accept)
          keep = max.diff<=cutoff
          #print(length(keep))
        }
        
        if(kernel == 'gaussian'){
          if(kernel.corr == TRUE){
            max.diff = (X1%*%W1%*%(invC1+invC2)*X1 + X2%*%W2%*%(invC1+invC2)*X2
                        - output%*%(invC1+invC2)*output)%*%matrix(1,p,1)
          }
          else{
            diff = (X1^2*repmat(t(w[,1]),N,1) + X2^2*repmat(t(w[,2]),N,1)
                    - output^2)/repmat(t(sigmah^2),N,1)
            max.diff = diff%*%matrix(1,p,1)
          }
          #find h
          minh = 0
          maxh = 1
          while(mean(exp(-max.diff/maxh))<accept){
            minh = maxh
            maxh = maxh*2
          }
          while(abs(maxh-minh)>0.01){
            h = (minh + maxh)/2
            if(mean(exp(-max.diff/h))<accept)
              minh = h
            else
              maxh = h
          }
          h = (minh + maxh)/2
          #print(sum(exp(-max.diff/h)))
          keep = (runif(N) - exp(-max.diff/h)) <=0
        }
        
        X = rbind(X, matrix(output[keep,], sum(keep), p))
      }
    }
    
    if(method == 'importance'){
      if(matching != FALSE){
        iter = ceiling(matching/N/accept)
        accept0 = accept
      }
      else{
        iter = 1
        accept0 = 1
      }
      for(k in 1:iter){
        permut = sample(1:N, N)
        X2 = matrix(X2[permut,],N,p)
        if(kernel.corr == TRUE){
          output = X1%*%W1 + X2%*%W2
          max.diff = (X1%*%W1%*%(invC1+invC2)*X1 + X2%*%W2%*%(invC1+invC2)*X2
                      - output%*%(invC1+invC2)*output)%*%matrix(1,p,1)
        }else{
          output = X1*repmat(t(w[,1]),N,1) + X2*repmat(t(w[,2]),N,1)
          output2 = X1^2*repmat(t(w[,1]),N,1) + X2^2*repmat(t(w[,2]),N,1)
          diff = (output2 - output^2)/repmat(t(sigmah), N, 1)
          max.diff = diff%*%matrix(1,p,1)
        }
        
        if(k==1){
          #find h
          minh = 0
          maxh = 1
          
          while(mean(exp(-max.diff/maxh))<accept){
            minh = maxh
            maxh = maxh*2
          }
          while(abs(maxh-minh)>0.01){
            h = (minh + maxh)/2
            if(mean(exp(-max.diff/h))<accept)
              minh = h
            else
              maxh = h
          }
        }
        w2 = exp(-max.diff/h)
        w2 = w2/sum(w2)
        keep = sample(1:N, size = ceiling(N*accept0), replace = TRUE, prob = w2)
        temp = matrix(output[keep,], length(keep), p)
        if(resampling == TRUE){
          require(MASS)
          if(resampling.corr == TRUE){
            for(it in 1:length(keep)){
              temp[it,] = mvrnorm(n=1,mu = temp[it,],Sigma = newC*sqrt(h)/8)
            }
          }
          else{
            temp = matrix(rnorm(p*length(keep), temp, repmat(t(sigmah)*sqrt(h)/8)), length(keep), p)
          }
        }
        X = rbind(X, temp)
      }
    }
    return(X)
  }
  else if(m==3){
    Comb1 = weierstrass(Samples[1:2], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept, weight = weight, average = average, matching = TRUE, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level + 1)
    
    twoSamples = list()
    if (is.list(Comb1)) {
      twoSamples[[1]] = Comb1$samples
    } else {
      twoSamples[[1]] = Comb1
    }
    twoSamples[[2]] = Samples[[3]]
    
    FinalComb = weierstrass(twoSamples, num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level)
    
    final <- proc.time()-pcm
    return(list('samples' = FinalComb, 
                'time' = final['elapsed']))
  }
  else{
    if(matching == FALSE)
      accept1 = sqrt(accept)
    else
      accept1 = accept
    if(display == TRUE)
      cat(paste('starting Level',toString(level),', Set',toString(sets),'-',toString(sets-1+ceiling(m/2)),'\n'))
    Comb1 = weierstrass(Samples[1:ceiling(m/2)], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level + 1, sets = sets)
    
    if(display == TRUE)
      cat(paste('Level',toString(level),', Set',toString(sets),'-',toString(sets-1+ceiling(m/2)),'completed. Starting set', toString(sets-1+ceiling(m/2)+1), '-', toString(sets-1+m),'\n'))
    Comb2 = weierstrass(Samples[(ceiling(m/2)+1):m], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr= resampling.corr, display = display, level = level + 1, sets = sets+ceiling(m/2))
    
    if(display == TRUE)
      cat(paste('Level',toString(level),', Set',toString(sets-1+ceiling(m/2)+1),'-',toString(sets-1+m),'completed. Combining two final sets.\n'))
    twoSamples = list()
    if (is.list(Comb1)) {
      twoSamples[[1]] = Comb1$samples
    } else {
      twoSamples[[1]] = Comb1
    }

    if (is.list(Comb2)) {
      twoSamples[[2]] = Comb2$samples
    } else {
      twoSamples[[2]] = Comb2
    }
    
    
    FinalComb = weierstrass(twoSamples, num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level)
    if(display == TRUE)
      cat(paste('Level',toString(level),'completed, retrieving upper level.\n'))
    
    final <- proc.time()-pcm
    return(list('samples' = FinalComb, 
                'time' = final['elapsed']))
  }
}

#' #' @export
#' plot_consensus_results <- function(plot_rows, plot_columns, full_post, consensus_mat, consensus_sca, 
#'                                    ylimit = NULL, dimensions = FALSE) {
#'   if (dimensions) {
#'     # first plot each of the dimensions 
#'     par(mfrow=c(1, 3))
#'     for (i in 1:ncol(full_post)) {
#'       plot(density(full_post[,i]), xlab = paste(expression(beta),i-1), main = 'full')
#'       abline(v=mean(full_post[,i]), lty = 4)
#'       plot(density(consensus_mat[,i]), xlab = paste(expression(beta),i-1), col = 'red', lty = 2, 'consensus mat')
#'       abline(v=mean(full_post[,i]), lty = 4)
#'       plot(density(consensus_sca[,i]), xlab = paste(expression(beta),i-1), col = 'blue', lty = 4, 'consensus scale')
#'       abline(v=mean(full_post[,i]), lty = 4)
#'     }
#'   }
#'   par(mfrow=c(plot_rows, plot_columns))
#'   for (i in 1:ncol(full_post)) {
#'     if (is.null(ylimit)) {
#'       plot(density(full_post[,i]), xlab = paste(expression(beta),i-1), main = '')
#'       lines(density(consensus_mat[,i]), col = 'red', lty = 2)
#'       lines(density(consensus_sca[,i]), col = 'blue', lty = 4)
#'     } else {
#'       if (length(ylimit)!=ncol(full_post)) {
#'         stop('plot_consensus_results: ylimit must be a vector of length ncol(full_post)')
#'       }
#'       print(ylimit[i])
#'       plot(density(full_post[,i]), xlab = paste(expression(beta),i-1), main = '', ylim = c(0, ylimit[i]))
#'       lines(density(consensus_mat[,i]), col = 'red', lty = 2)
#'       lines(density(consensus_sca[,i]), col = 'blue', lty = 4)
#'     }
#'   }
#'   par(mfrow=c(1,1))
#' }
