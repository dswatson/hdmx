# Register cores
library(doMC)
registerDoMC(4)

# Brute force
c_imp <- function(x, y, type = 'regression', weights = TRUE,
                  mtry = floor(sqrt(ncol(x))), B = 2000, 
                  n.cores = 4, seed = NULL) {
  
  # Define drop_1 function
  drop_1 <- function(j) {
    
    # Drop j
    d0 <- d[, -j]
    
    # Build models, calculate loss (and optionally wts)
    if (type == 'regression') {
      if (weights) {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                     num.threads = 1, seed = seed)
        wts0 <- 1 / predict(f0, data = d0, num.threads = 1,
                            type = 'se')$se^2
      } else {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, 
                     num.threads = 1, seed = seed)
      }
      loss0 <- (f0$predictions - d[, 'y'])^2
    } else if (type %in% c('classification', 'probability')) {
      if (weights) {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                     num.threads = 1, seed = seed, 
                     probability = TRUE) 
        wts0 <- 1 / (predict(f0, data = d0, num.threads = 1,
                             type = 'se')$se[, 1])^2
      } else {
        f0 <- ranger(data = d0, dependent.variable.name = 'y', 
                     mtry = mtry, num.trees = B, 
                     num.threads = 1, seed = seed, 
                     probability = TRUE) 
      }
      y_hat0 <- sapply(f0$predictions, '[[', 1)
      loss0 <- -(d[, 'y'] * log(y_hat0) + (1 - d[, 'y']) * log(1 - y_hat0))
    }
    
    # Perform paired one-sided t-test with optional precision weights
    if (weights) {
      df <- data.frame(Loss = c(loss, loss0),
                       Wts = c(wts, wts0),
                       Model = rep(c('Full', 'Null'), each = n),
                       Obs = rep(paste0('n', seq_len(n)), times = 2))
      s <- summary(lm(Loss ~ Model + Obs, weights = Wts, data = df))
      out <- c(coef(s)[2, 1:3], pt(coef(s)[2, 3], s$df[2], lower = FALSE))
    } else {
      t_test <- t.test(loss0, loss, paired = TRUE, alternative = 'greater')
      out <- c(t_test$estimate, t_test$estimate / t_test$statistic,
               t_test$statistic, t_test$p.value)
    }
    
    # Export
    return(out)
  }
  
  # Prelimz
  require(ranger)
  n <- nrow(x)
  p <- ncol(x)
  d <- cbind(x, 'y' = y)
  
  # Define test
  if (type == 'regression') {
    if (weights) {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                  num.threads = n.cores, seed = seed)
      # Calculate precision
      wts <- 1 / predict(f, data = d, num.threads = n.cores,
                         type = 'se')$se^2
    } else {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, 
                  num.threads = n.cores, seed = seed)
    }
    # Calculate sample-wise loss
    loss <- (f$predictions - y)^2
  } else if (type %in% c('classification', 'probability')) {
    if (weights) {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, keep.inbag = TRUE, 
                  num.threads = n.cores, seed = seed,
                  probability = TRUE)
      # Calculate sample-wise precision
      wts <- 1 / predict(f, data = d, num.threads = n.cores,
                         type = 'se')$se[, 1]^2
    } else {
      # Grow full forest
      f <- ranger(data = d, dependent.variable.name = 'y', 
                  mtry = mtry, num.trees = B, 
                  num.threads = n.cores, seed = seed,
                  probability = TRUE)
    }
    # Calculate sample-wise loss
    y_hat <- sapply(f$predictions, '[[', 1)
    loss <- -(y * log(y_hat) + (1 - y) * log(1 - y_hat))
  }
  
  # Execute in parallel
  delta <- foreach(j = seq_len(p), .combine = rbind) %dopar% drop_1(j)
  dimnames(delta) <- list(NULL, c('Delta', 'SE', 't', 'p.value'))
  return(data.frame(Feature = colnames(x), delta))
  
}

