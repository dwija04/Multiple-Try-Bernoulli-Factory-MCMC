library(LaplacesDemon)

# target is a gamma mixture of weibulls
C_f <- function(val, k = 10, mu )
{
  ret <- k/(exp(1)*val*dnorm(val, mean = mu, sd = sqrt(0.001)) )
  ret[ret<0] <- 0
  return(ret)
}

p_f <- function(val, k = 10, alpha = 10, beta = 0.01)
{
  m <- rgamma(1, shape = alpha, scale = beta)
  return(dweibull(val, shape = k, scale = m) / (k/(exp(1)*val)) )
}

# m is the number of tries
# k is the shape parameter of the weibull
# init is the initial value

bern_cat <- function(c, tries, pre)
{
  accept <- 0
  try <- 0
  while(!accept)
  {
    try <- try + 1
    if(max(c) == 0)
    {
      return(c(0, try))
    }
    c_new <- c / sum(c)
    c1 <- rcat(1, c_new)
    y <- tries[c1]
    
    c_theta <- c[c1]
    p_theta <- p_f(y, k = 10, alpha = 10, beta = 0.01)
    
    c2 <- rbern(1, p_theta)
    if(c2 == 1)
    {
      accept <- 1
      return(c(y,try))
    }
    
  }
  
}

mtm_bern <- function(n , m = 10, init = 0.1, k = 10, beta = 1)
{
  x <- numeric(n)
  x[1] = init
  loops <- numeric(n)
  loops.outer <- numeric(n)
  loops[1] <- 0
  loops.outer[1] <- 0
  
  aux <- numeric(m)
  for(i in 2:n)
  {
    #sampling y_j
    tries <- rnorm(m, mean = x[i-1], sd = sqrt(0.001))
    c <- C_f(tries, k = 10, mu = x[i-1])
    c[c < 0] <- 0
    
    #sorting the bounds, and the trials according to the indices of the bounds
    # sorted_indices <- order(c)
    # sorted_c <- c[sorted_indices]
    # rearranged_tries <- tries[sorted_indices]
    # tries <- rearranged_tries
    # c <- sorted_c
    
    # choosing y 
    foo <- bern_cat(c = c, tries = tries,  pre = x[i-1]) #choosing a proposal
    y <- foo[1]
    if(y == 0)
    {
      x[i] <- x[i-1]
      loops[i] <- foo[2]
      next
    }
    
    loops[i] <- foo[2] 
    
    #sampling x_j^*
    
    aux[1] = x[i-1]
    aux[2:m] <- rnorm(m-1, mean = y, sd = sqrt(0.001))
    
    
    #ratio = a/a+b
    c_tries <- c
    c_aux <- C_f(aux, k = 10, mu = y)
    c_aux[c_aux < 0] <- 0
    
    c_y <- m*max(c_tries)
    c_x <- m*max(c_aux)
    C <- c_y / (c_y + c_x)
    
    #portkey to accept or reject y
    accept <- 0
    loops.outer[i] <- 0
    while(!accept)
    {
      if(rbern(1, beta) == 0)
      {
        x[i] <- x[i-1]
        accept <- 1
        break
      }
      else
      {
        loops.outer[i] <- loops.outer[i] + 1
        C1 <- rbinom(1, 1, C)
        if(C1 == 1)
        {
          
          lam <- rgamma(1, shape = 10, scale = 0.01)
          A <- sum(dweibull(tries, shape = k, scale = lam) / dnorm(tries, mean = x[i-1], sd = sqrt(0.001)))
          p_y <- A / c_y
          C2 <- rbinom(1, 1, p_y)
          if(C2 == 1)
          {
            x[i] <- y
            accept <- 1
          }
          
        }
        else
        {
          lam <- rgamma(1, shape = 10, scale = 0.01)
          A <- sum(dweibull(aux, shape = k, scale = lam) / dnorm(aux, mean = y, sd = sqrt(0.001)))
          p_x <- A / c_x
          C2 <- rbinom(1, 1, p_x)
          if(C2 == 1)
          {
            x[i] <- x[i-1]
            accept <- 1
          }
        }
      }
    }
  }
  # avg_loops <- mean(loops)
  # avg_loops2 <- mean(loops2)
  print('done')
  return(cbind(x, loops, loops.outer))
}

foo <- mtm_bern(1e3, m = 100, beta = 1)
 
samples <- foo[, 1]
inner.l <- foo[, 2]
outer.l <- foo[, 3]
par(mfrow = c(2,2))
  
plot(density(samples))
acf(samples)
plot.ts(samples)

mean(samples)
sd(samples)
# plot(inner.l)
# plot(outer.l)


# z <- replicate(20, mtm_bern(1e4, m = 5, beta = 0.99))

mean_loops_1 <- mean(z[, 2, ])
mean_loops_2 <- mean(z[, 3, ])
max_loops_1 <- max(z[, 2, ])
max_loops_2 <- max(z[, 3, ])
    
print(mean_loops_1)  
print(mean_loops_2)
print(max_loops_1)
print(max_loops_2)

chain <- z[, 1, 13]
plot(density(chain))
acf(chain)
# mean(outer.l + inner.l)