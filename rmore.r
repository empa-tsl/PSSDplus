# Function for generating probability distributions based on more than 2 modal values and their 
# uncertainties
# 
# Yields a vector of length: number of iterations in the simulation
# 
# Arguments:   - values : vector of modal values
#              - coef : vector of uncertainty coefficients. Defaults to NULL. Either coef or max
#                and min need to be provided
#              - max : maximum of the probability distribution
#              - min : minimum of the probability distribution
#              - N : number of iterations in the simulation
#              - linf : lower truncation limit (useful for example to avoid negative concentrations).
#                Defaults to -Inf for no lower truncation
#              - lsup : upper truncation limit. Defaults to Inf for no upper truncation
# 
# Date of last modification: 03.07.2019
# 
# Associated publication: Systematic consideration of parameter uncertainty and variability in
#                         probabilistic species sensitivity distributions
# 
# Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and Véronique Adam
# 
# Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology,
#            Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
# 
# submitted to Integrated Environmental Assessment and Management in July 2019

# -------------------------------------------------------------------------------------------------



rmore <- function(values,
                  coef = NULL,
                  max = NULL,
                  min = NULL,
                  N,
                  linf = -Inf,
                  lsup = Inf){
  
  # Test input
  stopifnot(any(!is.null(max), !is.null(min), !is.null(coef)))
  
  if(!is.null(max) | !is.null(min)){
    stopifnot(!is.null(max),
              !is.null(min),
              !any(values > max),
              !any(values < min))
  }
  
  require(trapezoid)
  require(mc2d)
  
  # calculate the minimum and maximum of the distribution we are looking for
  if(!is.null(coef)){
    dist.min <- values[which.min(values)]*(1-coef[which.min(values)])
    dist.max <- values[which.max(values)]*(1+coef[which.max(values)])
  } else {
    dist.min <- min
    dist.max <- max
  }
  
  # store sorted values
  sort.values <- sort(values)
  # store unique values
  uni.values <- unique(sort.values)
  # store frequency of unique values
  freq.uni.values <- table(sort.values)
  # smallest value
  val.min <- sort.values[1]
  # largest value
  val.max <- sort.values[length(sort.values)]
  
  
  # calculate the heights of the probability distribution of each segment around each unique value
  height.values <- rep(NA, length(uni.values))
  for(i in 2:(length(uni.values)-1)){
    height.values[i] <- max( freq.uni.values[i]*N/(uni.values[i+1]-uni.values[i]),
                             freq.uni.values[i]*N/(uni.values[i]-uni.values[i-1]) )
  }
  # first
  height.values[1] <- freq.uni.values[1]*N/(uni.values[2]-uni.values[1])
  # last
  height.values[length(height.values)] <-
    freq.uni.values[length(height.values)]*N/( uni.values[length(height.values)]-
                                                 uni.values[length(height.values)-1] )
  
  
  ### LEFT TRIANGLE OF OVERALL SPECIES-SPECIFIC NOEC PROBABILITY DISTRIBUTION
  
  if((length(unique(sort.values)) < length(sort.values)) &
     (val.min == sort.values[2])){ 
    
    # how many values are identical
    n <- length(which(sort.values == val.min))
    # first value that is not identical to the previous one
    val.minP1 <- sort.values[n+1]
    # Height of the first value that is not identical to the previous ones
    h <- height.values[1]
    # Area of the triangular distribution, the first on the left-hand side
    A <- (val.min - dist.min)*h/2
    
    # Triangular distribution on the left-hand side
    # if the smallest value is smaller than the truncation limit, do nothing
    if(dist.min < linf){
      left <- NULL
      
    } else {
      left <- rtrunc("rtriang",
                     n = A,
                     min = dist.min,
                     mode = val.min,
                     max = val.min,
                     lsup = lsup,
                     linf = linf)
    }
    
    # else if unique value
  } else {
    
    # Area of the triangular distribution, the first on the left-hand side
    A <- N*(val.min-dist.min)/(2*(sort.values[2]-val.min))
    
    left <- rtrunc("rtriang",
                   n = A,
                   min = dist.min,
                   mode = val.min,
                   max = val.min,
                   lsup = lsup,
                   linf = linf)
  }
  
  
  
  ### DISTRIBUTIONS IN BETWEEN LEFT AND RIGHT TRIANGLES OF THE OVERALL SPECIES-SPECIFIC NOEC DISTRIBUTION 
  
  # Create an empty list to store the distribution parts
  mid <- list()
  
  # Creates a distribution between each pair of consecutive unique values
  for(i in 1:(length(uni.values)-1)){
    
    if(height.values[i] == 1 & height.values[i+1] == 1){
      # Calculate the uniform distributions in between
      
      mid[[i]] <- rtrunc("runif",
                         n = N,
                         min = uni.values[i],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = linf)
      
      
    } else if(isTRUE(all.equal(height.values[i], height.values[i+1]))){
      mid[[i]] <- rtrunc("runif",
                         n = N*freq.uni.values[i],
                         min = uni.values[i],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = linf)
      
    }else if(height.values[i] > height.values[i+1]){
      
      # height of shape
      h1 <- height.values[i]
      h2 <- height.values[i+1]
      # extreme edge of triangle
      max.trunc.distr <- (h1*uni.values[i+1]-h2*uni.values[i])/(h1-h2)
      # area of the triangle
      A <- ((max.trunc.distr - uni.values[i])*h1 - (max.trunc.distr - uni.values[i+1])*h2)/2
      
      mid[[i]] <- rtrunc("rtriang",
                         n = A,
                         min = uni.values[i],
                         mode = uni.values[i],
                         max = max.trunc.distr,
                         lsup = min(uni.values[i+1],lsup),
                         linf = linf)
      
      
      
    } else if(height.values[i] < height.values[i+1]){
      
      # height of shape
      h1 <- height.values[i]
      h2 <- height.values[i+1]
      # extreme edge of triangle
      min.trunc.distr <- (h2*uni.values[i] - h1*uni.values[i+1])/(h2 - h1)
      # area of the triangle
      A <- ((uni.values[i+1] - min.trunc.distr)*h2 - ( uni.values[i] - min.trunc.distr )*h1)/2
      
      mid[[i]] <- rtrunc("rtriang",
                         n = A,
                         min = min.trunc.distr,
                         mode = uni.values[i+1],
                         max = uni.values[i+1],
                         lsup = lsup,
                         linf = max(uni.values[i],linf))
      
    } 
    
    
  }
  
  
  
  ### RIGHT TRIANGLE OF THE OVERALL SPECIES-SPECIFIC NOEC PROBABILITY DISTRIBUTION
  
  if((length(unique(sort.values)) < length(sort.values)) &
     (val.max == sort.values[length(sort.values)-1])){ 
    
    # how many values are identical
    n <- length(which(sort.values == val.max))
    
    # first endpoint from the right that is not identical to the later one
    val.maxM1 <- sort.values[length(sort.values)-n]
    # Height of the last endpoint that is not similar to the next ones
    h <- height.values[length(height.values)]
    # Area of the triangular distribution, the last on the right-hand side
    A <- (dist.max - val.max)*h /2
    
    # Triangular distribution on the right-hand side
    # if the largest value is larger than the truncation limit, do nothing
    if(dist.max > lsup){
      right <- NULL
      
    } else {
      right <- rtrunc("rtriang",
                      n = A,
                      min = val.max,
                      mode = val.max,
                      max = dist.max,
                      lsup = lsup,
                      linf = linf)
    }
    
    # else if unique value
  } else {
    
    # Area of the triangular distribution, the first on the right-hand side
    A <- N*(dist.max-val.max)/(2*(val.max-sort.values[length(sort.values)-1]))
    
    right <- rtrunc("rtriang",
                    n = A, 
                    min = sort.values[length(sort.values)], 
                    mode = sort.values[length(sort.values)],
                    max = dist.max,
                    lsup = lsup,
                    linf = linf)
  }
  
  # Combine left handside, middle and right handside distributions.
  # Each distribution has the same weight (length(mid)=nb EP * N).
  step_distr <- c(left,do.call("c", mid),right)
  
  return(sample(step_distr, N))
  
}