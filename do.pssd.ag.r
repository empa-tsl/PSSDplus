# Function for generating NOEC distributions for each species
# --- EXAMPLE OF A SPECIAL CASE FOR SILVER
# 
# Yields a matrix with dimensions number of species * number of iterations in the simulation
# 
# Arguments:   - DP : matrix of data points
#              - UFt : matrix of uncertainty factors for the exposure time
#              - UFdd : matrix of uncertainty factors for the dose-descriptor
#              - SIM : number of iterations in the simulation
#              - CV.DP : coefficient of variation for the interlaboratory variation
#              - CV.UF : coefficient of variation for the use of non-substance-specific
#                uncertainty factors 
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



do.pSSD.Ag <- function(DP,
                       UFt,
                       UFdd,
                       SIM,
                       CV.DP,
                       CV.UF){
  
  # test if there is no data available for one species
  if(any(apply(DP,2,function(x) length(which(!is.na(x)))) == 0)){
    warning("No data is available for one or more species, it/they won't contribute to the PSSD calculation.")
    # find which species has no data
    ind.sp.rem <- which(apply(DP,2,function(x) length(which(!is.na(x)))) == 0)
    # remove those columns
    DP <- DP[,-ind.sp.rem]
    UFt <- UFt[,-ind.sp.rem]
    UFdd <- UFdd[,-ind.sp.rem]
  }
  
  # Create the step distributions (or triangular or trapezoidal) for each species
  # Create an empty matrix in which step distributions will be compiled
  NOEC_comb <- matrix(NA, ncol(DP), SIM,
                      dimnames = list(colnames(DP), NULL))
  
  # Fill in the matrix. If there is only one data point, NOEC stays the same. If  there are
  # 2 endpoints, a uniform distribution is produced. If  there are more than 2 endpoints, a step
  # distribution is produced. One line is for one species.
  require(trapezoid)
  require(mc2d)
  
  # store the corrected endpoints
  corr.endpoints <- DP/(UFdd*UFt)
  sort.endpoints <- apply(corr.endpoints,2, sort)
  
  for (sp in colnames(DP)){
    
    # store the indices of the minimal and maximal data point
    ind.min <- which.min(corr.endpoints[,sp])
    ind.max <- which.max(corr.endpoints[,sp])
    
    # calculate the theoretical minimum and maximum of the distribution we are looking for
    if (UFt[ind.min, sp] == 10){
      sp.min <- corr.endpoints[ind.min,sp]*(1-(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2)))/100
    } else {
      sp.min <- corr.endpoints[ind.min,sp]*(1-(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    }
    if (UFt[ind.max, sp] == 10){
      sp.max <- corr.endpoints[ind.max,sp]*(1+(sqrt((CV.DP/2.45)^2) + (CV.UF/2.45)^2))/1
    } else {
      sp.max <- corr.endpoints[ind.max,sp]*(1+(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    }
    
    # For species with one unique data point, NOEC stays the same:
    if(length(unique(sort.endpoints[[sp]])) == 1){
      NOEC_comb[sp,] <- rtrunc("rtriang", min = sp.min, 
                               mode = sort.endpoints[[sp]][1],
                               max = sp.max,
                               n = SIM, linf = 0)
      
      
      # For species with two endpoints:
    } else if(length(sort.endpoints[[sp]]) == 2){
      # Create a trapezoidal distribution including both endpoints
      NOEC_comb[sp,] <- rtrunc("rtrapezoid", SIM,
                               mode1 = sort.endpoints[[sp]][1],
                               mode2 = sort.endpoints[[sp]][2],
                               min = sp.min, max = sp.max,
                               linf = 0)
      
      
      # For species with three endpoints or more:
    } else {
      
      
      # Sample from this step distribution for each species
      NOEC_comb[sp,] <- rmore(values = sort.endpoints[[sp]], max = sp.max, min = sp.min, N = SIM, linf = 0)
      
    } 
  }
  # return the whole matrix
  return(NOEC_comb)
  
}

