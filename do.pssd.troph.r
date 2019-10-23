# Function for generating NOEC distributions for each species taking into account trophic levels
# 
# Yields a matrix with dimensions number of species * number of iterations in the simulation
# 
# Arguments:   - DP : matrix of endpoints
#              - UFt : matrix of uncertainty factors for the exposure time
#              - UFdd : matrix of uncertainty factors for the dose-descriptor
#              - SIM : number of iterations in the simulation
#              - CV.DP : coefficient of variation for the interlaboratory variation
#              - CV.UF : coefficient of variation for the use of non-substance-specific
#                uncertainty factors 
#              - TrophLevel : Trophic levels of each species. 1 for primary producers, 2 for herbivores and 3 for carnivores
#              - Troph.Perc.Typ : Typical fractions in which each trophic level is represented in freshwater communities. 
# 
# Date of last modification: 03.07.2019
# 
# Associated publication: Systematic consideration of parameter uncertainty and variability in
#                         probabilistic species sensitivity distributions
# 
# Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and Veronique Adam
# 
# Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology,
#            Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
# 
# submitted to Integrated Environmental Assessment and Management in July 2019

# -------------------------------------------------------------------------------------------------



do.pSSD.troph <- function(DP,
                          UFt,
                          UFdd,
                          SIM,
                          CV.DP,
                          CV.UF,
                          TrophLevel,
                          Troph.Perc.Typ = c(0.65, 0.25, 0.1)){
  
  # test if there is no data available for one species
  if(any(apply(DP,2,function(x) length(which(!is.na(x)))) == 0)){
    warning("No data is available for one or more species, it/they won't contribute to the PSSD calculation.")
    # find which species has no data
    ind.sp.rem <- which(apply(DP, 2, function(x) length(which(!is.na(x)))) == 0)
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
  sort.endpoints <- apply(corr.endpoints, 2, sort)
  
  for (sp in colnames(DP)){
    
    # store the indices of the minimal and maximal data point
    ind.min <- which.min(corr.endpoints[,sp])
    ind.max <- which.max(corr.endpoints[,sp])
    
    # calculate the theoretical minimum and maximum of the distribution we are looking for
    sp.min <- corr.endpoints[ind.min,sp]*(1-(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    sp.max <- corr.endpoints[ind.max,sp]*(1+(sqrt((CV.DP/2.45)^2 + (CV.UF/2.45)^2 + (CV.UF/2.45)^2)*2.45))
    
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
  
  # Sample species with typical proportions representative of a typical community (if terrestrial or marine, Troph.Perc need to be changed)
  Nb_Troph <- c(length(which(TrophLevel == 1)),
                length(which(TrophLevel == 2)),
                length(which(TrophLevel == 3)))
  # calculate the shares of the different trophic levels in the data
  Troph.Perc.Data <- Nb_Troph/ncol(DP)
  
  
  # identify the trophic level that will be unchanged
  if(length(which(Troph.Perc.Data < Troph.Perc.Typ)) == 1){
    lev.fix <- which(Troph.Perc.Data < Troph.Perc.Typ)
  } else {
    lev.fix <- which.max((Troph.Perc.Typ - Troph.Perc.Data)/Troph.Perc.Typ)
  }
  
  # create an empty vector to contain the corrected number of species for the species being corrected
  Nb_points <- rep(NA,3)
  # create an empty list with an element per trophic level
  NOEC_lev <- vector("list", 3)
  
  for(lev in 1:3){
    
    # for fix level, take as such
    if(lev == lev.fix){
      Nb_points[lev] <- Nb_Troph[lev.fix]
      
      # identify species corresponding to level lev
      ind.sp <- which(TrophLevel == lev) # vector of all indices for which the species level is the fixed one
      # store in list
      NOEC_lev[[lev]] <- NOEC_comb[ind.sp,]
      
    } else {
      
      # change the remaining levels
      Nb_points[lev] <- round(Nb_Troph[lev.fix] * Troph.Perc.Typ[lev] / Troph.Perc.Data[lev.fix]) # Nouveau nombre de lignes
      
      ind.sp <- which(TrophLevel == lev) # vector of all indices of the level
      # store in temporary matrix the NOEC distributions of all species of the level
      temp <- NOEC_comb[ind.sp,] 
      
      NOEC_lev[[lev]] <- matrix(NA, Nb_points[lev],SIM)
      for(i in 1:SIM){
        NOEC_lev[[lev]][,i] <- sample(temp[,i], Nb_points[lev])
      }
    }
    
  }
  
  # combine that stuff
  NOEC_TrophCorr <- do.call(rbind,NOEC_lev)
  
  # return the whole matrix
  return(NOEC_TrophCorr)
  
}
