# -------------------------------------------------------------------------------------------------
# Example of PSSD and PNEC calculations using the tool PSSD+ developed by Wigger et al. (2019)
# -------------------------------------------------------------------------------------------------
# 
# Associated publication: Systematic consideration of parameter uncertainty and variability in
#                         probabilistic species sensitivity distributions
# 
# Authors: Henning Wigger, Delphine Kawecki, Bernd Nowack and Véronique Adam
# 
# Institute: Empa, Swiss Federal Laboratories for Materials Science and Technology,
#            Technology and Society Laboratory, Lerchenfeldstrasse 5, 9014 St. Gallen, Switzerland
# 
# submitted to Integrated Environmental Assessment and Management in December 2018
# 
# -------------------------------------------------------------------------------------------------

###################################################################################################
##### MAIN PARAMETERS FOR THE CALCULATION #########################################################
###################################################################################################

# number of simulations for the triangular distributions of the data points and uncertainty factors
SIM <- 10000

# coefficient of variation for the data point distributions
CV.DP <- 0.3

# coefficient of variation for the uncertainty factor distributions
CV.UF <- 0.5

# read toxicity data
# data points from the literature
DP   <- as.matrix(read.table("data/DP_CNT.csv",   header = TRUE, sep = ","))
# uncertainty factors for conversion from acute to chronic values
UFt  <- as.matrix(read.table("data/UFt_CNT.csv",  header = TRUE, sep = ","))
# uncertainty factors for conversion of dose descriptors to NOEC
UFdd <- as.matrix(read.table("data/UFdd_CNT.csv", header = TRUE, sep = ","))


###################################################################################################
##### CALCULATION OF NOECs and PSSDs ##############################################################
###################################################################################################

# estimate the NOEC distributions for each species
NOEC <- do.pSSD(DP, UFt, UFdd, SIM, CV.DP, CV.UF)

# calculate the deterministic (modal) values of the NOEC
NOEC.det <- DP/(UFdd*UFt)


###################################################################################################
##### EXAMPLES OF RESULTS #########################################################################
###################################################################################################

### Plot the probability densities of some species

## For a species with 1 endpoint (T. thermophila)
# histogram
hist(NOEC[9,], freq = F, main = "T. thermophila - 1 endpoint", xlab = "NOEC [ug/L]")
# probability density
lines(density(NOEC[9,]), col = "red", lwd = 2) 
# Deterministic NOEC
abline(v = NOEC.det[1,9], col = "green", lty = 2, lwd = 2)

## For a species with 11 endpoints (C. dubia)
hist(NOEC[2,], freq = F, main = "C. dubia - 11 endpoints", xlab = "NOEC [ug/l]")
lines(density(NOEC[2,]), col = "red", lwd = 2)
for(i in 1:length(which(!is.na(NOEC.det[,2])))){
  abline(v = NOEC.det[i,2], col = "green", lty = 2, lwd = 2)
}


### Plot PSSDs

# calculate all PSSD curves
xvalues <- seq(-6,7,0.001)
PSSD <- matrix(NA, SIM, length(xvalues))
# calculate the ecdf function
for(i in 1:SIM){
  the.ecdf.f <- ecdf(log(NOEC[,i], base = 10))
  PSSD[i,] <- the.ecdf.f(xvalues)
}

# Plot the mean PSSD
plot(c(-2,6), c(0,1),
     main = "CNT in the freshwater compartment",
     xlab = "NOEC [µg/l]",
     ylab = "Proportion of Species Affected",
     type = "n",
     axes = F)
axis(1, at = c(-2,0,2,4,6), labels = expression(10^-2, 10^0, 10^2, 10^4, 10^6), lwd = 0.5)
axis(2, lwd = 0.5, las = 1)
abline(h = seq(0,1,0.2), lty = 1, lwd = 0.5, col = "gray90")
abline(v = seq(-6,5,1), lty = 1, lwd = 0.5, col = "gray90")
box(lwd = 0.5)
lines(xvalues, apply(PSSD,2,function(x) {mean(x)}), col = "red", lwd = 1.5)
legend("topleft", "Mean PSSD", col = "red", lwd = 1.5)


### PNEC histogram

# calculate the PNEC
PNEC <- apply(NOEC,2,function(x){ quantile(x, probs = 0.05) })

# plot a histogram of the PNEC
hist(PNEC, xlab = "PNEC [µg/l]")
