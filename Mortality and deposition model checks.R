###############################################################################
# Prep input data
library(plotrix)
library(rstan)

wreck.raw <- read.csv(here::here("data","CAAU Count Data COASST.csv"))

# Read in formatted data
mme.in <- readRDS(here::here("data","dep_model_data_WAOR.RDS"))

###############################################################################
# Read in parameter info
param.df <- read.csv(here::here("output","MCMC p950 d80 fl5 w2 GAM.csv"))

# Function for quick evaluation of stats
quanfind <- function(X){
	X1 <- mean(X, na.rm=T)
	X2 <- median(X, na.rm=T)
	X3 <- quantile(X, 0.05)
	X4 <- quantile(X,0.95)
	out <- c(X1, X2, X3, X4)
	return(out)
}

###############################################################################

###############################################################################
# Profile Functions
# Double logistic
dlf <- function(t_d, p.start, beta.start,p.end, beta.end){
	dval <- 1/((1 + exp(-(t_d-p.start)/beta.start)) * (1 + exp(-(p.end-t_d)/beta.end)))
	return(dval)
}

slf <- function(t_d, p.start, beta.start){
	dval <- 1/(1 + exp(-(t_d-p.start)/beta.start))
	return(dval)
}

# Look at spatial and temporal profiles
# Function for calculating profiles
# It will take a specific MCMC estimate and create/add to a plot
# and if needed return the profile values
Profeval <- function(i, add=FALSE, YMAX=50, spaceortime="time", returnvals=FALSE){
	params <- param.df[i,]
	
	if(spaceortime=="time"){
	
		nday <- 212
		days <- 1:nday
		Timeprof1 <- params$M1*dlf(days, params$T1s, params$B_T1s, params$T1e, params$B_T1e)/max(dlf(days, params$T1s, params$B_T1s, params$T1e, params$B_T1e))
		Timeprof2 <- params$M2*dlf(days, params$T2s, params$B_T2s, params$T2e, params$B_T2e)/max(dlf(days, params$T2s, params$B_T2s, params$T2e, params$B_T2e))
		Ttot <- Timeprof1 + Timeprof2		
		if(add==FALSE){
			plot(Timeprof1 ~ days, col=rgb(0,0,0,0.01), type="l", lwd=1, ylab="Temporal Deposition Profile", ylim=c(0, YMAX), xlab="Day")
			lines(Timeprof2 ~ days, col=rgb(0.9,0,0,0.01), type="l", lwd=1)
		} else {
			lines(Timeprof1 ~ days, col=rgb(0,0,0,0.01), type="l", lwd=1)
			lines(Timeprof2 ~ days, col=rgb(0.9,0,0,0.01), type="l", lwd=1)
		}
		if(returnvals){
			return(Ttot)
		}
	} else {
		Bnn <- 1:1143
		spprof1 <- dlf(Bnn, params$S1s, params$B_S1s, params$S1e, params$B_S1e)/max(dlf(Bnn, params$S1s, params$B_S1s, params$S1e, params$B_S1e))
		spprof2 <- dlf(Bnn, params$S2s, params$B_S2s, params$S2e, params$B_S2e)/max(dlf(Bnn, params$S2s, params$B_S2s, params$S2e, params$B_S2e))
		spprof3 <- params$M3*slf(Bnn, params$S3s, params$B_S3s)/max(slf(Bnn, params$S3s, params$B_S3s))
		spprof2 <- spprof2 + spprof3
		if(add==FALSE){
			plot(spprof1 ~ Bnn, col=rgb(0,0,0,0.01), type="l", lwd=1, ylab="Spatial Deposition Profile", ylim=c(0, YMAX), xlab="Distance alongshore (km)")
			lines(spprof2 ~ Bnn, col=rgb(0.9,0,0,0.01), type="l", lwd=1)
		} else {
			lines(spprof1 ~ Bnn, col=rgb(0,0,0,0.01), type="l", lwd=1)
			lines(spprof2 ~ Bnn, col=rgb(0.9,0,0,0.01), type="l", lwd=1)
		}
		if(returnvals){
			return(spprof1)
		}
	}
}

# Example Time profile
Samp <- sample(1:nrow(param.df), size=1000, replace=F)
Profeval(1, add=F, YMAX=50, spaceortime="time")
for(i in Samp){
  Profeval(i, add=T, YMAX=50, spaceortime="time")
}

# Example Spatial profile
Profeval(1, add=F, YMAX=1, spaceortime="space")
for(i in Samp){
  Profeval(i, add=T, YMAX=1, spaceortime="space")
}

###############################################################################

###############################################################################
# Create heatmap
# Sample a subset of the MCMC estimates (mostly to save time)
MCMCS <- sample(1:nrow(param.df), size=1000, replace=FALSE)

# Create a storage matrix for saving results
heat.mat <- matrix(0, nrow=1143, ncol=212)

for(i in 1:1000){
  # Define time (TT) and space (DD) vectors to be 
  # evaluated across
  TT <- 1:212
  DD <- 1:1143
  
  # subset parameters
  pp <- param.df[MCMCS[i],]
  
  # Calculate the time profiles for Ev 1 and Ev 2
  T1 <- pp$M1*dlf(TT, pp$T1s, pp$B_T1s, pp$T1e, pp$B_T1e)/max(dlf(TT, pp$T1s, pp$B_T1s, pp$T1e, pp$B_T1e))
  T2 <- pp$M2*dlf(TT, pp$T2s, pp$B_T2s, pp$T2e, pp$B_T2e)/max(dlf(TT, pp$T2s, pp$B_T2s, pp$T2e, pp$B_T2e))
  
  # Calculate the space profiles for components 1, 2, 3
  S1 <- dlf(DD, pp$S1s, pp$B_S1s, pp$S1e, pp$B_S1e)/max(dlf(DD, pp$S1s, pp$B_S1s, pp$S1e, pp$B_S1e))
  S2 <- dlf(DD, pp$S2s, pp$B_S2s, pp$S2e, pp$B_S2e)/max(dlf(DD, pp$S2s, pp$B_S2s, pp$S2e, pp$B_S2e))
  S3 <- pp$M3*slf(DD, pp$S3s, pp$B_S3s)/max(slf(DD, pp$S3s, pp$B_S3s))
  
  # Calculate the mean
  # as we selected 1000 parameter values I just multiply each by 1e-3
  # and sum them (this is to save on memory)
  HM <- 1e-3*(pp$D0 + (S1 %*% t(T1)) + (S2 %*% t(T2)) + (S3 %*% t(T2)))
  heat.mat <- heat.mat + HM
}

# Define function for assigning colours based on values
col.assign <- function(X){
  XSP <- c(0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,12.5,15,17.5,20,25,30,35)
  XN <- min(which(XSP > X))
  colout <- rev(rainbow(n=20, end=0.66, alpha=1))[XN]
  return(colout)
}

# Establish plot
plot(-1000, xlim=c(0,212),ylim=c(0,1143), xaxt="n", yaxt="n", ylab="", xlab="", bty="n")

# Loop over matrix and add polygons to create heatmap (theres probably a ggplot easier)
for(i in 1:1143){
  for(j in 1:212){
   polygon(y=c(i-0.5,i-0.5,i+0.5,i+0.5), x=c(j-0.5,j+0.5,j+0.5,j-0.5), 
           col=col.assign(heat.mat[i,j]), border=NA)
  }
}

# Now just plot the contour lines at given values
# plot(-1000, xlim=c(0,212),ylim=c(0,1143), xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
contour(x=1:212, y=1:1143, z=t(heat.mat), levels=c(0.5,1,2,5,10,15,20,30), drawlabels=FALSE)

#############################################################################################
