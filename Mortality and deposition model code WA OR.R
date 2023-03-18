#############################################################################################################
# Washington and Oregon Model
# Mortality and deposition
library(rstan)
library(plotrix)
library(here)
library(shinystan)

############################################################################################################

############################################################################################################
# Read in raw data
wreck.raw <- read.csv(here::here("data","CAAU Count Data COASST.csv"))
Nbeach <- dim(wreck.raw)[1] 

# ********* CHANGE VALUES HERE *********** #
# Define variable model parameters
psi <- rep(0.95,Nbeach)		# Persistence
p.val <- 0.80					      # Detection

# deposition matrix
DepM <- read.csv(here::here("data","DM GAM fl5 w2.csv"))

# Write out model formulation
model.fullname <- "p950 d80 fl5 w2 GAM"
# **************************************** #

############################################################################################################

############################################################################################################
# Read in data
mme.in <- readRDS(here::here("data","dep_model_data_WAOR.RDS"))

###################################################################################################################################################################

###############################################################################################################
# Model fitting
nchain = 4
nwarm = 100
niter = 200

mme.test <- stan(file="mortality_and_deposition_model_WAOR.stan",
                 data = mme.in,
                 chains = nchain, 
                 warmup = nwarm, 
                 iter = niter,
                 cores = 4, verbose = FALSE, 
                 control=list(adapt_delta=0.999), refresh = 20)

s1 <- summary(mme.test)$summary
# write.csv(s1, file=paste("Summary ", model.fullname, ".csv", sep=""), row.names=F) 

###############################################################################################################
###############################################################################################################

###############################################################################################################
# Convergence check
s1
traceplot(mme.test, pars=c("D0","nu"))
traceplot(mme.test, pars=c("T1s","T1e","B_T1e","B_T1s","M1"))                        
traceplot(mme.test, pars=c("T2s","T2e","B_T2e","B_T2s","M2"))
traceplot(mme.test, pars=c("S1s","S1e","B_S1e","B_S1s","S3s","B_S3s"))
traceplot(mme.test, pars=c("M3", "S2s", "S2e", "B_S2e", "B_S2s"))

launch_shinystan(mme.test)

################################################################################

################################################################################
################################################################################
# Model evaluation
param.list <- c("D0","T1s","T1e","B_T1e","B_T1s","M1","T2s","T2e","B_T2e","B_T2s",
                "M2","S1s","S1e","B_S1e","B_S1s", "S3s", "B_S3s", "M3", "S2s", "S2e", "B_S2e", "B_S2s","nu")
mme.coda <- extract(mme.test, pars=param.list, permuted=FALSE, inc_warmup=FALSE) 

# Now format
npar <- dim(mme.coda)[3]
nchain <- nchain
nperm <- niter-nwarm

param.df <- matrix(0, nrow=nchain*nperm, ncol=npar)

for(i in 1:npar){
	tmp1 <- mme.coda[,,i]
	parfill <- tmp1[,1]
	for(j in 2:nchain){
		parfill <- c(parfill, tmp1[,j])
	}
	param.df[,i] <- parfill
}

param.df <- data.frame(param.df)
names(param.df) <- param.list

# write.csv(param.df, file=paste("MCMC ", model.fullname, ".csv", sep=""), row.names=FALSE)

# Histograms of parameter inputs
plotblank <- function(){
	plot(-100,ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
}

par(mfcol=c(5,5), mar=c(3,3,1,1))

hist(param.df$T1s, 50, freq=FALSE, xlab="", main="T1s")
hist(param.df$T1e, 50, freq=FALSE, xlab="", main="T1e")
hist(param.df$B_T1s, 50, freq=FALSE, xlab="", main="B_T1s")
hist(param.df$B_T1e, 50, freq=FALSE, xlab="", main="B_T1e")
hist(param.df$M1, 50, freq=FALSE, xlab="", main="M1")

hist(param.df$S1s, 50, freq=FALSE, xlab="", main="S1s")
hist(param.df$S1e, 50, freq=FALSE, xlab="", main="S1e")
hist(param.df$B_S1s, 50, freq=FALSE, xlab="", main="B_S1s")
hist(param.df$B_S1e, 50, freq=FALSE, xlab="", main="B_S1e")
plotblank()

hist(param.df$T2s, 50, freq=FALSE, xlab="", main="T2s")
hist(param.df$T2e, 50, freq=FALSE, xlab="", main="T2e")
hist(param.df$B_T2s, 50, freq=FALSE, xlab="", main="B_T2s")
hist(param.df$B_T2e, 50, freq=FALSE, xlab="", main="B_T2e")
hist(param.df$M2, 50, freq=FALSE, xlab="", main="M2")

hist(param.df$S2s, 50, freq=FALSE, xlab="", main="S2s")
hist(param.df$S2e, 50, freq=FALSE, xlab="", main="S2e")
hist(param.df$B_S2s, 50, freq=FALSE, xlab="", main="B_S2s")
hist(param.df$B_S2e, 50, freq=FALSE, xlab="", main="B_S2e")
plotblank()

hist(param.df$S3s, 50, freq=FALSE, xlab="", main="S3s")
hist(param.df$B_S3s, 50, freq=FALSE, xlab="", main="B_S3s")
hist(param.df$M3, 50, freq=FALSE, xlab="", main="M3")
hist(param.df$D0, 50, freq=FALSE, xlab="", main="D0")
hist(param.df$nu, 50, freq=FALSE, xlab="", main="nu")

################################################################################

################################################################################
# Estimating deposition

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

#
seg <- read.csv(here::here("data",
                           "Beach Segment file for total deposition calculation.csv"))
seg$distvec <- (seg$CumuStart + seg$CumuEnd)/2

# Deposition function
TotDep <- function(i, segin){
	
	# define parameters to use
	params <- param.df[i,]
	dayvec <- 1:212
	distvec <- 1:1143
	BigDMat <- mme.in$BigDMat
	BigDMat <- BigDMat[1:212,]
	
	# Scaling
	TPmag1 <- dlf(dayvec, params$T1s, params$B_T1s, params$T1e, params$B_T1e)
	TPmag2 <- dlf(dayvec, params$T2s, params$B_T2s, params$T2e, params$B_T2e)
	SPmag1 <- dlf(distvec, params$S1s, params$B_S1s, params$S1e, params$B_S1e)
	SPmag2 <- dlf(distvec, params$S2s, params$B_S2s, params$S2e, params$B_S2e)
	SPmag3 <- slf(distvec, params$S3s, params$B_S3s)
	
	TPscale1 <- max(TPmag1)
	TPscale2 <- max(TPmag2)
	SPscale1 <- max(SPmag1)
	SPscale2 <- max(SPmag2)
	SPscale3 <- max(SPmag3)
	
	# Generate time signal
	Timeprof1 <- params$M1*dlf(dayvec, params$T1s, params$B_T1s, params$T1e, params$B_T1e)/TPscale1
	Timeprof2 <- params$M2*dlf(dayvec, params$T2s, params$B_T2s, params$T2e, params$B_T2e)/TPscale2
		
	# Now convert to deposition
	MortMat1 <- matrix(0, nrow=212, ncol=212) 
	MortMat2 <- matrix(0, nrow=212, ncol=212) 
	for(j in 1:212){
		MortMat1[,j] <- Timeprof1[j]*BigDMat[,j]
		MortMat2[,j] <- Timeprof2[j]*BigDMat[,j]
	}
	
	Temp1 <- rep(0, 212)
	Temp2 <- rep(0, 212)
	for(j in 1:212){
		Temp1[j] = sum(MortMat1[j,])
		Temp2[j] = sum(MortMat2[j,])
	}
	
	Spaceprof1 <- dlf(seg$distvec,params$S1s, params$B_S1s, params$S1e, params$B_S1e)/SPscale1
	Spaceprof2 <- dlf(seg$distvec,params$S2s, params$B_S2s, params$S2e, params$B_S2e)/SPscale2
	Spaceprof3 <- params$M3*slf(seg$distvec,params$S3s, params$B_S3s)/SPscale3
	
	Abval <- matrix(0, nrow=length(seg$distvec), ncol=length(dayvec))
	for(j in 1:length(seg$distvec)){
		Abval[j,] <- (params$D0 + (Temp1*Spaceprof1[j] + Temp2*(Spaceprof2[j] + Spaceprof3[j])))*segin$SegSize[j]*segin$SegProp[j]
	}
	
	TD <- sum(Abval)

	return(TD)
}

# repeat
nn <- nrow(param.df)
Depvec <- rep(0,nn)
for(i in 1:nn){
	Depvec[i] <- TotDep(i, segin=seg)
}

par(mfrow=c(1,1))
hist(Depvec, 70, main="Total Deposition", xlab="Deposition", ylab="Frequency")

mean(Depvec)
quantile(Depvec, c(0.5,0.025,0.05,0.25,0.75,0.95,0.975))

###############################################################################

###############################################################################
# Mortality function
TotMort <- function(i, segin){
	
	# define parameters to use
	params <- param.df[i,]
	dayvec <- 1:212
	distvec <- 1:1143
	
	# Scaling
	TPmag1 <- dlf(dayvec, params$T1s, params$B_T1s, params$T1e, params$B_T1e)
	TPmag2 <- dlf(dayvec, params$T2s, params$B_T2s, params$T2e, params$B_T2e)
	SPmag1 <- dlf(distvec, params$S1s, params$B_S1s, params$S1e, params$B_S1e)
	SPmag2 <- dlf(distvec, params$S2s, params$B_S2s, params$S2e, params$B_S2e)
	SPmag3 <- slf(distvec, params$S3s, params$B_S3s)
	
	TPscale1 <- max(TPmag1)
	TPscale2 <- max(TPmag2)
	SPscale1 <- max(SPmag1)
	SPscale2 <- max(SPmag2)
	SPscale3 <- max(SPmag3)
	
	# Generate time signal
	Temp1 <- params$M1*dlf(dayvec, params$T1s, params$B_T1s, params$T1e, params$B_T1e)/TPscale1
	Temp2 <- params$M2*dlf(dayvec, params$T2s, params$B_T2s, params$T2e, params$B_T2e)/TPscale2
		
	Spaceprof1 <- dlf(seg$distvec,params$S1s, params$B_S1s, params$S1e, params$B_S1e)/SPscale1
	Spaceprof2 <- dlf(seg$distvec,params$S2s, params$B_S2s, params$S2e, params$B_S2e)/SPscale2
	Spaceprof3 <- params$M3*slf(seg$distvec,params$S3s, params$B_S3s)/SPscale3
	
	Abval <- matrix(0, nrow=length(seg$distvec), ncol=length(dayvec))
	for(j in 1:length(seg$distvec)){
		Abval[j,] <- ((Temp1*Spaceprof1[j] + Temp2*(Spaceprof2[j] + Spaceprof3[j])))*segin$SegSize[j]*segin$SegProp[j]
	}
	
	TD <- sum(Abval)

	return(TD)
}

Mortvec <- rep(0,nn)
for(i in 1:nn){
	Mortvec[i] <- TotMort(i, segin=seg)
}

hist(Mortvec, 70, main="Total Mortality", xlab="Mortality", ylab="Frequency")

mean(Mortvec)
quantile(Mortvec, c(0.5,0.025,0.05,0.25,0.75,0.95,0.975))

###############################################################################






