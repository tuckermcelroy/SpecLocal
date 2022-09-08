######### Script for local spectral estimation

library(VGAM)

# Alter path to local directory!
setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SpecEstLM")

## load files
source("spec.taper.r")
source("taper.flat.r")
source("parzen.taper.r")
source("parzen.window.r")
source("spec.arma.r")
source("spec.exp.r")
source("spec.local.r")
source("delta.optimal.r")
source("spec.global.r")
source("ma.bandwidth.r")
source("parzen.consts.r")
source("bandwidth.quadkernel.r")
source("ARMAauto.r")
source("polymult.r")
source("spec.localmod.r")

## load data: alter path to data directory!
gdp <- read.csv(file="GDP.csv")
temp <- read.table(file="Climate.dat")

###############################################
## PART I:  Simulations

## NOTES: generate various stationary processes.
#	Test four methods for estimating spectrum at 0 and pi:
#	(1) Parzen tapered estimate, with optimal bandwidth 
#		determined by plug-in formula from Politis (2003),
#		which relies on flat-top tapered sample autocovariances.
#		This flat-top has optimal bandwidth from the "empirical rule".
#	(2) Flat-top tapered estimate, with optimal bandwidth
#		determined by the "empirical rule"
#	(3) Local quadratic estimator, with bandwidth "delta" determined
#		by plug-in formula from the paper,
#		which relies on flat-top tapered sample autocovariances.
#		This flat-top has optimal bandwidth from the "empirical rule".
#		We also consider a range of other delta values,
#		including the optimal delta of the process
#	(4) Positive local quadratic estimator, which is as (3) but
#		is constructed so as to guarantee positivity.

# Alter path to directory desired for output files
setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SpecEstLM\\NumericalUni")

###################################################################################
### generate simulation results; case of comparing constant, quadratic, and quartic

## Instructions: 
#   To simulate either Gaussian, Lalpace, or Student t ARMA,
#     use Process 1 below.  Set j to select the jth value of phis
#     for the AR parameter, and h to select the hth value of thetas
#     for the MA parameter.  Below, within the main loop (lines 108-120)
#     comment/uncomment as appropriate for the corresponding case.
#   To simulate the polynomial Gaussian process, run the code for
#     the Process 2 below, and comment/uncomment within the main loop
#     for the fourth case (lines 108-120).

## Process 1: ARMA
phis <- c( -.9, -.5, 0, .5, .9)
thetas <- c(-.8, -.4, 0, .4, .8)
j <- 5
h <- 4
phi <- phis[j]
theta <- thetas[h]
gamma <- ARMAauto(phi,theta,1000)
spec.true <- spec.arma(phi,theta,1,2) 

## Process 2: order 2 polynomial Gaussian process 
phi1s <- c( -.9, -.5, 0, .5, .9)
phi2s <- c( -.9, -.5, 0, .5, .9)
j <- 5
h <- 4
phi1 <- phi1s[j]
phi2 <- phi2s[h]
coeff1.mat <- toeplitz(phi1^seq(0,999))
coeff1.mat[upper.tri(coeff1.mat)] <- 0
coeff2.mat <- toeplitz(phi2^seq(0,999))
coeff2.mat[upper.tri(coeff2.mat)] <- 0
spec.true <-  spec.arma(phi1,NULL,1,2) + spec.arma(phi2,NULL,2,2)

## settings
# Here set the sample size: T = 50, 100, 200, 400, 800 are used in the paper
T <- 800
Monte <- 10^4
c <- .5
parzen.const <- parzen.consts(10000)
delta.range <- seq(1,50)/200

## process calculations
gamma.mat <- t(chol(toeplitz(gamma[1:T])))
delta.opt <- delta.optimal(gamma,T)
delta.opt0 <- delta.opt[1]
delta.optpi <- delta.opt[2] 

## simulation
spec.ests0 <- NULL
spec.estspi <- NULL
for(i in 1:Monte)
{
  
  # Default case: Gaussian ARMA
  x.sim <- gamma.mat %*% rnorm(T)
  
  # Second case: Laplace ARMA
#  x.sim <- gamma.mat %*% rlaplace(T,scale=1/sqrt(2))
  
  # Third case: Student t ARMA
#  x.sim <- gamma.mat %*% rt(T,df=6)/sqrt(3/2)
  
  # Fourth case: order 2 polynomial Gaussian process
#  eps <- rnorm(1000)
#  x.sim <- coeff1.mat[(1001-T):1000,] %*% eps +
#    coeff2.mat[(1001-T):1000,] %*% (eps^2 -1)
  
	M <- min(T,max(1,ma.bandwidth(x.sim,.05,4)/c))
	gamma.hat <- acf(x.sim,type="covariance",plot=FALSE,lag.max=T-1)$acf
	gamma.hat <- gamma.hat*apply(array(seq(0,T-1)/M,c(1,1,T)),3,function(x){taper.flat(x,c)})
	delta.hat <- delta.optimal(gamma.hat,T)
	delta.hat0 <- delta.hat[1]
	delta.hatpi <- delta.hat[2] 
	parzen.band <- bandwidth.quadkernel(gamma.hat,T,parzen.const[1],parzen.const[2])
	parzen.band0 <- parzen.band[1]
	parzen.bandpi <- parzen.band[2]

	# Parzen estimate
	spec.est0 <- spec.taper(x.sim,parzen.taper,parzen.band0,T,TRUE)[T/2]
	# Flat-top estimate
	spec.est0 <- c(spec.est0,
		spec.taper(x.sim,function(x){taper.flat(x,c)},M,T,TRUE)[T/2])
	# Local quadratic at estimated delta, regress on constant
	spec.est0 <- c(spec.est0,
		spec.localmod(x.sim,delta.hat0,0,T,0,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at theoretical delta, regress on constant
	spec.est0 <- c(spec.est0,
		spec.localmod(x.sim,delta.opt0,0,T,0,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at range of delta,  regress on constant
	for(k in 1:length(delta.range))
	{
		spec.est0 <- c(spec.est0,
			spec.localmod(x.sim,delta.range[k],0,T,0,TRUE,FALSE)[[2]][T/2])
	}
	# Local positive quadratic at estimated delta,  regress on constant
	spec.est0 <- c(spec.est0,
		spec.localmod(x.sim,delta.hat0,0,T,0,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at theoretical delta,  regress on constant
	spec.est0 <- c(spec.est0,
		spec.localmod(x.sim,delta.opt0,0,T,0,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at range of delta,  regress on constant
	for(k in 1:length(delta.range))
	{
		spec.est0 <- c(spec.est0,
			spec.localmod(x.sim,delta.range[k],0,T,0,TRUE,TRUE)[[2]][T/2])
	}
	# Local quadratic at estimated delta, regress on quadratic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.hat0,0,T,1,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at theoretical delta, regress on quadratic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.opt0,0,T,1,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at range of delta,  regress on quadratic
	for(k in 1:length(delta.range))
	{
	  spec.est0 <- c(spec.est0,
	                 spec.localmod(x.sim,delta.range[k],0,T,1,TRUE,FALSE)[[2]][T/2])
	}
	# Local positive quadratic at estimated delta,  regress on quadratic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.hat0,0,T,1,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at theoretical delta,  regress on quadratic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.opt0,0,T,1,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at range of delta,  regress on quadratic
	for(k in 1:length(delta.range))
	{
	  spec.est0 <- c(spec.est0,
	                 spec.localmod(x.sim,delta.range[k],0,T,1,TRUE,TRUE)[[2]][T/2])
	}
	# Local quadratic at estimated delta, regress on quartic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.hat0,0,T,2,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at theoretical delta, regress on quartic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.opt0,0,T,2,TRUE,FALSE)[[2]][T/2])
	# Local quadratic at range of delta,  regress on quartic
	for(k in 1:length(delta.range))
	{
	  spec.est0 <- c(spec.est0,
	                 spec.localmod(x.sim,delta.range[k],0,T,2,TRUE,FALSE)[[2]][T/2])
	}
	# Local positive quadratic at estimated delta,  regress on quartic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.hat0,0,T,2,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at theoretical delta,  regress on quartic
	spec.est0 <- c(spec.est0,
	               spec.localmod(x.sim,delta.opt0,0,T,2,TRUE,TRUE)[[2]][T/2])
	# Local positive quadratic at range of delta,  regress on quartic
	for(k in 1:length(delta.range))
	{
	  spec.est0 <- c(spec.est0,
	                 spec.localmod(x.sim,delta.range[k],0,T,2,TRUE,TRUE)[[2]][T/2])
	}
	
	# Parzen estimate
	spec.estpi <- spec.taper(x.sim,parzen.taper,parzen.bandpi,T,TRUE)[T]
	# Flat-top estimate
	spec.estpi <- c(spec.estpi,
		spec.taper(x.sim,function(x){taper.flat(x,c)},M,T,TRUE)[T])
	# Local quadratic at estimated delta, regress on constant
	spec.estpi <- c(spec.estpi,
		spec.localmod(x.sim,delta.hatpi,pi,T,0,TRUE,FALSE)[[2]][T])
	# Local quadratic at theoretical delta, regress on constant
	spec.estpi <- c(spec.estpi,
		spec.localmod(x.sim,delta.optpi,pi,T,0,TRUE,FALSE)[[2]][T])
	# Local quadratic at range of delta, regress on constant
	for(k in 1:length(delta.range))
	{
		spec.estpi <- c(spec.estpi,
			spec.localmod(x.sim,delta.range[k],pi,T,0,TRUE,FALSE)[[2]][T])
	}
	# Local positive quadratic at estimated delta, regress on constant
	spec.estpi <- c(spec.estpi,
		spec.localmod(x.sim,delta.hatpi,pi,T,0,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at theoretical delta, regress on constant
	spec.estpi <- c(spec.estpi,
		spec.localmod(x.sim,delta.optpi,pi,T,0,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at range of delta, regress on constant
	for(k in 1:length(delta.range))
	{
		spec.estpi <- c(spec.estpi,
			spec.localmod(x.sim,delta.range[k],pi,T,0,TRUE,TRUE)[[2]][T])
	}
	# Local quadratic at estimated delta, regress on quadratic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.hatpi,pi,T,1,TRUE,FALSE)[[2]][T])
	# Local quadratic at theoretical delta, regress on quadratic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.optpi,pi,T,1,TRUE,FALSE)[[2]][T])
	# Local quadratic at range of delta, regress on quadratic
	for(k in 1:length(delta.range))
	{
	  spec.estpi <- c(spec.estpi,
	                  spec.localmod(x.sim,delta.range[k],pi,T,1,TRUE,FALSE)[[2]][T])
	}
	# Local positive quadratic at estimated delta, regress on quadratic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.hatpi,pi,T,1,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at theoretical delta, regress on quadratic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.optpi,pi,T,1,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at range of delta, regress on quadratic
	for(k in 1:length(delta.range))
	{
	  spec.estpi <- c(spec.estpi,
	                  spec.localmod(x.sim,delta.range[k],pi,T,1,TRUE,TRUE)[[2]][T])
	}
	# Local quadratic at estimated delta, regress on quartic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.hatpi,pi,T,2,TRUE,FALSE)[[2]][T])
	# Local quadratic at theoretical delta, regress on quartic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.optpi,pi,T,2,TRUE,FALSE)[[2]][T])
	# Local quadratic at range of delta, regress on quartic
	for(k in 1:length(delta.range))
	{
	  spec.estpi <- c(spec.estpi,
	                  spec.localmod(x.sim,delta.range[k],pi,T,2,TRUE,FALSE)[[2]][T])
	}
	# Local positive quadratic at estimated delta, regress on quartic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.hatpi,pi,T,2,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at theoretical delta, regress on quartic
	spec.estpi <- c(spec.estpi,
	                spec.localmod(x.sim,delta.optpi,pi,T,2,TRUE,TRUE)[[2]][T])
	# Local positive quadratic at range of delta, regress on quartic
	for(k in 1:length(delta.range))
	{
	  spec.estpi <- c(spec.estpi,
	                  spec.localmod(x.sim,delta.range[k],pi,T,2,TRUE,TRUE)[[2]][T])
	}
	

	spec.ests0 <- rbind(spec.ests0,spec.est0)
	spec.estspi <- rbind(spec.estspi,spec.estpi)

	if(i %% 1000 == 0) { print(i) }
}

bias.0 <- colMeans(spec.ests0 - spec.true[1])
var.0 <- sqrt(rowMeans((t(spec.ests0) - colMeans(spec.ests0))^2))
rmse.0 <- sqrt(bias.0^2 + var.0^2)
round(cbind(bias.0,var.0,rmse.0),digits=3)

bias.pi <- colMeans(spec.estspi - spec.true[2])
var.pi <- sqrt(rowMeans((t(spec.estspi) - colMeans(spec.estspi))^2))
rmse.pi <- sqrt(bias.pi^2 + var.pi^2)
round(cbind(bias.pi,var.pi,rmse.pi),digits=3)

# Write statements are commented out so files are not accidentally overwritten!
#   To write output, change output file name and uncomment

#write(t(cbind(round(cbind(bias.0,var.0,rmse.0),digits=3),
#	round(cbind(bias.pi,var.pi,rmse.pi),digits=3))),file="dgpARMA54modT800.dat",ncol=6)
 
  

###############################################
## PART II: Data

## Section 1: economic data

# Alter path as desired to directory for the figures
setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SpecEstLM\\Figures")

# loading
gdp <- ts(gdp[,2],start=1947,frequency=4)
plot(gdp,xlab="Year")
# restrict to years up through 2018, omitting 2019 and 2020
gdp.gr <- diff(log(gdp[1:288]))

# plot of growth rate; comment out next line to write the file 
#pdf(file="gdp_growth.pdf",height=8,width=10)
plot(ts(gdp.gr,start=1947,frequency=4),xlab="Year",ylab="Growth Rate")
dev.off()

# exploratory
mean(gdp.gr[1:287])
spec.ar(ts(gdp.gr[1:287],frequency=4))
acf(gdp.gr[1:287])
pacf(gdp.gr[1:287])
ar.ols(gdp.gr[1:287])
 
# plot of acf of growth rate; comment out next line to write the file 
#pdf(file="gdp_acf.pdf",height=8,width=10)
acf(gdp.gr[1:287], main = "GDP growth")
dev.off()

# restrict to last 20 years
span <- 20
period <- 4
data <- gdp.gr[(length(gdp.gr)-period*span+1):length(gdp.gr)]
T <- length(data)

# settings
c <- .5
#delta.const <- .5931
parzen.const <- parzen.consts(10000)
delta.range <- seq(1,50)/200

# acf and bandwidth estimates
M <- min(T,max(1,ma.bandwidth(data,.05,4)/c))
gamma.hat <- acf(data,type="covariance",plot=FALSE,lag.max=T-1)$acf
gamma.hat <- gamma.hat*apply(array(seq(0,T-1)/M,c(1,1,T)),3,function(x){taper.flat(x,c)})
delta.hat <- delta.optimal(gamma.hat,T)
delta.hat0 <- delta.hat[1]
#delta.hatpi <- delta.hat[2] 
parzen.band <- bandwidth.quadkernel(gamma.hat,T,parzen.const[1],parzen.const[2])
parzen.band0 <- parzen.band[1]
#parzen.bandpi <- parzen.band[2]

# Parzen estimate
spec_parzen <- spec.taper(data,parzen.taper,parzen.band0,T,TRUE)
spec.est0 <- spec_parzen[T/2]
# Flat-top estimate
spec_flat <- spec.taper(data,function(x){taper.flat(x,c)},M,T,TRUE)
spec.est0 <- c(spec.est0,spec_flat[T/2])
# Local quadratic at estimated delta
spec_local <- spec.local(data,delta.hat0,0,T,TRUE,FALSE)[[2]]
#spec_local <- spec.global(data,function(x){taper.flat(x,c)},M,delta.hat0,T,demean=TRUE,pos=FALSE,plot=FALSE)
spec.est0 <- c(spec.est0,spec_local[T/2])
# Local positive quadratic at estimated delta
spec_pos <- spec.local(data,delta.hat0,0,T,TRUE,TRUE)[[2]]
#spec_pos <- spec.global(data,function(x){taper.flat(x,c)},M,delta.hat0,T,demean=TRUE,pos=TRUE,plot=FALSE)
spec.est0 <- c(spec.est0,spec_pos[T/2])
# AR OLS competitor
spec_ar <- spec.ar(data,method="ols",n.freq=floor(T/2)+1,plot=FALSE)
spec.est0 <- c(spec.est0,spec_ar$spec[1])

# t statistic 
mu.null <- .02/period
mu.null <- .03/period
mu.null <- .04/period
tstat <- sqrt(T)*(mean(data)-mu.null)/sqrt(spec.est0)

# output results
output <- rbind(spec.est0,tstat)
rownames(output) <- c("spectrum","t-statistic")
colnames(output) <- c("Parzen","flat-top","Local","Local Pos","AR OLS")
round(output,digits=8)

## confidence intervals

# lower limit
period*(mean(data) - qnorm(.975)*sqrt(spec.est0/T))
# upper limit
period*(mean(data) - qnorm(.025)*sqrt(spec.est0/T))

# spectral plots
plot(ts(c(spec_parzen[T],spec_parzen),start=-1/2,frequency=T),ylim=c(-.00005,.00020),
     ylab="Spectral Density",xlab="Frequency/2 pi")
lines(ts(c(spec_flat[T],spec_flat),start=-1/2,frequency=T),col=2)
lines(ts(c(spec_local[T],spec_local),start=-1/2,frequency=T),col=3)
lines(ts(c(spec_pos[T],spec_pos),start=-1/2,frequency=T),col=4)
lines(ts(c(rev(spec_ar$spec),spec_ar$spec[-1]),start=-1/2,frequency=T),col=5)


## Section 2: Climate Data

# Alter path as desired to directory for the figures
setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SpecEstLM\\Figures")

# loading
temp <- ts(temp[,2],start=1880,frequency=1)
plot(temp,xlab="Year")
temp.gr <- diff(temp)

# plot of growth rate; comment out next line to write the file
#pdf(file="temp_growth.pdf",height=8,width=10)
plot(ts(temp.gr,start=1881,frequency=1),xlab="Year",ylab="Growth Rate")
dev.off()

# exploratory
mean(temp.gr)
spec.ar(ts(temp.gr,frequency=1))
acf(temp.gr)
pacf(temp.gr)
ar.ols(temp.gr)

# plot of acf of growth rate; comment out next line to write the file
#pdf(file="temp_acf.pdf",height=8,width=10)
acf(temp.gr, main = "GLOTI growth")
dev.off()

# global settings
c <- .5
#delta.const <- .5931
parzen.const <- parzen.consts(10000)
delta.range <- seq(1,50)/200

# split the data into two spans
period <- 1
cut <- 70
sample1 <- temp.gr[1:(cut-1)]
sample2 <- temp.gr[(cut+1):140]

# compute f(0) 
data <- c(temp.gr[1:cut]-mean(sample1),
          temp.gr[(cut+1):140]-mean(sample2))
T <- length(data)

# acf and bandwidth estimates
#  Note: code is written to subtract out sample mean in acf calculations.
#   Here we don't want to do this, since we have already subtracted out
#   means for the two subsamples. 
#M <- min(T,max(1,ma.bandwidth(data,.05,4,demean=FALSE)/c))
M <- 4
gamma.hat <- acf(data,type="covariance",plot=FALSE,lag.max=T-1,demean=FALSE)$acf
gamma.hat <- gamma.hat*apply(array(seq(0,T-1)/M,c(1,1,T)),3,function(x){taper.flat(x,c)})
delta.hat <- delta.optimal(gamma.hat,T)
delta.hat0 <- delta.hat[1]
#delta.hatpi <- delta.hat[2] 
parzen.band <- bandwidth.quadkernel(gamma.hat,T,parzen.const[1],parzen.const[2])
parzen.band0 <- parzen.band[1]
#parzen.bandpi <- parzen.band[2]

# Parzen estimate
spec_parzen <- spec.taper(data,parzen.taper,parzen.band0,T,FALSE)
spec.est0 <- spec_parzen[T/2]
# Flat-top estimate
spec_flat <- spec.taper(data,function(x){taper.flat(x,c)},M,T,FALSE)
spec.est0 <- c(spec.est0,spec_flat[T/2])
# Local quadratic at estimated delta
spec_local <- spec.local(data,delta.hat0,0,T,FALSE,FALSE)[[2]]
#spec_local <- spec.global(data,function(x){taper.flat(x,c)},M,delta.hat0,T,demean=FALSE,pos=FALSE,plot=FALSE)
spec.est0 <- c(spec.est0,spec_local[T/2])
# Local positive quadratic at estimated delta
spec_pos <- spec.local(data,delta.hat0,0,T,FALSE,TRUE)[[2]]
#spec_pos <- spec.global(data,function(x){taper.flat(x,c)},M,delta.hat0,T,demean=FALSE,pos=TRUE,plot=FALSE)
spec.est0 <- c(spec.est0,spec_pos[T/2])
# AR OLS competitor
spec_ar <- spec.ar(data,method="ols",n.freq=floor(T/2)+1,plot=FALSE,demean=FALSE)
spec.est0 <- c(spec.est0,spec_ar$spec[1])

# t statistic 
tstat <- (mean(sample1)-mean(sample2))/
  sqrt(spec.est0 * ( 1/length(sample1) + 1/length(sample2) ))

# output results
output <- rbind(spec.est0,tstat)
rownames(output) <- c("spectrum","t-statistic")
colnames(output) <- c("Parzen","flat-top","Local Quad","Log-periodogram","AR OLS")
round(output,digits=8)

# spectral plots
plot(ts(c(spec_parzen[T],spec_parzen),start=-1/2,frequency=T),
     ylab="Spectral Density",xlab="Frequency/2 pi")
lines(ts(c(spec_flat[T],spec_flat),start=-1/2,frequency=T),col=2)
lines(ts(c(spec_local[T],spec_local),start=-1/2,frequency=T),col=3)
lines(ts(c(spec_pos[T],spec_pos),start=-1/2,frequency=T),col=4)
lines(ts(c(rev(spec_ar$spec),spec_ar$spec[-1]),start=-1/2,frequency=T),col=5)

