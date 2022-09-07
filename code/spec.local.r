spec.local <- function(data,delta,omega,mesh,demean=TRUE,pos=TRUE)
{

	##########################################################################
	#
	#	spec.local
	# 	    Copyright (C) 2020  Tucker McElroy
	#
	#    This program is free software: you can redistribute it and/or modify
	#    it under the terms of the GNU General Public License as published by
	#    the Free Software Foundation, either version 3 of the License, or
	#    (at your option) any later version.
	#
	#    This program is distributed in the hope that it will be useful,
	#    but WITHOUT ANY WARRANTY; without even the implied warranty of
	#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#    GNU General Public License for more details.
	#
	#    You should have received a copy of the GNU General Public License
	#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	#
	############################################################################

	################# Documentation #####################################
	#
	#	Purpose: compute a spectral density estimate by local linear regression
	#	Background:	
	#		The spectral density f(w) = sum_h gamma(h) e^{-iwh},
	#		can be estimated by 
	#			\hat{f}(w) = solve(t(X) %* X) %*% t(X) %*% Y,
	#		where X has a column of ones and a column of (w-wj)^2,
	#		and wj are frequencies within 2*pi*delta of w,
	#		and Y is column of periodogram at wj.
	#	Inputs:
	#		data: univariate time series, should be stationary
 	#		delta : neighborhood (omega-2pidelta,omega+2pidelta) 
	#			of frequencies are considered for the regression;
	#			0 < delta <= .5
	#		omega : frequency of interest, between -pi and pi
	#		mesh: total number of frequencies, w = 2 pi j/mesh for
	#		 [mesh/2]-mesh+1 <= j <= [mesh/2]
	#		demean: TRUE if sample mean is to be subtracted.	
	#		pos: TRUE if a positive estimate is desired.
	#	Outputs:
	#		\hat{f}(w) at mesh number of frequencies
	#
	############################################

	euler <- 0.5772156649
	n <- length(data)
	if(demean) data <- data - mean(data)
	xcov <- rep(0,n+1)
	for(j in 0:(n-1)) { xcov[j+1] <- (data[1:(n-j)] %*% data[(j+1):n])/n }
	freqs <- 2*pi*seq(floor(mesh/2)-mesh+1,floor(mesh/2))/mesh 
	perg <- cos(0*freqs)*xcov[1]
	for(h in 1:(n-1))
	{ perg <- perg + 2*cos(h*freqs)*xcov[h+1] }

	# determine frequency band
	mid <- n - floor(n/2)
	m <- ceiling(delta*n)
	if(omega == 0)
	{
		low <- mid + 1
		hi <- mid + m
	}
	if((omega == pi) || (omega == -1*pi))
	{
		hi <- n
		low <- n - m + 1
	}
	freqs.sub <- freqs[low:hi]
 
	# compute regressors	
	reg <- (abs(omega) - freqs.sub)^2
	y <- perg[low:hi]
	if(pos) { y <- log(perg[low:hi]) + euler }
	if(m > 1) {
	 	lm.fit <- lm(y ~ reg)
		coefs <- lm.fit$coef[1:2]
	} else { coefs <- c(y,0) }
	
	# get spectral estimate
	fspec <- coefs[1] + coefs[2]*(omega-freqs)^2
	if(pos) { fspec <- exp(fspec) }

	return(list(coefs,fspec))
}

