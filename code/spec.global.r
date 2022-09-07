spec.global <- function(data,taper,M,delta,mesh,demean=TRUE,pos=TRUE,plot=TRUE)
{

	##########################################################################
	#
	#	spec.global
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
	#	Purpose: compute a spectral density estimate by combining
	#		tapered spectral estimation and local linear regression
	#	Background:	
	#		see spec.taper and spec.local
	#	Inputs:
	#		data: univariate time series, should be stationary
	#		taper: function Lambda(x) 
	#		M: bandwidth fraction, between 1 and sample size
 	#		delta : neighborhood (omega-2pidelta,omega+2pidelta) 
	#			of frequencies are considered for the regression;
	#			0 < delta <= .5
	#		omega : frequency of interest, between -pi and pi
	#		mesh: total number of frequencies, w = 2 pi j/mesh for
	#		 [mesh/2]-mesh+1 <= j <= [mesh/2]
	#		demean: TRUE if sample mean is to be subtracted.	
	#		pos: TRUE if a positive estimate is desired.
	#		plot: TRUE if a plot is desired
	#	Outputs:
	#		\hat{f}(w) at mesh number of frequencies
	#	Requires: spec.taper, spec.local
	#
	############################################

	kappa <- function(w,delta,omega)
	{
		val <- pmax(2*pi*delta - abs(w-omega),0)
		return(val/(2*pi*delta))
	}

	freqs <- 2*pi*seq(floor(mesh/2)-mesh+1,floor(mesh/2))/mesh 
	spec.middle <- spec.taper(data,taper,M,mesh,demean)
	spec.zero <- spec.local(data,delta,0,mesh,demean,pos)[[2]]
	spec.pi <- spec.local(data,delta,pi,mesh,demean,pos)[[2]]
	spec.mpi <- spec.local(data,delta,-pi,mesh,demean,pos)[[2]]
 	spec.est <- (1 - kappa(freqs,delta,0) - kappa(freqs,delta,pi) - 
		kappa(freqs,delta,-pi))*pmax(spec.middle,0) + 
		kappa(freqs,delta,0)*spec.zero + kappa(freqs,delta,pi)*pmax(spec.pi,0) + 
		kappa(freqs,delta,-pi)*pmax(spec.mpi,0)

	if(demean) { data <- data - mean(data) }
	v2 <- mean(data^2)
	normalize <- mean(spec.est)/v2
	spec.est <- spec.est/normalize

	if(plot)
	{
		plot(ts(spec.est,frequency=floor(mesh/2),start=-1),
			xlab="Cycles",ylab="Spectrum")
		lines(ts(spec.middle,frequency=floor(mesh/2),start=-1),col=2)
	}

	return(spec.est)
}

