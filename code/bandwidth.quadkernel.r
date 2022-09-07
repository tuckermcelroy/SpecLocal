bandwidth.quadkernel <- function(gamma,T,const1,const2)
{

	##########################################################################
	#
	#	bandwidth.quadkernel
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
	#	Purpose: compute optimal spectrum bandwidth for quadratic kernel
	#	Background:	
	#		The spectral density f(w) = sum_h gamma(h) e^{-iwh}
	#		can be estimated using a quadratic kernel as a
	#		taper of autocovariances.  The optimal bandwidth
	#		is obtained by the method of Politis (2003).
	#	Inputs:
	#		gamma: autocovariance of a univariate time series, 
	#			either from a process or from a sample.
	#			vector has H+1 elements, for lags 0 through H
	#		T: sample size used
	#		const1: bandwidth constant for \int x^2 Lambda(x) dx
	#		const2: bandwidth constant for \int {lambda(t)}^2 dt
	#	Outputs:
	#		band.0: bandwidth for frequency 0
	#		band.pi: bandwidth for frequency pi
	#
	############################################

	epsilon <- .03
	H <- length(gamma)-1
	gamma.alt <- gamma*(-1)^seq(0,H)
	gamma2 <- gamma*seq(0,H)^2
	gamma2.alt <- gamma2*(-1)^seq(0,H)
	spec.zero <- sum(c(rev(gamma),gamma[-1]))
	spec.pi <- sum(c(rev(gamma.alt),gamma.alt[-1]))
	spec2.zero <- -1*sum(c(rev(gamma2),gamma2[-1]))
	spec2.pi <- -1*sum(c(rev(gamma2.alt),gamma2.alt[-1]))
	# frequency zero
	little.c <- (-1/2)*spec2.zero*const1
	little.c <- max(abs(little.c),epsilon)
	big.c <- spec.zero^2*const2
	band.0 <- (big.c/(4*little.c^2*T))^(1/5)
	band.0 <- min(1/band.0,T)
	# frequency pi
	little.c <- (-1/2)*spec2.pi*const1
	little.c <- max(abs(little.c),epsilon)
	big.c <- 2*spec.pi^2*const2
	band.pi <- (big.c/(4*little.c^2*T))^(1/5)
	band.pi <- min(1/band.pi,T)

	band <- c(band.0,band.pi)
	return(band)
}
