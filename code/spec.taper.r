spec.taper <- function(data, taper, M, mesh, demean=TRUE)
{

	##########################################################################
	#
	#	spec.taper
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
	#	Purpose: compute a tapered spectral density estimate
	#	Background:	
	#		The spectral density f(w) = sum_h gamma(h) e^{-iwh},
	#		can be estimated by 
	#			\hat{f}(w) = sum_h \hat{gamma}(h) Lambda(h/bn) e^{-iwh},
	#		where \hat{gamma}(h) are the sample autocovariances,
	#		Lambda is a taper, n is sample size, and b is bandwidth fraction.	
	#	Inputs:
	#		data: univariate time series, should be stationary
	#		taper: function Lambda(x) 
	#		M: bandwidth fraction, between 1 and sample size
	#		mesh: total number of frequencies, w = 2 pi j/mesh for
	#		 [mesh/2]-mesh+1 <= j <= [mesh/2]
	#		demean: TRUE if sample mean is to be subtracted.	
	#	Outputs:
	#		\hat{f}(w) at mesh number of frequencies
	#
	#####################################################################

	n <- length(data)
	if(demean) data <- data - mean(data)
	xcov <- rep(0,n+1)
	for(j in 0:(n-1)) { xcov[j+1] <- (data[1:(n-j)] %*% data[(j+1):n])/n }
	freqs <- 2*pi*seq(floor(mesh/2)-mesh+1,floor(mesh/2))/mesh 
	fspec <- 0
	fspec <- taper(0)*xcov[1]*cos(0*freqs)
	if (M > 1)
	{
		for (h in 2:M)
		{
			fspec <- fspec + 2*taper((h-1)/M)*xcov[h]*cos((h-1)*freqs)
		}
	}		
	return(fspec)
}	