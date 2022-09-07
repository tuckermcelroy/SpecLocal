spec.arma <- function(phi, theta, innovar, mesh)
{

	##########################################################################
	#
	#	spec.arma
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
	#	Purpose: compute the ARMA spectral density 
	#	Background:	
	#		The ARMA spectral density f(w) = sum_h gamma(h) e^{-iwh}
	#		is given by 
	#				theta(e^{-iw}) theta(e^{iw}) sigma^2
	#		  f(w) =    ----------------------------
	#				  phi(e^{-iw}) phi(e^{iw})
	#		where phi(z) = 1 - phi[1] z ... - phi[p] z^p
	#		and theta(z) = 1 + theta[1] z ... + theta[q] z^q,
	#		sigma is the innovation standard deviation.
	#	Inputs:
	#		phi: vector of p elements, AR coefficients (can be NULL)
	#		theta: vector of q elements, MA coefficients (can be NULL)
	#		innovar: sigma^2, the innovation variance
	#		mesh: total number of frequencies, w = 2 pi j/mesh for
	#		 [mesh/2]-mesh+1 <= j <= [mesh/2]
	#	Outputs:
	#		f(w) at mesh number of frequencies
	#
	#####################################################################

	p <- length(phi)
	q <- length(theta)
	freqs <- 2*pi*seq(floor(mesh/2)-mesh+1,floor(mesh/2))/mesh 
	fnum <- exp(-1i*0*freqs)
	if(q > 0) 
	{ for(j in 1:q) { fnum <- fnum + theta[j]*exp(-1i*j*freqs) } }
	fdenom <- exp(-1i*0*freqs)
	if(p > 0) 
	{ for(j in 1:p) { fdenom <- fdenom - phi[j]*exp(-1i*j*freqs) } }
	fspec <- innovar*Mod(fnum)^2/Mod(fdenom)^2	
	return(fspec)
}