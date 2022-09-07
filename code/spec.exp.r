spec.exp <- function(tau, innovar, mesh)
{

	##########################################################################
	#
	#	spec.exp
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
	#	Purpose: compute the EXP spectral density 
	#	Background:	
	#		The EXP spectral density f(w) = sum_h gamma(h) e^{-iwh}
	#		is given by 
	#		  f(w) = exp { tau(e^{-iw}) + tau(e^{iw}) } sigma^2,
	#		where tau(z) = tau[1] z ... + tau[m] z^m,
	#		sigma is the innovation standard deviation.
	#	Inputs:
	#		tau: vector of m elements, cepstral coefficients (can be NULL)
	#		innovar: sigma^2, the innovation variance
	#		mesh: total number of frequencies, w = 2 pi j/mesh for
	#		 [mesh/2]-mesh+1 <= j <= [mesh/2]
	#	Outputs:
	#		f(w) at mesh number of frequencies
	#
	#####################################################################

	m <- length(tau)
	freqs <- 2*pi*seq(floor(mesh/2)-mesh+1,floor(mesh/2))/mesh 
	fspec <- 0*freqs
	if(m > 0) 
	{ for(j in 1:m) { fspec <- fspec + tau[j]*2*cos(freqs*j) } }
	fspec <- innovar*exp(fspec)	
	return(fspec)
}