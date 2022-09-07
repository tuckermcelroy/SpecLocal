parzen.consts <- function(mesh)
{

	##########################################################################
	#
	#	parzen.consts
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
	#	Purpose: compute constants for Parzen kernel
	#	Background:	
	#		Parzen taper lambda and Parzen window Lambda;
	#		computes constants \int x^2 Lambda(x) dx and
	#		  \int {lambda(t)}^2 dt
	#	Inputs:
	#		mesh: Riemann mesh
	#	Outputs:
	#		const1: bandwidth constant for \int x^2 Lambda(x) dx
	#		const2: bandwidth constant for \int {lambda(t)}^2 dt
	#	Requires: parzen.taper.r and parzen.window.r
	#
	############################################

	extent <- 100
	y <- apply(array(seq(1,extent*mesh)/mesh,c(1,1,extent*mesh)),3,parzen.window)
	y <- y*(seq(1,extent*mesh)/mesh)^2
	const1 <- 2*extent*mean(y)

	y <- apply(array(seq(1,mesh)/mesh,c(1,1,mesh)),3,parzen.taper)
	const2 <- 2*mean(y^2)

	return(c(const1,const2))
}

