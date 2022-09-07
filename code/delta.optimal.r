delta.optimal <- function(gamma,TT)
{
  
  ##########################################################################
  #
  #	delta.optimal
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
  #	Purpose: compute optimal spectrum bandwidth for local linear regression
  #	Background:	
  #		The spectral density f(w) = sum_h gamma(h) e^{-iwh},
  #		can be estimated by 
  #			\hat{f}(w) = solve(t(X) %* X) %*% t(X) %*% Y,
  #		where X has a column of ones and a column of (w-wj)^2,
  #		and wj are frequencies within 2*pi*delta of w,
  #		and Y is column of periodogram at wj.
  #		delta is a neighborhood (omega-2pidelta,omega+2pidelta),
  #		  of frequencies that are considered for the regression;
  #	        omega is 0 or pi, and	0 < delta <= .5
  #	Inputs:
  #		gamma: autocovariance of a univariate time series, 
  #			either from a process or from a sample.
  #		TT: sample size being used, should be <= length(gamma)
  #	Outputs:
  #		delta.0: for frequency zero
  #		delta.pi: for frequency pi
  #
  ############################################
  
  mse <- function(delta,gamma,omega,TT)
  {
    n <- length(gamma)
    m <- ceiling(delta*TT)
    freqs <- 2*pi*seq(floor(TT/2)-TT+1,floor(TT/2))/TT
    fspec <- gamma[1]*cos(0*freqs)
    if(n > 1)
    {
      for(h in 2:n)
      {
        fspec <- fspec + gamma[h]*2*cos((h-1)*freqs)
      }
    }
    mid <- TT - floor(TT/2)
    if(omega == 0)
    {
      low <- mid + 1
      hi <- mid + m
      spec.true <- sum(c(rev(gamma),gamma[-1]))
    }
    if(omega == pi)
    {
      hi <- TT
      low <- TT - m + 1
      gamma.alt <- gamma*(-1)^seq(0,n-1)
      spec.true <- sum(c(rev(gamma.alt),gamma.alt[-1]))
    }
    freqs.sub <- freqs[low:hi]
    fspec.sub <- fspec[low:hi]
    f0 <- sum(fspec.sub^2)/m
    f2 <- sum((omega-freqs.sub)^2*fspec.sub^2)/m
    f4 <- sum((omega-freqs.sub)^4*fspec.sub^2)/m
    g0 <- sum(fspec.sub)/m
    g2 <- sum((omega-freqs.sub)^2*fspec.sub)/m
    c2 <- sum((omega-freqs.sub)^2)/m
    c4 <- sum((omega-freqs.sub)^4)/m
    q <- c4 - c2^2
    avar <- f0*c4^2 - 2*c4*c2*f2 + c2^2*f4
    avar <- avar/(m*q^2)
    abias <- c4*g0 - c2*g2 
    abias <- abias/q - spec.true
    if(m == 1) 
    {
      avar <- f0 
      abias <- g0 - spec.true
    }
    val <- sqrt(avar + abias^2)
    return(val)
  }
  
  # case omega = 0
  fit <- optimize(mse,lower=0,upper=.5,gamma=gamma,omega=0,TT=TT)
  delta.0 <- fit$minimum
  delta.0 <- min(delta.0,.5)

  # case omega = pi
  fit <- optimize(mse,lower=0,upper=.5,gamma=gamma,omega=pi,TT=TT)
  delta.pi <- fit$minimum
  delta.pi <- min(delta.pi,.5)
  
  delta <- c(delta.0,delta.pi)
  return(delta)
}
