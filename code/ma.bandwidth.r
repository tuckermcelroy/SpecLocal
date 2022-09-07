ma.bandwidth <- function(data,alpha,seasonal,demean=TRUE)
{
	n <- length(data)
	gamma.hat <- acf(data,lag=n-1,type="covariance",plot=FALSE,demean=demean)$acf[,,1]
	crit <- qnorm(1-alpha/2)
	K.n <- 1 + floor(3*sqrt(log(n,base=10)))
	K.n <- max(K.n,seasonal)
	rho.hat <- gamma.hat[-1]/gamma.hat[1]
	rho.test <- rho.hat/sqrt(log(n,base=10)/n)
	k <- 1
	while(k < (n-K.n))
	{
		if(max(abs(rho.test[k:(k+K.n-1)])) < crit) 
		{ 
			q.hat <- k-1
			k <- n-K.n
		} else { k <- k+1 }
	}
	return(q.hat)
}
