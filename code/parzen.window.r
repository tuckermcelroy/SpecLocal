parzen.window <- function(z)
{
	z <- abs(z)
	val <- (3/(8*pi))*(4*sin(z/4)/z)^4
	return(val)
}

