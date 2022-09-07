parzen.taper <- function(x)
{
	x <- abs(x)
	val <- 0
	if(x < .5)
	{
		val <- 1 - 6*x^2 + 6*x^3
	} else
	{
		if(x < 1)
		{
			val <- 2*(1-x)^3
		}
	}
	return(val)
}
