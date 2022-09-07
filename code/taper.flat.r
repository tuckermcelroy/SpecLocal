taper.flat <- function(x,c)
{
	if (x <= c) { taperVal <- 1 }
	else { taperVal <- (1-x)/(1-c) }
	taperVal <- max(taperVal,0)
	return(taperVal)		
}
