estimate <- 
function(xx, x, y)
{
	n <- length(x)
	nn <- length(xx)
 	nlower <- rowSums(t(array(x <= rep(xx, each = n), c(n, nn))))
 	ngreater <- rowSums(t(array(x >= rep(xx, each = n), c(n, nn))))
	equal <- nlower + ngreater - n
	y <- c(y[1], y, y[length(y)])
	x <- c(-1/.Machine$double.eps, x, 1/.Machine$double.eps)
	out <- (y[nlower + 1]*(x[nlower + 2] - xx) + y[nlower + 2]*(xx - x[nlower + 1]))/(x[nlower + 2] - x[nlower + 1])
	if(sum(equal) > 0)
		out[equal > 0] <- mapply(function(z1, z2) mean(y[z1 + seq_len(z2) + 1]), n - ngreater[equal > 0], equal[equal > 0])  
  	out
}
