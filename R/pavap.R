pavap <- 
function(y, wmat, ngroups)
{
	n <- length(y)/ngroups
	sol <- c()
	w <- w2 <- list()
	for(i in 1:ngroups)
  	{
		w[[i]] <- diag(wmat)[1:n + (i - 1)*n]
		w2[[i]] <- array(colSums(wmat[1:n + (i - 1)*n, ]), c(n, ngroups))[, -i, drop = FALSE]
		sol <- c(sol, pava(y[1:n + (i - 1)*n], w[[i]]))
	}  
	dif <- (y - sol) %*% wmat %*% (y - sol)
	ii <- 0
	repeat
	{
		i <- ii%%ngroups + 1
		ii <- ii + 1
		ww <- rowSums((y - sol)[-(1:n + (i - 1)*n)] * w2[[i]])
		sol[1:n + (i - 1)*n] <- pava(y[1:n + (i - 1)*n] + ww/w[[i]], w[[i]])
		if(i == ngroups)
		{
			dif2 <- (y - sol) %*% wmat %*% (y - sol)
			if(dif - dif2 < .Machine$double.eps)
				break
			dif <- dif2
		}
	}
	return(list(x = sol))
}
