amilb <-
function(xlearn, ...) UseMethod("amilb")


amilb.formula <-
function(formula, data, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("amilb")
    output$call <- call
    if(missing(data))
        stop("There is no xlearn set.\n\n")
    if(is.null(data))
        stop("There is no xlearn set.\n\n")
    if(class(data) != "data.frame")
       stop("data parameter is not a data.frame.\n\n")
    dataset <- model.frame(formula = formula, data = data)
    ylearn <- dataset[, 1]
    xlearn <- dataset[, -1, drop = FALSE]
    output <- amilb.default(xlearn, ylearn, ...)
    output$call <- call
    output
}


amilb.default <- 
function(xlearn, ylearn, xtest = xlearn, mfinal = 100, monotone_constraints = rep(0, dim(xlearn)[2]), prior = NULL, ...)
{
	output <- list()
    output$call <- match.call()
    if(missing(xlearn))
        stop("xlearn parameter is missing")
    if(missing(ylearn))
        stop("ylearn vector is missing.\n\n")
    check <- checks(xlearn, ylearn, xtest, mfinal, monotone_constraints, prior, ndigits = 5)
    if(is.null(check))
        return(output)
    xlearn <- check$xlearn
    ylearn <- check$ylearn
    xtest <- check$xtest
    mfinal <- check$mfinal
    ordering <- check$monotone_constraints
    ndigits <- check$ndigits
    trainset <- check$trainset
    prior <- check$prior
    rm(check)

	K <- length(levels(as.factor(ylearn)))
	Y <- 1*outer(ylearn, 2:K, ">=") 
	orig_y <- cbind(1:length(ylearn), ylearn)
	xlearn <- as.matrix(xlearn)
	xtest <- as.matrix(xtest)
	rownames(xlearn) <- rownames(xtest) <- NULL
	if(sum(ordering == 0) > 0)
	{
		xlearn <- cbind(xlearn, 1)
		xtest <- cbind(xtest, 1)
		ordering <- c(ordering, 0)
	}
	xlearnn <- xlearn
	xlearnn[, (ordering == -1)] <- -xlearnn[, (ordering == -1)]
	N <- dim(xlearnn)[1]
	d <- dim(xlearnn)[2]
	plearn <- array(1/K, c(N, K))
	regg <- reg <- c()
	numord <- sum(ordering != 0)
	numberord <- c()
	numnoord <- sum(ordering == 0)
	ord <- (1:d)[ordering != 0]
	noord <- (1:d)[ordering == 0]
	Xord <- Xnoord  <- c()
	for(k in 1:(K - 1))
	{
		Xord <- c(Xord, ord + d*(k - 1))
		Xnoord <- c(Xnoord, noord + d*(k - 1))
	}
	Xordering <- rep(ordering, K - 1)
	YY <- c(Y)
	XX <- array(0, c((K - 1)*N, (K - 1)*d))
	for(j in 1:(K - 1))
		XX[N*(j - 1) + 1:N, d*(j - 1) + 1:d] <- xlearnn
	if(numord > 0 && K > 2)
	{
		ranges <- c()
		for(j in 1:numord)
			ranges <- c(ranges, 0.5 + diff(range(xlearnn[, ord[j]])))
		xlearnord <- xlearnn[, ord, drop = FALSE]
		for(k in 3:K)
			xlearnord <- rbind(xlearnord, xlearnn[, ord, drop = FALSE] + rep(ranges, each = N)*(k - 2))
	}  
	M <- dim(xtest)[1]
	Xtest <- array(0, c((K - 1)*M, (K - 1)*d))
	for(j in 1:(K - 1))
		Xtest[(M*(j - 1) + 1):(M*j), (d*(j - 1) + 1):(d*j)] <- xtest
	W <- array(0, c(N*(K - 1), N*(K - 1)))
	dfl <- numord + 1*(numnoord > 0)
	flearn <- array(0, c(N*(K - 1), dfl))
	Flearn <- array(0, c(N, K))
	FF <- c(Flearn[, -1])
	log_lik <- sum(log(plearn[orig_y]))
	error <- 0
	for(m in 1:mfinal)
	{
		Plearn <- plearn
		for(k in 1:(K - 1))
			Plearn[, k] <- rowSums(plearn[, k:K])
		PP <- c(Plearn[ , -1])
		for(i in 1:(K - 1))
			for(j in 1:i)
			{
				if(j == i)
					diag(W[N*(j - 1) + 1:N, N*(j - 1) + 1:N]) <- Plearn[, i + 1]*(1 - Plearn[, i + 1])
				else  
					diag(W[N*(i - 1) + 1:N, N*(j - 1) + 1:N]) <- diag(W[N*(j - 1) + 1:N, N*(i - 1) + 1:N]) <- Plearn[, i + 1]*(1 - Plearn[, j + 1])
			}
		C <- tryCatch(chol(W), error = function(cond) return(NA))
		if(is.na(C[1]))
		{
			error <- 1
			break
		}
		direction <- tryCatch(solve(W, YY - PP), error = function(cond) return(NA))
		if(is.na(direction[1]))
		{
			error <- 2
			break
		}
		CXX <- C %*% XX        
		flearnn <- flearn*0
		step <- 2
		repeat
		{    
			step <- step/2
			if(step < 1e-10)
			{
				error <- 3
				log_lik2 <- NULL
				break
			}
			z <- FF + step * direction
			dev <- sum(z^2)
			if(numord == 0)
			{
				regg <- lm(C%*%z ~ CXX - 1)
				Fsum <- XX %*% regg$coefficients
			}
			else
			{
				if(K == 2)
				{
					#Backfitting:
					for(iii in 1:500)
					{
						for(j in 1:numord)
						{
							zz <- z - rowSums(flearnn[, -j, drop = FALSE])
							flearnn[, j] <- gpava(xlearnn[, ord[j]], zz, weights = diag(W))$x
						}
						if(numnoord > 0)
						{
							zz <- z - rowSums(flearnn[, -dfl, drop = FALSE])
							regg <- lm(C%*%zz ~ CXX[, Xnoord] - 1)
							flearnn[, dfl] <- c(XX[, Xnoord] %*% regg$coefficients)
						}           
						dev <- c( t(z - rowSums(flearnn)) %*% W %*% (z - rowSums(flearnn)), dev)[1:2]
						if(dev[1] < .Machine$double.eps || dev[1]/dev[2] > 0.999 ||
							floor(dev[1] * 10^(ndigits - ceiling(abs(pmax(log10(dev[1]), -1))))) == 
							floor(dev[2] * 10^(ndigits - ceiling(abs(pmax(log10(dev[2]), -1))))))
							break
					}
				}
				else
				{
					#Backfitting:
					for(iii in 1:500)
					{
						for(j in 1:numord)
						{
							zz <- z - rowSums(flearnn[, -j, drop = FALSE])
							ordsort <- order(xlearnord[, j], zz)
							flearnn[ordsort, j] <- pavap(zz[ordsort], W[ordsort, ordsort], K - 1)$x
						}
						if(numnoord > 0)
						{
							zz <- z - rowSums(flearnn[, -dfl, drop = FALSE])
							regg <- lm(C%*%zz ~ CXX[, Xnoord] - 1)
							flearnn[, dfl] <- c(XX[, Xnoord] %*% regg$coefficients)
						}
						dev <- c( t(z - rowSums(flearnn)) %*% W %*% (z - rowSums(flearnn)), dev)[1:2]
						if(dev[1] < .Machine$double.eps || dev[1]/dev[2] > 0.999 || 
							floor(dev[1] * 10^(ndigits - ceiling(abs(pmax(log10(dev[1]), -1))))) == 
							floor(dev[2] * 10^(ndigits - ceiling(abs(pmax(log10(dev[2]), -1))))))
							break
					}
				}
				Fsum <- rowSums(flearnn)
			}
			Flearn <- cbind(0, array(F, c(N, K - 1)))
			for(k in 2:K)
				Flearn[, k] <- Flearn[, k - 1] + Flearn[, k]
			Flearn <- Flearn - apply(Flearn, 1, max)    
			plearnn <- exp(Flearn)/rowSums(exp(Flearn))
			log_lik2 <- sum(log(plearnn[orig_y]))
			if(log_lik2 >= log_lik)
			{
				flearn <- flearnn
				plearn <- plearnn
				reg <- regg
				FF <- Fsum
				break
			}  
		}
		if(is.null(log_lik2) || floor(-log_lik * 10^(8 - ceiling(abs(pmax(log10(-log_lik), -1))))) == 
 			floor(-log_lik2 * 10^(8 - ceiling(abs(pmax(log10(-log_lik2), -1))))))
			break
		log_lik <- log_lik2
	}    
	Ftest <- rep(0, M*(K - 1))
	if(numnoord > 0)
		Ftest <- Xtest[, Xnoord] %*% reg$coefficients
	if(numord > 0)
		for(k in 1:(K - 1))
			for(j in 1:numord)
			{
				if(ordering[ord[j]] > 0)
					Ftest[(k - 1)*M + 1:M] <- Ftest[(k - 1)*M + 1:M] + 
					estimate(xtest[, ord[j]], sort(xlearn[, ord[j]]), sort(flearn[(k - 1)*N + 1:N, j]))
				else
					Ftest[(k - 1)*M + 1:M] <- Ftest[(k - 1)*M + 1:M] + 
					estimate(xtest[, ord[j]], sort(xlearn[, ord[j]]), rev(sort(flearn[(k - 1)*N + 1:N, j])))
			}
	Ftest <- cbind(0, array(Ftest, c(M, K - 1)))
	Fout <- Ftest    
	for(k in 2:K)
		Ftest[, k] <- Ftest[, k - 1] + Ftest[, k]
	# No inf/inf:
	Ftest <- Ftest - apply(Ftest, 1, max)
	ptest <- exp(Ftest)/rowSums(exp(Ftest))
	ptest <- t(t(ptest) * prior)
	ptest <- ptest/rowSums(ptest)
	apparent <- 100 * (1 - mean(apply(t(t(plearn) * prior), 1, which.max) == ylearn))
	output$trainset <- trainset
	output$prior <- prior
	output$apparent <- apparent
	output$mfinal <- m
	output$loglikelihood <- log_lik
	output$posterior <- ptest
	output$class <- apply(ptest, 1, which.max)
	output
}
