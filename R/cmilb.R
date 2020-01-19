cmilb <-
function(xlearn, ...) UseMethod("cmilb")


cmilb.formula <-
function(formula, data, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("cmilb")
    output$call <- call
    if(missing(data))
        stop("There is no xlearn set.\n\n")
    if(is.null(data))
        stop("There is no xlearn set.\n\n")
    dataset <- model.frame(formula = formula, data = data)
    ylearn <- dataset[, 1]
    xlearn <- dataset[, -1, drop = FALSE]
    output <- cmilb.default(xlearn, ylearn, ...)
    output$call <- call
    output
}


cmilb.default <- 
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
	cumplearn <- array(0, c(N, K + 1))
	regg <- reg <- c()
	numord <- sum(ordering != 0)
	numberord <- c()
	numnoord <- sum(ordering == 0)
	ord <- (1:d)[ordering != 0]
	noord <- (1:d)[ordering == 0]
	M <- dim(xtest)[1]
	alpha0 <- alpha <- -log((2:K - 1)/(K - 2:K + 1))
	H <- array(0, c(K - 1, K - 1))  
	score <- rep(0, K - 1)
	FF <- Fsum <- rep(0, N)
	dfl <- numord + 1*(numnoord > 0)
	flearn <- array(0, c(N, dfl))
	log_lik <- sum(log(plearn[orig_y]))
	error <- 0
	for(m in 1:mfinal)
	{
		for(k in 1:K)
			cumplearn[, k] <- rowSums(plearn[, k:K, drop = FALSE])
		cc <- cumplearn*(1 - cumplearn)
		w <- (cc[, -(K + 1)] + cc[, -1])[orig_y]
		if(sum(w == 0) > 0)
			break      
		flearnn <- flearn*0
		Direction <- (1 - cumplearn[, -(K + 1)] - cumplearn[, -1])[orig_y]/w
		Step <- 2
		repeat
		{
			Step <- Step/2
			if(Step < 1e-10)
			{
				error <- 1
				log_lik2 <- NULL
				break
			}    
			z <- FF + Step*Direction
			dev <- sum(z^2)
			if(numord == 0)
			{
				regg <- lm(z ~ xlearnn - 1)
				Fsum <- xlearnn %*% regg$coefficients
			}
			else
			{
				for(iii in 1:500)
				{
					for(j in 1:numord)
					{
						zz <- z - rowSums(flearnn[, -j, drop = FALSE])
						flearnn[, j] <- gpava(xlearnn[, ord[j]], zz, weights = w)$x
					}
					if(numnoord > 0)
					{
						zz <- z - rowSums(flearnn[, -dfl, drop = FALSE])
						regg <- lm(zz ~ xlearnn[, noord] - 1)
						flearnn[, dfl] <- c(xlearnn[, noord] %*% regg$coefficients)
					}       
					dev <- c(sum(w*((z - rowSums(flearnn))^2)), dev)[1:2]

					if(dev[1] < .Machine$double.eps || dev[1]/dev[2] > 0.999 ||
						floor(dev[1] * 10^(ndigits - ceiling(abs(pmax(log10(dev[1]), -1))))) == 
						floor(dev[2] * 10^(ndigits - ceiling(abs(pmax(log10(dev[2]), -1))))))
						break
				}
				Fsum <- rowSums(flearnn)
			}
			if(K > 2)
			{
				f <- function(x)
				{
					fF <- t(t(array(rep(F, K - 1), c(N, K - 1))) + x)
					fpcum <- cbind(1, 1/(1 + exp(-fF)), 0)
					fplearn <- fpcum[, -(K + 1)] - fpcum[, -1]
					fplearn <- pmax(fplearn, exp(-745))
					res <- -sum(log(fplearn[orig_y]))
					res
				}
				alphaa <- suppressWarnings(nlm(f, alpha)$estimate)
				Falpha <- t(t(array(rep(F, K - 1), c(N, K - 1))) + alphaa)
				pcum <- cbind(1, 1/(1 + exp(-Falpha)), 0)
				plearnn <- pcum[, -(K + 1)] - pcum[, -1]
				plearnn <- pmax(plearnn, exp(-745))
				log_lik2 <- sum(log(plearnn[orig_y]))
				if(log_lik2 >= log_lik)
				{
					alpha <- alphaa
					FF <- Fsum
					flearn <- flearnn
					plearn <- plearnn
					reg <- regg
					break
				}   
			}
			else
			{
				Falpha <- array(rep(F, K - 1), c(N, K - 1))
				pcum <- cbind(1, 1/(1 + exp(-Falpha)), 0)
				plearnn <- pcum[, -(K + 1)] - pcum[, -1]
				log_lik2 <- sum(log(plearnn[orig_y]))      
				if(log_lik2 >= log_lik)
				{
					FF <- Fsum
					flearn <- flearnn
					plearn <- plearnn
					reg <- regg
					break
				}
			}      
		}    
		if(is.null(log_lik2) || floor(-log_lik * 10^(8 - ceiling(abs(pmax(log10(-log_lik), -1))))) == 
			floor(-log_lik2 * 10^(8 - ceiling(abs(pmax(log10(-log_lik2), -1))))))
			break
		log_lik <- log_lik2
	}    
	Ftest <- rep(0, M)
	if(numnoord > 0)
		Ftest <- xtest[, noord] %*% reg$coefficients
	if(numord > 0)
		for(j in 1:numord)
		{
			if(ordering[ord[j]] > 0)
				Ftest <- Ftest + estimate(xtest[, ord[j]], sort(xlearn[, ord[j]]), sort(flearn[, j]))
			else
				Ftest <- Ftest + estimate(xtest[, ord[j]], sort(xlearn[, ord[j]]), rev(sort(flearn[, j])))
		}
	Ftest <- t(t(array(rep(Ftest, K - 1), c(M, K - 1))) + alpha)
	Fout <- cbind(0, Ftest)
	pcum <- cbind(1, 1/(1 + exp(-Ftest)), 0)
	ptest <- pcum[, -(K + 1)] - pcum[, -1]
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