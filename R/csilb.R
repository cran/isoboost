csilb <-
function(xlearn, ...) UseMethod("csilb")


csilb.formula <-
function(formula, data, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("csilb")
    output$call <- call
    if(missing(data))
		stop("There is no xlearn set.\n\n")
    if(!is(class(data),"data.frame"))
		stop("data parameter is not a data.frame.\n\n")
    dataset <- model.frame(formula = formula, data = data)
    ylearn <- dataset[, 1]
    xlearn <- dataset[, -1, drop = FALSE]
    output <- csilb.default(xlearn, ylearn, ...)
    output$call <- call
    output
}


csilb.default <- 
function(xlearn, ylearn, xtest = xlearn, mfinal = 100, monotone_constraints = rep(0, dim(xlearn)[2]), prior = NULL, ...)
{
	output <- list()
    output$call <- match.call()
    if(missing(xlearn))
        stop("xlearn parameter is missing")
    if(missing(ylearn))
        stop("ylearn vector is missing.\n\n")
    check <- checks(xlearn, ylearn, xtest, mfinal, monotone_constraints, prior)
    if(is.null(check))
        return(output)
    xlearn <- check$xlearn
    ylearn <- check$ylearn
    xtest <- check$xtest
    mfinal <- check$mfinal
    ordering <- check$monotone_constraints
    trainset <- check$trainset
    prior <- check$prior
    rm(check)

	K <- length(levels(as.factor(ylearn)))
	cntrl <- rpart.control(maxdepth = 1, minsplit = 0, minbucket = 1, maxcompete = 0, cp = -1, xval = 0,
				maxsurrogate = 0, usesurrogate = 0)    
	orig_y <- cbind(1:length(ylearn), ylearn) 
	xlearn <- as.matrix(xlearn)
	rownames(xlearn) <- NULL
	xlearnn <- xlearn
	xlearnn[, (ordering == -1)] <- -xlearnn[, (ordering == -1)]
	N <- dim(xlearn)[1]
	d <- dim(xlearn)[2]
	plearn <- array(1/K, c(N, K))
	cumplearn <- array(0, c(N, K + 1))
	M <- dim(xtest)[1]
	H <- array(0, c(K - 1, K - 1))
	alpha <- -log((2:K - 1)/(K - 2:K + 1))
	score <- rep(0, K - 1)
	flearn <- Flearn <- rep(0, N)
	ff <- array(0, c(N, d))
	ftest <- Ftest <- rep(0, M)
	for(k in 1:K)
		cumplearn[, k] <- rowSums(plearn[, k:K, drop = FALSE])
	for(m in 1:mfinal)
	{
		cc <- cumplearn*(1 - cumplearn)
		w <- pmax((cc[, -(K + 1)] + cc[, -1])[orig_y] , 2*.Machine$double.eps)
		z <- (1 - cumplearn[, -(K + 1)] - cumplearn[, -1])[orig_y]/w
		z <- pmax(pmin(z, 4), -4)
		for(j in 1:d)
		{
			if(ordering[j] == 0)
			{
				xx <- xlearn[, j]
				fit <- rpart(z ~ xx, weights = w, control = cntrl, method = "anova")
				ff[, j] <- predict(fit)
			}
			else
				ff[, j] <- gpava(xlearnn[, j], z, weights = w)$x
		}
		min <- which.min(colSums(w*(ff - z)^2))
		flearn <- ff[, min]  
		if(ordering[min] == 0)
		{
			xx <- xlearn[, min]
			fit <- rpart(z ~ xx, weights = w, control = cntrl, method = "anova")
			xx <- xtest[, min]
			ftest <- predict(fit, newdata = data.frame(xx))
		}
		else if(ordering[min] > 0)
			ftest <- estimate(xtest[, min], sort(xlearn[, min]), sort(flearn))
		else
			ftest <- estimate(xtest[, min], sort(xlearn[, min]), rev(sort(flearn)))
		Ftest <- Ftest + ftest
		Flearn <- Flearn + flearn
		if(K > 2)
		{
			f <- function(x)
			{
				fF <- t(t(array(rep(Flearn, K - 1), c(N, K - 1))) + x)
				fpcum <- cbind(1, 1/(1 + exp(-fF)), 0)
				fplearn <- fpcum[, -(K + 1)] - fpcum[, -1]
				fplearn <- pmax(fplearn, exp(-745))
				res <- -sum(log(fplearn[orig_y]))
				res
			}
			alpha <- suppressWarnings(nlm(f, alpha)$estimate)
		}
		Falpha <- t(t(array(rep(Flearn, K - 1), c(N, K - 1))) + alpha)
		cumplearn <- cbind(1, 1/(1 + exp(-Falpha)), 0)
		plearn <- cumplearn[, -(K + 1)] - cumplearn[, -1]
		plearn <- pmax(plearn, exp(-745))
	}

	Ftest <- t(t(array(rep(Ftest, K - 1), c(M, K - 1))) + alpha)
	pcum <- cbind(1, 1/(1 + exp(-Ftest)), 0)
	ptest <- pcum[, -(K + 1)] - pcum[, -1]
	ptest <- t(t(ptest) * prior)
	ptest <- ptest/rowSums(ptest)
	apparent <- 100 * (1 - mean(apply(t(t(plearn) * prior), 1, which.max) == ylearn))
	output$trainset <- trainset
	output$prior <- prior
	output$apparent <- apparent
	output$mfinal <- m
	output$loglikelihood <- sum(log(plearn[orig_y]))
	output$posterior <- ptest
	output$class <- apply(ptest, 1, which.max)
	output
}
