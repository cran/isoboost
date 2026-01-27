asilb <-
function(xlearn, ...) UseMethod("asilb")


asilb.formula <-
function(formula, data, ...)
{
    output <- list()
    call <- match.call()
    call[[1L]] <- as.name("asilb")
    output$call <- call
    if(missing(data))
        stop("There is no xlearn set.\n\n")
    if(is.null(data))
        stop("There is no xlearn set.\n\n")
    if(!is(class(data),"data.frame"))
       stop("data parameter is not a data.frame.\n\n")
    dataset <- model.frame(formula = formula, data = data)
    ylearn <- dataset[, 1]
    xlearn <- dataset[, -1, drop = FALSE]
    output <- asilb.default(xlearn, ylearn, ...)
    output$call <- call
    output
}


asilb.default <- 
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
	Y <- 1*outer(ylearn, 2:K, ">=") 
	colnames(Y) <- 2:K
	xlearn <- as.matrix(xlearn)
	rownames(xlearn) <- NULL
	xlearnn <- xlearn
	xlearnn[, (ordering == -1)] <- -xlearnn[, (ordering == -1)]
	N <- dim(xlearn)[1]
	d <- dim(xlearn)[2]
	plearn <- array(1/K, c(N, K))
	flearn <- Flearn <- array(0, c(N, K))
	ff <- array(0, c(N, d))
	ftest <- Ftest <- array(0, c(dim(xtest)[1], K))
	for(m in 1:mfinal)
	{
		pp <- 0 
		for(k in K:2)
		{
			pp <- pp + plearn[, k]
			w <- pmax(pp*(1 - pp), 2*.Machine$double.eps)
			z <- Y[, k - 1]
			z[Y[, k - 1] == 1] <- 1/pp[Y[, k - 1] == 1]
			z[Y[, k - 1] == 0] <- -1/(1 - pp[Y[, k - 1] == 0])
			z <- pmin(pmax(z, -4), 4)      
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
			flearn[, k] <- ff[, min]
			if(ordering[min] == 0)
			{
				xx <- xlearn[, min]
				fit <- rpart(z ~ xx, weights = w, control = cntrl, method = "anova")
				xx <- xtest[, min]
				ftest[, k] <- predict(fit, newdata = data.frame(xx))
			}
			else
			{
				if(ordering[min] > 0)
					ftest[, k] <- estimate(xtest[, min], sort(xlearn[, min]), sort(flearn[, k]))
				else
					ftest[, k] <- estimate(xtest[, min], sort(xlearn[, min]), rev(sort(flearn[, k])))
			}
		}
		Ftest <- Ftest + ftest
		Flearn <- Flearn + flearn
		Flearnn <- Flearn
		for(k in 2:K)
			Flearnn[, k] <- Flearnn[, k - 1] + Flearnn[, k]
		Flearnn <- Flearnn - apply(Flearnn, 1, max)
		plearn <- exp(Flearnn)/rowSums(exp(Flearnn))
	}
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
	output$loglikelihood <- sum(log(plearn[cbind(1:N, ylearn)]))
	output$posterior <- ptest
	output$class <- apply(ptest, 1, which.max)
	output
}
