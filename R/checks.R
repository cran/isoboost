checks <-
function(xlearn, ylearn, xtest, mfinal, monotone_constraints, prior, ndigits = 5)
{   
    if(is.null(xlearn))
        stop("xlearn is NULL.\n\n")
    if(sum(!sapply(as.list(xlearn), typeof) %in% 
       c("integer", "double", "complex", "logical")) > 0)
        stop("Invalid xlearn set.\n\n")
    if(is.null(xtest))
        stop("xtest is NULL.\n\n")
    if(sum(!sapply(as.list(xtest), typeof) %in% 
       c("integer", "double", "complex", "logical")) > 0)
        stop("Invalid xtest set.\n\n")
    if(is.null(ylearn))
        stop("ylearn is NULL.\n\n")
    if(class(ylearn) == "factor")
    {
        if(sum(suppressWarnings(is.na(as.numeric(levels(ylearn))))) > 0)
            stop("If ylearn is a factor, its levels must be numbers.\n\n")
        ylearn <- as.numeric(levels(ylearn)[ylearn])
    }
    ylearn <- c(ylearn)
    if(length(ylearn) != dim(xlearn)[1])
        stop("xlearn and ylearn have different number of observations.\n\n")
    dimension <- dim(xlearn)[2]
    if(is.null(colnames(xlearn)))
    {
	    if(dim(xtest)[2] != dimension)
            stop("Missing variables in xtest set.\n\n")
    }
    if(!is.null(colnames(xlearn)) & !is.null(colnames(xtest)))
    {
	    if(sum(colnames(xlearn) %in% colnames(xtest)) < dimension)
            stop("Missing variables in xtest set.\n\n")
        else
	        xtest <- xtest[, colnames(xlearn)]
    }
    trainset <- cbind(as.matrix(xlearn), ylearn)
    if(sum(is.na(trainset)) > 0)
    {
        trainset <- trainset[complete.cases(trainset), , drop = FALSE]
        stop("Missing values in the training or grouping set have been deleted.\n\n")
    }
    ylearn <- trainset[, dim(trainset)[2]]
    xlearn <- trainset[, -dim(trainset)[2], drop = FALSE]
    numgroups <- length(levels(as.factor(ylearn)))
    if(sum(abs(as.numeric(levels(as.factor(ylearn))) - {1:numgroups})) > 0)
        stop("Classes must be identified by natural numbers varying from 1 to the number of classes.\n\n")
    if(numgroups == 1)
        stop("xlearn set has only one class.\n\n")
    if(is.null(monotone_constraints))
        stop("monotone_constraints vector must be specified.\n\n")
    if(length(monotone_constraints) != dimension)
        stop("monotone_c vector must have the same length as the number of explanatory variables in xlearn.\n\n")
    if(sum(!(monotone_constraints %in% c(-1, 0, 1))) > 0)
        stop("Invalid monotone_c vector.\n\n")  
    if(is.null(ndigits) | floor(ndigits) != ndigits | length(ndigits) > 1 | ndigits < 1)
        stop("ndigits must be a whole number.")
    if(!is.null(prior))
    {
        if(length(prior) != numgroups)
            stop("Wrong number of classes in a priori probabilities.\n\n")
        if(sum(prior > 1 | prior < 0) > 0 || abs(sum(prior) - 1) > 1e-12)
            stop("prior values must be in the interval [0, 1] and sum 1.\n\n")
    }
    if(is.null(prior))
        prior <- rep(1/numgroups, numgroups)
    colnames(xlearn) <- NULL
    colnames(xtest) <- NULL
    output <- list()
    output$xlearn <- xlearn
    output$ylearn <- ylearn
    output$trainset <- trainset
    output$xtest <- xtest
    output$mfinal <- mfinal
    output$monotone_constraints <- monotone_constraints
    output$ndigits <- ndigits
    output$prior <- prior
    output
}
