#'
#' function to split a vector in pieces, needed for e.g. cross validation
#' 
#' @param x the vector to split
#' @param n the number of splits
#' 
SplitVec <- function(x,n) split(x, factor(sort(rank(x)%%n)))


#'
#' Helper function for parallel processing to compute prediction error measure for different cuts of the variable
#'
#' @param c cutpoint for variable
#' @param response the response, see \code{link{OptimalCutPoint}
#' @param var the variable, see \code{link{OptimalCutPoint}
#' @param treat treatment, see \code{link{OptimalCutPoint}
#' @param covars covariates, see \code{link{OptimalCutPoint} 
#' @param data see \code{link{OptimalCutPoint}
#' @param K see \code{link{OptimalCutPoint}  
#' @param repeats see \code{link{OptimalCutPoint}
#' @param seed see \code{link{OptimalCutPoint}
#'
CutPointCV <- function(c, response, var, treat, covars, data, K, repeats, seed) {

    ## create dichotomized variable
    data.c <- data
    data.c$varc <- factor(ifelse(data[[var]] < c, "low", "high"), levels=c("low", "high"))

    ## set seed to garantee same splits for every cut
    set.seed(seed)

    ## matrix for prediction
    predictions <- matrix(NA, nrow=nrow(data), ncol=repeats)

    ## over all repeats
    for (r in 1:repeats) {

        ## get splits
        splits <- SplitVec(sample(1:nrow(data)), K)

        ## CV
        for (k in 1:K) {

            ## get train data and test data
            train.data <- splits[[k]]
            test.data <- setdiff(1:nrow(data.c), train.data)

            ##OptimalCutPoint get formula
            ## 1) Logistic Regression
            ## Note: use * to account for the interaction between the variable and the treatment
            if (!is.Surv(response)) {
                formula <- paste(response, paste("varc", treat, sep="*"), sep="~") 
            } else {
                formula <- paste("response~varc", treat, sep="*")
            }

            ## add covars
            if (!is.null(covars)) {
                formula <- paste(formula, paste(covars, collapse="+"), sep="+")
            } 
            formula <- as.formula(formula)

            ## fit model and predict
            if (!is.Surv(response)) {
                ## logistic regression 
                if (length(unique(na.omit(data.c$varc[train.data]))) > 1) {
                    lrm <- try(glm(formula, data=data.c, subset=train.data, family=binomial()))
                    if (!inherits(lrm, "try-error")) {
                        predictions[test.data, r] <- predict(lrm, newdata=data.c[test.data,], type="response")
                    }
                }
            } else {
                cox <- coxph(formula, data=data.c, subset=train.data)
            }

        }
    }

    if (!is.Surv(response)) {
        ## return median of AUCs: we have so many AUCs like we have repeats
        return(median(ComputeAUC(predictions, labels=data[[response]])))
    }

}


#'
#' Computes odds ratio and p-value of interaction term of variable and treatment
#'
#' Helper function for parallel processing
#'
#' @param c cutpoint for variable
#' @param response the response, see \code{link{OptimalCutPoint}
#' @param var the variable, see \code{link{OptimalCutPoint}
#' @param treat treatment, see \code{link{OptimalCutPoint}
#' @param covars covariates, see \code{link{OptimalCutPoint} 
#' @param data see \code{link{OptimalCutPoint}
#'
CutPointTests <- function(c, response, var, treat, covars, data) {

    ## create dichotomized variable
    data.c <- data
    data.c$varc <- factor(ifelse(data[[var]] < c, "low", "high"), levels=c("low", "high"))

    ##OptimalCutPoint get formula
    ## 1) Logistic Regression
    ## Note: use * to account for the interaction between the variable and the treatment
    if (!is.Surv(response)) {
        formula <- paste(response, paste("varc", treat, sep="*"), sep="~") 
    } else {
        formula <- paste("response~varc", treat, sep="*")
    }

    ## add covars
    if (!is.null(covars)) {
        formula <- paste(formula, paste(covars, collapse="+"), sep="+")
    } 
    formula <- as.formula(formula)

    ## fit model and predict
    if (!is.Surv(response)) {
        ## logistic regression 
        model <- try(glm(formula, data=data.c, family=binomial()))
    } else {
        model <- try(coxph(formula, data=data.c))
    }

    if (!inherits(model, "try-error")) {
        ## Extract OR and p-value from model for interaction effect
        # get summary
        ttt <- summary(model)$coefficients
        # get exponents of odds-ratios
        ttt[,1] <- exp(ttt[,1])
        # get interaction effect and return odds-ratio and p-value
        m <- intersect(grep(treat, rownames(ttt)), grep("varc", rownames(ttt)))[1] 
        return(ttt[m,c(1,4)])
    } else {
        ## if the model could not be fitted return two NAs
        return(rep(NA, 2))
    }

}


#'
#' Cut point search via cross-validation
#' 
#' Searches an optimal cut-point for a continuous variable 'var' in interaction with treatmen 'treat' for a binary response. The cut point is optimized via cross-validated AUC
#'
#' @param response the name of the binary response variable
#' @param var the name of continuous variable
#' @param treat the name of the treatment variable
#' @param covars optional vector with covariate names
#' @param data the data frame with the data
#' @param num.cuts the number of cuts to be tested. If the number of cuts is larger than the number of unique values...
#' @param K the number of folds of the cross-validation
#' @param repeats the nu
#' @param cores
#' @param seed the seed to be used for cross-validation
#'
#' @return matrix with cut points, median AUC and AUCs for every repeat
#'
#' @import parallel
#' @export
#'
OptimalCutPoint <- function(response, var, treat, covars=NULL, data, num.cuts=10000, K=10, repeats=5, cores=1, seed=123) {

    ## libraries
    require(parallel)

    ## get sorted uniqe var
    num.uni <-  sort(unique(na.omit(data[[var]])))

    ## get cuts
    ## if number of cuts > length(var) take every var as one cut except the first one
    if (num.cuts >= length(num.uni)) {
        cuts <- num.uni[-1]
    } else { 
        cuts <- quantile(num.uni, probs=seq(0, 1, length.out=(num.cuts+2)))
        cuts <- cuts[-c(1, length(cuts))]
    }   

    ## main loop parallel excuted over all posisble cuts
    cuts.errors <- mclapply(cuts, 
                            function(xx) {
                                
                                ## get test results and prediction error
                                coefs <- CutPointTests(c=xx, response=response, var=var, treat=treat, covars=covars, data=data)
                                error <- CutPointCV(c=xx, response=response, var=var, treat=treat, covars=covars, data=data, K=K, repeats=repeats, seed=seed)

                                return(c(coefs, error))

                            },mc.cores=cores)

    cuts.errors <- do.call("rbind", cuts.errors)

    if(!is.Surv(response)) {
        colnames(cuts.errors) <- c("OR", "P", "AUC")
        rownames(cuts.errors) <- NULL
        cuts.errors <- cbind(cuts, cuts.errors)
        return(cuts.errors)
    }

}

#'
#' Helper function to compute the AUCs from a prediction matrix
#'
#' @param predictions matrix with predictions, rows are the individuals, cols the repeats
#' @param labels true class labels
#' @param correct
#'
#' @return vector with AUCs, length equal to the number of cols
#'
ComputeAUC <- function(predictions, labels, correct=T) {

    ## libraries
    require(caTools)

    ## get AUCs
    aucs <- try(colAUC(X=predictions, y=labels))
    if (!inherits(aucs, "try-error")) {
        ## correct AUCs < 0.5
        if (correct) {
            aucs <- pmax(aucs, 1-aucs)
        }
    } else {
        aucs <- rep(NA, ncol(predictions))
    }

    return(aucs)
}


#'
#' Helper function to get a tail-orinted (upwards: starting from few to all patients) subpopulation pattern
#'
#' @param x vector with covariate
#' @param g desired number of groups
#'
#' @return list with k elements, each holding the indices of the elements from x belonging to that group
#'
#' @examples
#' set.seed(123)
#' x <- runif(200)
#' STEPPTailUp(x, 10)
#' 
STEPPTailUp <- function(x,g) {

    ## needs Hmisc
    require(Hmisc)

    ## get intervals
    cuts <- cut2(x, g=g)
    i <- lapply(levels(cuts), function(xx) which(cuts==xx))
    names(i) <- levels(cuts)

    ## get sucessive union of intervals as final grouping
    groups <- lapply(seq_along(i), function(xx) unique(unlist(i[1:xx])))

    ## construct names of groups: range if x in this groups
    int <- sapply(groups, 
                  function(xx) {
                      paste("<=", signif(max(x[xx], na.rm=T) ,3), sep="")
                  })
    ## set label of last window that contains all samples to 'All'
    int[length(int)] <- "All"
    names(groups) <- int 

    # return groups
    return(groups)
}

#'
#' Helper function to get a tail-orinted (downwards: starting from all to few patients) subpopulation pattern
#'
#' @param x vector with covariate
#' @param g desired number of groups
#'
#' @return list with k elements, each holding the indices of the elements from x belonging to that group
#'
#' @examples
#' set.seed(123)
#' x <- runif(200)
#' STEPPTailDown(x, 10)
#' 
STEPPTailDown <- function(x,g) {

    ## needs Hmisc
    require(Hmisc)

    ## get intervals
    cuts <- cut2(x, g=g)
    i <- lapply(levels(cuts), function(xx) which(cuts==xx))
    names(i) <- levels(cuts)

    ## get sucessive union of intervals as final grouping
    groups <- rev(lapply(seq_along(i), function(xx) unique(unlist(rev(i)[1:xx]))))

    ## construct names of groups: range if x in this groups
    int <- sapply(groups, 
                  function(xx) {
                      paste(">=", signif(min(x[xx], na.rm=T), 3), sep="")
                  })
    ## set label of first window that contains all samples to 'All'
    int[1] <- "All"
    names(groups) <- int 

    # return groups
    return(groups)
}

#'
#' Helper function to get a tail-oriented (both parts) subpopulation pattern
#' 
#' @param x vector with covariate
#' @param k number of groups. Should be odd (2*g+1), if even (k+1) is used
#'
STEPPTail <- function(x,k) {

    ## check if k is an odd number
    ## if even take k+1
    if (bitwAnd(k, 1)==0) {
        warning("STEPPTail: number of groups k is even, taking k+1 instead!")
        k <- k + 1
    }

    ## get first tail with rising values of x
    groups.up <- STEPPTailUp(x, as.integer((k+1)/2))
    groups.down <- STEPPTailDown(x, as.integer((k+1)/2))[-1]

    ## return list
    return(list(groups.up=groups.up, groups.down=groups.down))

}


#'
#' Helper function to get a sliding-window subpopulation pattern
#'
#' @param x vector with covariate
#' @param n2 size of each window (number of patients in one window)
#' @param n1 min difference between windows ('slide width')
#'
#' @return list with k elements, each holding the indices of the elements from x belonging to that group
#'
#' @examples
#' set.seed(123)
#' x <- runif(200)
#' STEPPSliding(x, n2=10, n1=5)
#'
#' @references Bonetti2000_StatInMed
#' 
STEPPSliding <- function(x,n2,n1) {

    ## needs Hmisc
    require(Hmisc)

    ## get x ordered and tabled
    x.tab <- table(sort(x))
    eta <- sort(unique(na.omit(x)))
    eta.max <- max(eta, na.rm=T)
    
    ## define n1 and n2
    n <- length(x)
    # new n1 is the number of overlaps between two consecutive windows
    n1 <- n2-n1
    if(n1 >= n | n2 >=n) {
        stop("n2 and (n1-n2) must be < number of samples")
    }

    ## get estimate of number of groups
    ## for allocation of vectors
    k <- ceiling(1+(1-n2/n)/(n2/n-n1/n))


    ## allocate list with groups and character vector holding the names of the groups
    groups <- vector("list", k)
    intervals <- rep("", k)

    ## initialize window borders w.min, w.max and counter b
    w.min <- min(x, na.rm=T)-1
    w.max <- w.min
    b <- 0
    while (w.max<eta.max) {

        ## increase counter b
        b <- b + 1

        ## get window maximum
        # get orderd unique x greater than window min
        eta.b <- eta[eta>w.min]
        # get cumulated numbers of items greater than window min
        sum.b <- cumsum(x.tab[eta>w.min])
        # define maximum as minimum of items over window size
        # if there is no such number take eta.max
        w.max <- min(eta.b[sum.b>=n2][1], eta.max, na.rm=T)

        ## get groups and definition of the interval
        groups[[b]] <- which(x > w.min & x <= w.max)
        intervals[b] <- paste("(", signif(w.min,3), ",", signif(w.max,3), "]", sep="")

        ## update window minimum for next step
        # get all etas smaller than window maximum
        eta.b <- eta[eta<=w.max]
        # get cumulated sum (reverse) of counts where eta is \leq window max
        sum.b <- rev(cumsum(rev(x.tab[eta <= w.max])))
        w.min <- min(eta.b[sum.b <= n1])
    }
    

    ## get groups and set names
    ## Note: group was defined as list with length k, might be different
    ##       take first b elements
    groups <- groups[1:b]
    names(groups) <- intervals[1:b]

    ## return groups
    return(groups)

}



#'
#' Coef plot from a STEPP analysis
#'
#' @param coefs matrix with coefficients and CI
#' @param line for 'smooth' a lowess curve is plotted through the coefficients, 'simple' connects all the points with a line, 'none' omits the line
#' @param hline intercept for help lines, if NULL no lines will be drawn. If missing the log parameter is used to draw a line either at 0 (log=TRUE) or 1 (log=FALSE)
#' @param xlab the x axis label
#' @param xlab the y axis label
#' @param log indicates wheter the coefficients are on a log scale.coef
#' @param log.scale indicates wheter the y axis should use a log scale. Will be ignored if the data are already on a log scale 
#'
#' @return ggplot object with the coefficient plot
#'
STEPPCoefPlot <- function(coefs, groups=NULL, line=c("smooth", "simple", "none"), hline, xlab="x", ylab="OR", main="", log, log.scale=T, ylim=NULL) {

    ## load ggplot2
    require(ggplot2)
    require(scales)

    ## match arguments
    line <- match.arg(line)

    ## check groups
    if (!is.null(groups) && nrow(coefs) != length(groups)) {
        warning("STEPPCoefPlot: length of group vector must match number of rows of coefficient matrix! Ignoring group vector!")
        groups <- NULL
    }

    ## create data frame for plotting
    data.plot <- as.data.frame(data.frame(coefs))
    if (!is.null(groups)) {
        data.plot$groups <- groups
    } else {
        data.plot$groups <- 1:nrow(data.plot)
    }

    ## main plot coefficients and confidende band
    pl <- ggplot(data.plot, aes(groups, betas)) + geom_point() + geom_ribbon(data=data.plot, aes(ymin=lower, ymax=upper), alpha=0.3)

    ## add line according to 'line' argument
    if (line=="smooth") {
        pl <- pl + geom_smooth(color="blue", se=F, method="loess")
    } else if (line=="simple") {
        pl <- pl + geom_line(color="blue")
    }

    ## add help lines if given
    if (missing(hline)) {
        hline <- if (log) 0 else 1 
    }
    if (!is.null(hline)) {
        pl <- pl + geom_hline(yintercept=hline, linetype=2, alpha=0.5)
    }

    ## x axis: add group labels and axis label
    if (is.null(groups)) {
        pl <- pl + scale_x_discrete(breaks=c(1:nrow(coefs)), labels=rownames(coefs)) + theme(axis.text.x=element_text(angle=90, vjust=0.5))
    }
    pl <- pl + xlab(xlab)

    ## y axis: log scale and axis label
    if (log.scale && !log) {
        pl <- pl + scale_y_continuous(trans=log2_trans(), breaks=trans_breaks("log2", function(xx) 2^xx))
        ylab <- paste(ylab, "(log2 scale)")
    } else if (log.scale && log) {
        warning("Coefficients in STEPP analyis already in log scale, ignoring 'log.scale' argument.")
    }
    pl <- pl + ylab(ylab)

    ## add ylim if set
    if (!is.null(ylim)) {
       pl <- pl + coord_cartesian(ylim=ylim) 
    }

    ## add title
    pl <- pl + ggtitle(main)

    ## adjust margins
    pl <- pl + theme(plot.margin = unit(c(1,1,1.5,0.5), "lines"))

    ##  return plot
    return(pl)

}

#'
#' Helper function to produce a plot that displayes the subpopulations (groups) and number of patients
#'
#' @param groups list with groups
#' @param intervals two-dimensional matrix with lower and upper bound of the groups. If none is given (default) the names of the group list are parsed
#' @param xlab the x axis label
#' @param ylab the y axis label
#' @param main the plot title. Default to NULL (no title set)
#'
#' @return ggplot object with the group plot
#' 
STEPPGroupPlot <- function(groups, intervals=NULL, xlab="Subpopulations [No]", ylab="y", main=NULL) {

    ## load ggplot
    require(ggplot2)

    ## create data frame for plotting
    nums <- sapply(groups, length)
    data.plot <- data.frame(groups=factor(paste("G", seq_along(groups), sep="")),
                            id=seq_along(groups))
    #data.plot$x1 <- data.plot$as.numerocgroups-0.1
    #data.plot$x2 <- data.plot$groups+0.1

    ## add intervals:
    ## 1) parse names of group list (default)
    ## 2) take intervals if given
    if (is.null(intervals)) {
        ttt <- gsub("[() \\]\\[]", "", names(groups), perl=T)
        ttt <- strsplit(ttt, ",")
        intervals <- matrix(as.numeric(unlist(ttt)), nrow=length(groups), ncol=2, byrow=T)
    }
    else {
        if (!all(dim(intervals)==c(length(groups),2))) {
            stop("intervals must have two columns and the same number of rows as groups!")
        }
    }
    colnames(intervals) <- c("lo", "up")
    data.plot <- cbind(data.plot, as.data.frame(intervals))

    ## main plot
    pl <- ggplot(data.plot, aes(groups)) + geom_rect(aes(x=groups, xmin=id-0.2, xmax=id+0.2, ymin=lo, ymax=up))
    
    ## x axis: add numbers in groups at x axis
    pl <- pl + scale_x_discrete(breaks=levels(data.plot$groups), labels=as.character(nums)) +  theme(axis.text.x=element_text(angle=90, vjust=0.5)) + xlab(xlab)

    ## y axis:
    pl <- pl + ylab(ylab)

    ## plot title
    if (!is.null(main)) {
        pl <- pl + ggtitle(main)
    }

    ## adjust margins
    pl <- pl + theme(plot.margin = unit(c(0,1,0.5,0.5),"lines"))

    ## return plot
    return(pl)
}



#'
#' Helper function to produce a boxplot of the continuous variable
#' 
#' @param var the continuous variable as numeric vector
#' @param ylab the label for the y axis
#'
#' @return a ggplot object with the boxplot
#'
STEPPBoxPlot <- function(var, ylab) {

    ## get data frame for plotting
    data.plot <- data.frame(var)

    ## create main plot
    pl <- ggplot(data.plot, aes(x=factor(1), y=var)) + geom_boxplot() + xlab("") + ylab(ylab)

    ## return plot
    return(pl)
}

#'
#' Helper function to produce a density plot of the continuous variable
#' 
#' @param var the continuous variable as numeric vector
#' @param xlab the label for the x axis
#'
#' @return a ggplot object with the density plot
#'
STEPPDensityPlot <- function(var, xlab) {

    ## get data frame for plotting
    data.plot <- data.frame(var)

    ## create main plot
    pl <- ggplot(data.plot, aes(x=var)) + geom_histogram(aes(y= ..density..), alpha=0.5) + geom_line(stat="density", color="blue") + xlab(xlab)

    ## return plot
    return(pl)
}



#'
#' Helper function to combine 
#'
#' @param pl1 plot number 1 (above)
#' @param pl2 plot number 2 (below)
#'
#' @return combined plot with aligned x axis
#'
STEPPCombinePlots <- function(pl1, pl2) {

    ## create gtables from the objects
    g1 <- ggplot_gtable(ggplot_build(pl1))
    g2 <- ggplot_gtable(ggplot_build(pl2))

    ## get maximal width
    width.max = unit.pmax(g1$widths[2:3], g2$widths[2:3])

    ## set maximal with for all objects
    g1$widths[2:3] <- width.max
    g2$widths[2:3] <- width.max

    ## arrange objects
    pl <- arrangeGrob(g1, g2, nrow=2, ncol=1, clip=T, heights = unit(c(2,.50),c("null", "null")))
    return(pl)
}


#'
#' Plot function for an STEPP object
#'
#' @param x STEPP object as returned by STEPP function
#' @param type the type of the plot to be produced: can be 'coef' for just the coefficients, 'group' for just the groupings, 'coef_group' for both plots aligned, or 'coef_box' and 'coef_density' for the coefficient plot with either a boxplot or density plot displaying the distribution of the variable of interest
#' @param coef.xaxis type of the x-axis of the coefficient plot, either 'medians' for a continuous scale using the medians of the variable of interest in each subgroup. 'groups' uses a discrete (and equidistant) x-axis with group names (limits of the variable of interest in that subgroup). The parameter has only an effect when type is 'coef', for combined plots it is set automatically.
#' @param coef.ylab the y-axis label for the coefficient plot. It missing the label will be guessed based on the type of the analysis: cox, logistic, survrate or rate
#' @param plot should the resulting plot be plotted or only returned ?
#' @param ... further arguments to  STEPPCoefPlot or STEPPGroupPlot if type is 'group'
#'
#' @export
#'
plot.STEPP <- function(x, type=c("coef", "group", "coef_group", "coef_box", "coef_density"), coef.xaxis=c("medians", "groups"), coef.ylab, plot=T, ...) {

    ## match argumnets
    type <- match.arg(type)
    coef.xaxis <- match.arg(coef.xaxis)

    ## coef plot if needed
    if (type %in% c("coef", "coef_group", "coef_box", "coef_density")) {
       
        ## set coef.xaxis right
        coef.xaxis <- switch(type,
                             coef=coef.xaxis,
                             coef_group="groups",
                             coef_box="medians",
                             coef_density="medians")

        ## get groupwise medians
        ## and xlab
        if (coef.xaxis=="medians") {
            groups <- sapply(x$groups, function(xx) median(x$model[xx,x$var],na.rm=T))
            xlab <- paste(x$var, "(groupwise median)")
        } else {
            groups <- NULL
            xlab <- x$var
        }
    
        ## get ylab
        if (missing(coef.ylab)) {
            coef.ylab <- switch(x$type,
                           cox="HR",
                           survrate="Survival",
                           logistic="OR",
                           rate="Rate")
            coef.ylab <- paste(ifelse(x$log, "log2 ", ""), coef.ylab, sep="")
        }

        
        ## get plot
        pl.coef <- STEPPCoefPlot(x$coefs, groups=groups, log=x$log, xlab=xlab, ylab=coef.ylab, ...)
    }

    ## add group plot if needed
    ## if only group plot is demanded pass ... argument to STEPPGroupPlot
    if (type=="group") {

        ## get intervals, better to calculate them than to use rownames of coefs matrix (would be possible only for sliding window anyway)
        intervals <- t(sapply(x$groups, function(xx) range(x$model[xx,x$var],na.rm=T)))

        ## get plot
        pl.gr <- STEPPGroupPlot(groups = x$groups, intervals=intervals, ylab = x$var, ...)

        ## adjust margins
        pl.gr <- pl.gr + theme(plot.margin = unit(c(1,1,1.5,0.5), "lines"))

    } else if (type=="coef_group") {

        ## get intervals, better to calculate them than to use rownames of coefs matrix (would be possible only for sliding window anyway)
        intervals <- t(sapply(x$groups, function(xx) range(x$model[xx,x$var],na.rm=T)))

        ## get plot
        pl.gr <- STEPPGroupPlot(groups = x$groups, intervals=intervals, ylab = x$var)
    }
    
    ## create boxplot or density plot if needed
    if (type=="coef_box") {
        
        ## get plot
        pl.bx <- STEPPBoxPlot(x$model[[x$var]], ylab = x$var) 

        ## flip coordinates and
        ## get same x axis as coef plot
        xlimits <- ggplot_build(pl.coef)$panel$ranges[[1]]$x.range
        pl.bx <- pl.bx + coord_flip(ylim=xlimits)

        ## get rid of axes
        pl.bx <- pl.bx + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) 
    } else if (type=="coef_density") {
        
        ## get plot
        pl.dens <- STEPPDensityPlot(var = x$model[[x$var]], xlab = x$var)

        ## get same x axis as coef plot
        xlimits <- ggplot_build(pl.coef)$panel$ranges[[1]]$x.range
        pl.dens <- pl.dens + coord_cartesian(xlim=xlimits)
            
        ## delete x axis
        pl.dens <- pl.dens + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
    }
    

    ## combine plots
    pl <- switch(type,
                 coef=pl.coef,
                 group=pl.gr,
                 coef_group=STEPPCombinePlots(pl.coef, pl.gr),
                 coef_box=STEPPCombinePlots(pl.coef, pl.bx),
                 coef_density=STEPPCombinePlots(pl.coef, pl.dens))

    ## get return list: has also single plots
    ret <- switch(type,
                 coef=list(pl=pl),
                 group=list(pl=pl),
                 coef_group=list(pl=pl, pl.coef=pl, pl.gr=pl.gr),
                 coef_box=list(pl=pl, pl.coef=pl.coef, pl.bx=pl.bx),
                 coef_density=list(pl=pl, pl.coef=pl.coef, pl.dens=pl.dens))

    ## plot the combined plot
    if (plot) {
        print(pl)
    }

    ## return combined plot
    invisible(ret)
}

#'
#' Bonferroni adaption of the confidence interval
#' 
#' @param coefs matrix with coefficients and confidence intervals. The first column holds the coeffecients, the second the lower and the third the upper bound of the interval
#' @param alpha the \alpha from the CI
#' @param k number of subgroups. Default to number of rows from the coefs matrix
#'
AdaptCI <- function(coefs, alpha, k=nrow(coefs)) {

    ## get adaption factor
    mu <- qnorm(1-alpha/(2*k))/qnorm(1-alpha/2)

    ## adapt CI
    coefs[,2] <- coefs[,1]-(coefs[,1]-coefs[,2])*mu
    coefs[,3] <- coefs[,1]+(coefs[,3]-coefs[,1])*mu

    ## return results
    return(coefs)

}


#'
#' ComputeTreatSigma
#'
#' Helper function to calculate the sigma matrix for the treatment variable and the different subgroups
#'
#' @param treat 
#' @param groups
#' 
ComputeTreatSigma <- function(treat, groups) {

    ## get treatment variable vector
    t.mat <- matrix(rep((as.integer(treat)-1), length(groups)), ncol=length(groups))
    colnames(names(groups))
    t.mat <- sapply(1:length(groups), function(xx) {t.mat[!(1:nrow(t.mat)) %in% groups[[xx]],xx] <- NA; t.mat[,xx]})

    ## get sigma matrix
    sigma <- cov(t.mat, use="pairwise.complete.obs")

    ## replace NAs with 0
    sigma[is.na(sigma)] <- 0

    ## return sigma
    return(sigma)

}


#'
#' Helper function to peform an Omnibus test for the equality of coefficients across subpopulations (groups)
#'
#' @param coefs matrix with coefficients and variances
#' @param treat the treatment variable as two level factor
#' @param groups the list with groups
#'
#' @return list with 
#'
STEPPOmnibusTest <- function(coefs, treat, groups) {

    ## get number of groups
    k <- length(groups)

    ## get difference of coefs delta
    delta <- diff(coefs[,"betas"])

    ## get variances
    A=matrix(c(-1, 1, rep(0, (nrow(coefs)-1))), 
             nrow=(nrow(coefs)-1), 
             ncol=nrow(coefs), byrow=T)

    ## get sigma matrix
    sigma <- ComputeTreatSigma(treat, groups)

    ## invert product of sigma and A
    sigma.inv <- solve(A %*% sigma %*% t(A))

    ## compute G
    G <- t(delta) %*% sigma.inv %*% delta

    ## compute p value based on chisq with k-1 degrees of freedom
    return(list(G=as.numeric(G), df=(k-1), p.value=pchisq(as.numeric(G), df=(k-1), lower.tail=F)))
 
}



#'
#' STEPPSurvModel
#'
#' A helper function for the STEPP analysis in case of a surival response: if a treatment is given a Cox model is used with optional covariates. If no treatment is given a survival rate (or more precisely the survival probability) is estimated using the Kaplan-Meier estimate at a given timepoint
#'
#' @param response the response, must be a Surv object
#' @param treat the name of the treatment variable in the data
#' @param covars optinal names of covariates
#' @param data the data frame used for evalution
#' @param the timepoint for evaluating the KM estimate, only used if no treatment is given
#' @param alpha the alpha level used for confidence intervals
#'
#' @import survival
#'
STEPPSurvModel <- function(response, treat, covars=NULL, data, timepoint=60, alpha) {

    ## prepare vector with return values
    ret <- rep(NA, 3)

    if (!missing(treat)) {
        ## fit cox model if treat is give

        formula <- paste("response", treat, sep="~")

        ## ceck for covariates if given and create formula
        if(!is.null(covars)) {
            formula <- paste(formula, paste(covars, collapse="+"), sep="+")
        }
        formula <- as.formula(formula)

        ## fit cox model
        cox <- coxph(formula, data=data)
            
        # extract hazard ratio
        ret[1] <- coef(cox)[1]

        # get confint
        ci <- confint(cox, level=(1-alpha))
        ret[2:3] <- ci[1:2]

    }  else {

        ## if no treat is given ignore covars
        ## just uses KM estimate of survival at given timepoint
        km <- survfit(response~1, conf.int=(1-alpha))
        km.sum <- summary(km, time=timepoint, extend=T)

        ## extract coefs -> here the survival
        ret[1] <- km.sum$surv
        ret[2] <- km.sum$lower
        ret[3] <- km.sum$upper
    }

    return(ret)

}

#'
#' STEPPLogisticModel
#'
#' A helper function for the STEPP analysis in case of a binary response: fits a logistic regression model if a treatment and optional covars are given. If no treatment is given a rate is calculated. In the latter case confidence intervals are calculates based on the biomial distribution.
#'
#' @param response the name of the response in the data
#' @param treat the name of the treatment variable in the data
#' @param covars optinal names of covariates
#' @param data the data frame used for evalution
#' @param alpha the alpha level used for confidence intervals
#' 
STEPPLogisticModel <- function(response, treat, covars=NULL, data, alpha, verbose=T) {

    ## prepare vector with return values
    ret <- rep(NA, 3)
    
    if (!missing(treat)) {
        ## Logistic Regression
        if (verbose) cat("Treatment ", treat, ", performing logistic regression with ", nrow(data), " patients at alpha ", alpha, "\n")

        formula <- paste(response, treat, sep="~")

        ## ceck for covariates if given and create formula
        if (!is.null(covars)) {
            formula <- paste(formula, paste(covars, collapse="+"), sep="+")
        }
        formula <- as.formula(formula)

        # fit model
        lrm <- glm(formula, data=data, family=binomial())


        # extract odds ratio
        # note, first comes intercept
        ret[1] <- coef(lrm)[2]

        # get confint
        # again, first comes intercept
        ci <- try(confint(lrm, level=(1-alpha)))
        if (!inherits(ci, "try-error")) {
            ret[2] <- ci[2,1]
            ret[3] <- ci[2,2]
        }     

    } else {
            
        ## we have only a simple rate, use this for a binom test with CI
        # number of sucesses, we define second level of the vector as success
        s <- table(data[[response]])[2]
        # number of tries
        n <- sum(!is.na(data[[response]]))

        if (verbose) cat("No treatment, performing binomial test with", n, "patients at alpha", alpha, "\n")

        ## do exact binomial test
        bin.test <- binom.test(s, n, conf.level=(1-alpha))

        ret[1] <- s/n
        ret[2] <- bin.test$conf.int[1]
        ret[3] <- bin.test$conf.int[2]

    }

    return(ret)
}


#'
#' STEPP analysis
#'
#' @param response the response of the model, either a Surv object or the name of a factor variable in 'data'
#' @param var the name of the continuous variable
#' @param treat the name of the treatment variable, either a numerical vector or a factor. First level is fail, second is success
#' @param covars other variables to include in the model. Default to NULL
#' @param data the data frame the variables are evaluated in
#' @param mode the method to compute subrgoups. Either 'sw' for sliding-window or 'to' for tail-oriented
#' @param alpha the alpha for testing. Default to 0.05
#' @param adapt.ci indicates whether a pointwise confidence iterval is drawn or an confidence band
#' @param log indicates wheter the Odds-ratios/Hazard ratios will be given at a log scale or not. Success rates (binary outcome) and survival rates are not affected. Default to FALSE
#' @param ... further arguments to 'STEPPSliding' or 'STEPPTail' depending on 'mode'
#'
#' @return object from class 'STEPP'
#'
#' @export
#'
STEPP <- function(response, var,  treat, covars=NULL, data=parent.frame(), mode=c("sw", "to"), alpha=0.05, adapt.ci=c("none", "normal", "bonferroni"), log=F, timepoint,  ...) {

    
    ## check for var, and treatment
    if (!missing(treat)) {
        args <- c(var=var, treat=treat) 
    } else { 
        args <- c(var=var)
    }
    m <- match(args, names(data))
    if(any(is.na(m))) {
        stop(paste(names(args)[is.na(m)], collapse=","), "cannot be found in the data")
    }

    ## check response
    if (!is.Surv(response) && is.character(response)) {
        if (!response %in% names(data)) {
            stop("response is not given in the data")
        }
    } else if (!is.Surv(response) && !is.character(resonse)) {
        stop("Response must be either a Surv object or the name of a factor in data")
    }

    ## get type: cox, survrate, logistic, rate
    if(is.Surv(response)) { 
        type <- if (!missing(treat)) "cox" else "survrate"
    } else {
        type <- if (!missing(treat)) "logistic" else "rate"
    }

    ## check if treat is numeric of factor with two levels
    if (!missing(treat) && !is.numeric(data[[treat]]) && (is.factor(data[[treat]]) && nlevels(data[[treat]])!=2)) {
        warning("treat is neither numeric or a factor with two levels.")
    }

    ## check if x in data is numeric
    if (!is.numeric(data[[var]])) {
        warning(var, "is converted to numeric!")
        data[[var]] <- as.numeric(as.character(data[[var]]))
    }
    if (length(unique(data[[var]])) < 10) {
        warning("x does not seem to be a continual variable!")
    }

    
    ## check mode
    mode <- match.arg(mode)

    ## check adapt.ci
    adapt.ci <- match.arg(adapt.ci)

    if (mode=="sw") {

        ## get groups
        groups <- STEPPSliding(data[[var]], ...)

        ## get number of groups
        k <- length(groups)
    } else if (mode=="to") {

        ## get groups
        ttt <- STEPPTail(data[[var]], ...)
        groups.up <- ttt$groups.up
        groups.down <- ttt$groups.down
        groups <- c(groups.up, groups.down)

        ## get number of groups
        k <- length(groups)
    }

    ## adapt alpha if needed
    if (adapt.ci=="bonferroni") {
        alpha <- alpha/k
    }
    
    ## allocate three vectors for coefficients and the confidence bounds
    coefs <- matrix(NA, nrow=k, ncol=3, dimnames=list(names(groups), c("betas", "lower", "upper")))
    ## main loop over subgroups
    for (b in 1:k) {

        data.b <- data[groups[[b]],]

        coefs[b,] <- switch(type, 
                            cox=STEPPSurvModel(response=response[groups[[bb]],], treat=treat, covars=covars, data=data.b, alpha=alpha),
                            survrate=STEPPSurvModel(response=response[groups[[b]],], timepoint=timepoint, alpha=alpha),
                            logistic=STEPPLogisticModel(response=response, treat=treat, covars=covars, data=data.b, alpha=alpha),
                            rate=STEPPLogisticModel(response=response, data=data.b, alpha=alpha)
                            )

    }

    ## adapt CI if needed
    ## use Bonferoni correction for CIs
    if (adapt.ci=="normal") {
        coefs <- AdaptCI(coefs, alpha=alpha, k=k)
    }


    ## test if type is "cox" or "logistic"
    if (type %in% c("cox", "logistic")) {

        test <- list()

        if (mode=="sw") {

            ## perform one overall test
            test$overall <-  STEPPOmnibusTest(coefs, data[[treat]], groups)

        } else if (mode=="to") {

            ## perform overall tests
            ## test both sides separately
            ## on log coefficients
            test$up <- STEPPOmnibusTest(coefs[1:((k+1)/2),], data[[treat]], groups.up)
            test$down <- STEPPOmnibusTest(coefs[((k+3)/2):nrow(coefs),], data[[treat]], groups.down)

        }
    } else {
        test <- NULL
    }

    ## in case of HR or OR (type 'cox' or 'logistic') 
    ## check log parameter
    if (type %in% c("cox", "logistic") && !log) {
        coefs <- exp(coefs)
    }


    ## return result
    ret <- list(coefs=coefs, groups=groups, test=test, model=data, var=var, mode=mode, type=type, log=log)
    attr(ret, "class") <- "STEPP"
    return(ret)

}

