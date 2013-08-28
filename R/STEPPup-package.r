#' STEPPup
#'
#' @name STEPPup
#' @docType package
NULL


#'
#' Helper function for parallel processing
#'
CutPointCV <- function(c, response, var, treat, covars, data, K, repeats, seed) {

    ## create dichotomized variable
    data.c <- data
    data.c$varc <- factor(ifelse(data[[var]]<=c, "low", "high"), levels=c("low", "high"))

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

            ## get formula
            ## 1) Logistic Regression
            if (!is.Surv(response)) {
                formula <- paste(response, paste("varc", treat, sep="+"), sep="~") 
            } else {
                formula <- paste("response~varc", treat, sep="+")
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
        ComputeAUC(predictions, labels=data[[response]])
    }

}


#'
#' Cut point search via cross-validation
#' 
#' Searches an optimal cut-point for a continuous variable 'var' in interaction with treatmen treat. The cut point is optimized via cross-validated AUC
#'
#' @param response
#' @param var
#' @param treat
#' @param covars 
#' @param data
#' @param num.cuts
#' @param K
#' @param repeats
#' @param cores
#' @param
#'
#' @return matrix with cut points, median AUC and AUCs for every repeat
#'
#' @import parallel
#' @export
#'
OptimalCutPoint <- function(response, var, treat, covars=NULL, data, num.cuts=10000, K=10, repeats=5, cores=1, seed=123) {

    ## libraries
    require(parallel)

    ## parallel executtion
    #if (cores > 1) {
    #    registerDoMC(cores)
    #}


    ## get sorted uniqe var
    num.uni <-  sort(unique(na.omit(data[[var]])))

    ## get cuts
    ## if number of cuts > length(var) take every var as one cut except the lastg one
    if (num.cuts >= length(num.uni)) {
        cuts <- num.uni[-length(num.uni)]
    } else { 
        cuts <- seq(min(data[[var]], na.rm=T), max(data[[var]], na.rm=T), length.out=(num.cuts+2))
        cuts <- cuts[-c(1, length(cuts))]
    }   

    ## main loop parallel excuted over all posisble cuts
    cuts.errors <- mclapply(cuts, FUN=CutPointCV, response=response, var=var, treat=treat, covars=covars, data=data, K=K, repeats=repeats, seed=seed, mc.cores=cores)
    cuts.errors <- do.call("rbind", cuts.errors)

    if(!is.Surv(response)) {
        colnames(cuts.errors) <- paste("AUC", 1:ncol(cuts.errors))
        rownames(cuts.errors) <- NULL
        cuts.errors <- data.frame(cuts=cuts, "median AUC"=apply(cuts.errors, 1, median), cuts.errors)
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
                      paste("<=", max(x[xx], na.rm=T), sep="")
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
                      paste(">=", min(x[xx], na.rm=T), sep="")
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
#' @rerferences Bonetti2000_StatInMed
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
#' @param xlab the x axis label
#' @param xlab the y axis label
#' @param log indicates wheter the coefficients are on a log scale.coef
#' @param log.scale indicates wheter the y axis should use a log scale. Will be ignored if the data are already on a log scale 
#'
#' @return ggplot object with the coefficient plot
#'
STEPPCoefPlot <- function(coefs, groups=NULL, line=c("smooth", "simple", "none"), xlab="x", ylab="OR", main="", log, log.scale=T, ylim=NULL) {

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

    ## add help line for OR/HR=1
    intercept <- ifelse(log, 0, 1)
    pl <- pl + geom_hline(yintercept=intercept, linetype=2, alpha=0.5)

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
    pl <- ggplot(data.plot, aes(x=var)) + geom_line(stat="density") + xlab(xlab)

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
#' @param ylim if given the y axis of the coef plot is limited to this window (zoom). Default is NULL
#' @param log.scale indicates wheter a log scale for the y axis of the coef plot should be used. Default is FALSE. If the coefficients are already given in a log2 scale this parameter is igrnored and a warning is issued
#'
plot.STEPP <- function(x, type=c("coef", "group", "coef_group", "coef_box", "coef_density"), coef.xaxis=c("medians", "groups"),  plot=T, ...) {

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
        ylab <- paste(ifelse(x$log, "log2 ", ""), ifelse(x$type=="cox", "HR", "OR"), sep="")

        ## get plot
        pl.coef <- STEPPCoefPlot(x$coefs, groups=groups, log=x$log, xlab=xlab, ylab=ylab, ...)
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
#' @param log indicates wheter the coefficients (OR or HR) will be given at a log scale or not. Default to FALSE
#' @param ... further arguments to 'STEPPSliding' or 'STEPPTail' depending on 'mode'
#'
#' @return object from class 'STEPP'
#'
STEPP <- function(response, var,  treat, covars=NULL, data=parent.frame(), mode=c("sw", "to"), alpha=0.05, adapt.ci=T, log=F, ...) {

    
    ## check for var, and treatment
    args <- c(var=var, treat=treat)
    m <- match(args, names(data))
    if(any(is.na(m))) {
        stop(paste(names(args)[is.na(m)], collapse=","), "cannot be found in the data")
    }

    
    ## check response
    if(!is.Surv(response) && is.character(response)) {
        if (!response %in% names(data)) {
            stop("response is not given in the data")
        }
    } else {
        stop("Response must be either a Surv object or the name of a factor in data")
    }

    ## create formula
    if (is.Surv(response)) {
        formula <- paste("response", treat , sep="~")
    } else {
        formula <- paste(response, treat , sep="~")
    }

    ## ceck for covariates if given and create formula
    if(!is.null(covars)) {
        m <- match(covars, names(data))
        if(any(is.na(m))) {
            stop("Covariates", paste(covars[is.na(m)], collapse=","), "cannot be found in the data")
        }
        formula <- paste(formula, paste(covars, collapse="+"), sep="+")
    }
    formula <- as.formula(formula)


    ## check if treat is numeric of factor with two levels
    if(!is.numeric(data[[treat]]) && (is.factor(data[[treat]]) && nlevels(data[[treat]])!=2)) {
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

    if(mode=="sw") {

        ## get groups
        groups <- STEPPSliding(data[[var]], ...)

        ## get number of groups
        k <- length(groups)
    } else if(mode=="to") {

        ## get groups
        ttt <- STEPPTail(data[[var]], ...)
        groups.up <- ttt$groups.up
        groups.down <- ttt$groups.down
        groups <- c(groups.up, groups.down)

        ## get number of groups
        k <- length(groups)
    }

    ## allocate three vectors for coefficients and the confidence bounds
    betas <- lower <- upper <- rep(NA, k)

    ## main loop over subgroups
    for (b in 1:k) {

        if (is.Surv(response)) {
            ## Cox Model

            # fit model
            cox <- coxph(formula, data=data, subset=groups[[b]])

            # extract hazard ratio
            betas[b] <- coef(cox)[1]

            # get confint
            ci <- confint(cox, level=(1-alpha))
            lower[b] <- ci[1]
            upper[b] <- ci[2]

        } else {
            ## Logistic Regression

            # fit model
            lrm <- glm(formula, data=data, subset=groups[[b]], family=binomial())

            # extract odds ratio
            # note, first comes intercept
            betas[b] <- coef(lrm)[2]

            # get confint
            # again, first comes intercept
            ci <- try(confint(lrm, level=(1-alpha)))
            if (!inherits(ci, "try-error")) {
                lower[b] <- ci[2,1]
                upper[b] <- ci[2,2]
            } else {
                lower[b] <- NA
                upper[b] <- NA
            } 
        }
    }

    ## construct matrix of coefs and ci
    coefs <- cbind(betas, lower, upper)
    rownames(coefs) <- names(groups)

    ## adapt CI if needed
    ## use Bonferoni correction for CIs
    if (adapt.ci) {
        coefs <- AdaptCI(coefs, alpha=alpha, k=k)
    }

    
    ## get type: cox or logistic
    ## for STEPP object
    type <- if(is.Surv(response)) "cox" else "logistic"

    ## perform test and return results
    if (mode=="sw") {

        ## perform test
        ## on log coefficients
        test.overall <- STEPPOmnibusTest(coefs, data[[treat]], groups)

        if(!log) {
            coefs <- exp(coefs)
        }

        ## return result
        ret <- list(coefs=coefs, groups=groups, test.overall=test.overall, model=data, var=var, mode=mode, type=type, log=log)
        attr(ret, "class") <- "STEPP"
        return(ret)

    } else if (mode=="to") {

        ## perform overall tests
        ## test both sides separately
        ## on log coefficients
        test.overall.up <- STEPPOmnibusTest(coefs[1:((k+1)/2),], data[[treat]], groups.up)
        test.overall.down <- STEPPOmnibusTest(coefs[((k+3)/2):nrow(coefs),], data[[treat]], groups.down)

        if(!log) {
            coefs <- exp(coefs)
        }

        ## return result
        ret <- list(coefs=coefs, groups=groups, test.overall.up=test.overall.up, test.overall.down=test.overall.down, model=data, var=var, mode=mode, type=type, log=log)
        attr(ret, "class") <- "STEPP"
        return(ret)

    }

}

