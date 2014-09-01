#
# This is a generalized but slow replacement for limma
#
# Running multi-core may require running R with:
#
#     OPENBLAS_NUM_THREADS=1 R
#
#

library(Matrix)
library(limma)  # EList class, p.adjust
library(parallel)

##################################
# Utility functions
##################################

optim.positive <- function(initial, func) {
    result <- optim(rep(0,length(initial)), function(x) func(exp(x)*initial))
    result$par <- exp(result$par) * initial
    result
}


as.diagonal.matrix <- function(vec) diag(vec, length(vec))

invert.matrix <- function(A) {
    if (length(A)) 
        solve(A)
    else
        A #Empty matrix. As usual R+ fails to handle an edge case correctly.
}

##################################
# Multivariate distribution classes
##################################


good <- function(x, ...) UseMethod('good') #Which dimensions are usable (non-Inf covariance)
expect <- function(x, ...) UseMethod('expect')
p.value <- function(x, ...) UseMethod('p.value')
log.density <- function(x, ...) UseMethod('log.density')
random <- function(x, ...) UseMethod('random')
transformed <- function(x, ...) UseMethod('transformed')
marginal <- function(x, ...) UseMethod('marginal')
conditional <- function(x, ...) UseMethod('conditional')

weights.vector <- function(x, ...) UseMethod('weights.vector')

#
# Sample p values, should be a uniform distribution and plot a straight line
#
verify.p.values <- function(dist) {
    data <- mapply(function(i) p.value(dist, random(dist)), 1:1000)
    plot(sort(data))
}


##################################
# Multivariate Normal
##################################

mvnormal <- function(mean,covar) {
    object <- list(
        mean=as.matrix(mean),
        covar=as.matrix(covar)
        )
    class(object) <- 'mvnormal'
    object    
}

good.mvnormal <- function(self) {
    good <- is.finite(self$mean)
    for(i in seq_len(nrow(self$covar))) {
        good <- good & is.finite(self$covar[i,])
    }
    good
}

expect.mvnormal <- function(self) c(self$mean)

weights.vector.mvnormal <- function(self) {
    1.0 / diag(self$covar)
}

random.mvnormal <- function(self) {
    A <- chol(self$covar)
    c(self$mean + t(A) %*% rnorm(ncol(A)))
}

p.value.mvnormal <- function(self, x) {
   offset <- as.matrix(x) - self$mean
   df <- nrow(self$covar)
   q <- t(offset) %*% solve(self$covar, offset)
   pchisq(c(q), df, lower.tail=FALSE)
}

log.density.mvnormal <- function(self, x) {
   offset <- as.matrix(x) - self$mean
   
   -0.5*( 
     log(2*pi)*nrow(self$covar)
     + log(det(self$covar))
     + t(offset) %*% solve(self$covar, offset)
   )
}

transformed.mvnormal <- function(self, A, offset=c(0)) {
    mvnormal(       
        A %*% self$mean + offset,
        A %*% self$covar %*% t(A)
    )
}

marginal.mvnormal <- function(self, i) {
    mvnormal(
        self$mean[i],
        self$covar[i,i]
    )
}

conditional.mvnormal <- function(self, i1, i2, x2) {
    x2 <- as.matrix(x2)
    
    mean1 <- self$mean[i1]
    mean2 <- self$mean[i2]
    offset2 <- x2-mean2
    covar11 <- self$covar[i1,i1,drop=FALSE]
    covar12 <- self$covar[i1,i2,drop=FALSE]
    covar21 <- self$covar[i2,i1,drop=FALSE]
    covar22 <- self$covar[i2,i2,drop=FALSE]
    covar22inv <- invert.matrix(covar22)
    covar12xcovar22inv <- covar12 %*% covar22inv
    
    mvnormal(
        mean1 + covar12xcovar22inv %*% offset2,
        covar11 - covar12xcovar22inv %*% covar21
    )
}


##################################
# Multivariate t
##################################

mvt <- function(mean,covar,df) {
    object <- list(
        mean=as.matrix(mean),
        covar=as.matrix(covar),
        df=df
        )
    class(object) <- 'mvt'
    object    
}

good.mvt <- function(self) {
    good <- is.finite(self$mean)
    for(i in seq_len(nrow(self$covar))) {
        good <- good & is.finite(self$covar[i,])
    }
    good
}

expect.mvt <- function(self) c(self$mean)

weights.vector.mvt <- function(self) {
    1.0 / diag(self$covar)
}

random.mvt <- function(self) {
    A <- chol(self$covar)
    c( self$mean + (t(A)%*%rnorm(ncol(A))) * sqrt(self$df/rchisq(1,df=self$df)) )
}

p.value.mvt <- function(self, x) {
   offset <- as.matrix(x) - self$mean
   p <- nrow(self$covar)
   q <- (t(offset) %*% solve(self$covar, offset)) / p
   pf(c(q), p, self$df, lower.tail=FALSE)
}

log.density.mvt <- function(self, x) {
   offset <- as.matrix(x) - self$mean
   p <- nrow(self$covar)
   v <- self$df
   (
      lgamma(0.5*(v+p))
      - lgamma(0.5*v) 
      - (0.5*p)*log(pi*v)
      - 0.5*log(det(self$covar))
      - (0.5*(v+p))*log(1+c( t(offset) %*% solve(self$covar, offset) )/v)
   )
}

transformed.mvt <- function(self, A, offset=c(0)) {
    mvt(       
        A %*% self$mean + offset,
        A %*% self$covar %*% t(A),
        self$df
    )
}

marginal.mvt <- function(self, i) {
    mvt(
        self$mean[i],
        self$covar[i,i],
        self$df
    )
}

conditional.mvt <- function(self, i1, i2, x2) {
    x2 <- as.matrix(x2)
    p2 <- length(i2)
    
    mean1 <- self$mean[i1]
    mean2 <- self$mean[i2]
    offset2 <- x2-mean2
    covar11 <- self$covar[i1,i1,drop=FALSE]
    covar12 <- self$covar[i1,i2,drop=FALSE]
    covar21 <- self$covar[i2,i1,drop=FALSE]
    covar22 <- self$covar[i2,i2,drop=FALSE]
    covar22inv <- invert.matrix(covar22)
    covar12xcovar22inv <- covar12 %*% covar22inv
    df <- self$df
    
    mvt(
        mean1 + covar12xcovar22inv %*% offset2,
        
        (covar11 - covar12xcovar22inv %*% covar21) *
            ((df + c(t(offset2)%*%covar22inv%*%offset2)) / (df + p2)),
        
        df + p2
    )
}


##################################
# Core routines
##################################


#
# Input:
#
# get.dist(i, param)
# - i is feature number (row number of data)
# - param is a parameter vector
# - optimization is constrained to all param positive,
#   starting from initial
#
# Output:
#
# Some dist may not be all(good( ))
#
fit.noise <- function(data, design, get.dist, initial, cores=1) {
    #decomp <- qr(design)
    #i2 <- ncol(design)+seq_len(nrow(design)-ncol(design))
    #Q <- qr.Q(decomp, complete=TRUE)
    #Q2 <- Q[, i2, drop=FALSE]
    #tQ2 <- t(Q2)
    #
    #z2 <- data %*% Q2
    #
    #is.good <- rep(TRUE, nrow(data))
    #for(i in seq_len(nrow(data)))
    #    if (!all(is.finite(data[i,])) ||
    #        !all(good(get.dist(i, initial))))
    #        is.good[i] <- FALSE    
    #
    #scorer <- function(param) {
    #    total <- 0.0
    #    for(i in seq_len(nrow(data))) 
    #        if (is.good[i]) {
    #            dist <- transformed(get.dist(i, param), tQ2)
    #            total <- total + log.density(dist, z2[i,])
    #        }
    #    -total
    #}     
    #
    #result <- optim.positive(initial, scorer)
        
    indicies <- c()
    retains <- list()
    z2s <- list()
    tQ2s <- list()    

    for(i in seq_len(nrow(data))) {
        retain <- seq_len(ncol(data))[ is.finite(data[i,]) & good(get.dist(i, initial)) ]
        
        # Does enough data remain?
        if (length(retain) > ncol(design)) {
            decomp <- qr(design[retain,,drop=FALSE])
            i2 <- ncol(design)+seq_len(length(retain)-ncol(design))
            Q <- qr.Q(decomp, complete=TRUE)
            tQ2 <- t(Q[, i2, drop=FALSE])            
            z2 <- tQ2 %*% data[i,retain]
        
            indicies[length(indicies)+1] <- i
            retains[[i]] <- retain
            z2s[[i]] <- z2
            tQ2s[[i]] <- tQ2
        }
    }

    if (length(initial) == 0) {
        result <- list(par=initial)
    } else {
        if (cores <= 1) {
            scorer <- function(param) {
                total <- 0.0
                for(i in indicies) {
                    dist <- transformed(
                        marginal(
                            get.dist(i, param), 
                            retains[[i]]), 
                        tQ2s[[i]])
                    total <- total + log.density(dist, z2s[[i]])
                }
                -total
            }     
        } else {
            scorer <- function(param) {
                score.one <- function(i) {
                    dist <- transformed(
                        marginal(
                            get.dist(i, param), 
                            retains[[i]]), 
                        tQ2s[[i]]
                        )
                    log.density(dist, z2s[[i]])
                }
                
                -sum(unlist(mclapply(indicies, score.one, mc.cores=cores)))
            }        
        }

        result <- optim.positive(initial, scorer)
    }    
    
    dist <- lapply(seq_len(nrow(data)), function(i) get.dist(i, result$par))

    
    noise.p.values <- rep(NA,nrow(data))
    for(i in indicies) {
        this.dist <- transformed(
            marginal(dist[[i]], retains[[i]]), 
            tQ2s[[i]]
            )
        noise.p.values[i] <- p.value(this.dist, z2s[[i]])
    }
    
    good.noise.p.values <- noise.p.values[!is.na(noise.p.values)]
    
    #Bonferroni combined p value
    noise.combined.p.value <- min(1, min(good.noise.p.values) * length(good.noise.p.values))
    
    #Stouffer combined p value
    #noise.combined.p.value <- pnorm(sum(qnorm(good.noise.p.values)) / sqrt(length(good.noise.p.values)), lower.tail=T)
    
    
    noise.description <- sprintf('noise p-value = %g  (small is bad)', noise.combined.p.value)

    
    object <- list(
        design = design,
        optim.output = result,
        param = result$par,
        data = data,
        dist = dist,
        noise.p.values = noise.p.values,
        noise.combined.p.value = noise.combined.p.value,
        noise.description = noise.description
    )
    class(object) <- 'fitnoise'
    
    object
}

#
# Input:
#
# fit - result from fit.noise
#
# Output:
#
# Some coef.dist may be NULL
#
fit.coef <- function(fit, design) {
    if (is.null(colnames(design))) {
        coef.names <- character()
        for(i in seq_len(ncol(design)))
            coef.names[i] <- sprintf('coef%d', i)
        colnames(design) <- coef.names
    }

    data <- fit$data
    dist <- fit$dist

    coef.dist <- list()
    length(coef.dist) <- nrow(data)
    
    coef <- matrix(as.double(NA),nrow=nrow(data),ncol=ncol(design))
    colnames(coef) <- colnames(design)
    
    for(i in seq_len(nrow(data))) {        
        retain <- seq_len(ncol(data))[ is.finite(data[i,]) & good(dist[[i]]) ]
        
        # Does enough data remain?
        if (length(retain) >= ncol(design) &&
            ncol(design) > 0 &&
            rankMatrix(design[retain,,drop=FALSE]) >= ncol(design)) {
            
            this.dist <- marginal(dist[[i]], retain)
            
            decomp <- qr(design[retain,,drop=FALSE])    
            i1 <- seq_len(ncol(design))
            i2 <- ncol(design)+seq_len(length(retain)-ncol(design))
            Q <- qr.Q(decomp, complete=TRUE)
            R <- qr.R(decomp)
            Rinv <- invert.matrix(R)
            tQ <- t(Q)
            Q2 <- Q[, i2, drop=FALSE]
            tQ2 <- t(Q2)
            
            z <- tQ %*% data[i,retain]
            
            cond.dist <- conditional(transformed(this.dist,tQ), i1, i2, z[i2])
            coef.dist[[i]] <- transformed(cond.dist, -Rinv, Rinv %*% z[i1])
            coef[i,] <- expect(coef.dist[[i]])
        }
    }

    fit$coef <- coef
    fit$coef.dist <- coef.dist
    
    fit
}



fit <- function(data, design, get.dist, initial, noise.design=NULL, cores=1) {
    if (is.null(noise.design))
        noise.design <- design
    
    fit.coef(fit.noise(data, noise.design, get.dist, initial, cores=cores), design)
}



print.fitnoise <- function(self) {
    if (!is.null(self$noise.description))
        cat(self$noise.description,'\n')  
    else
        cat('Fit',self$param,'\n')
}


#
# Output:
#
# Some pvalues may be NA,
# these shouldn't be counted toward multiple testing.
#
p.values.contrasts <- function(fit, contrasts) {
    contrasts <- as.matrix(contrasts)

    p.values <- rep(as.double(NA), nrow(fit$coef))
    for(i in seq_len(length(p.values)))
        if (!is.null(fit$coef.dist[[i]])) {
            p.values[i] <- p.value( 
                transformed( fit$coef.dist[[i]], t(contrasts) ), 
                rep(0, ncol(contrasts))
                )
        }
    p.values
}

p.values.coefs <- function(fit, coefs) {
    contrasts <- matrix(0, nrow=ncol(fit$coef), ncol=length(coefs))
    for(i in seq_len(length(coefs)))
        contrasts[coefs[i],i] <- 1

    p.values.contrasts(fit, contrasts)
}

p.values.noise <- function(fit) {
    fit$noise.p.values
}


weights.matrix <- function(fit) {
    weights <- matrix(0, nrow=nrow(fit$data), ncol=ncol(fit$data))

    for(i in seq_len(nrow(fit$data)))
        if (!is.null(fit$dist[[i]]))
            weights[i,] <- weights.vector(fit$dist[[i]])
    
    weights
}


##################################
# Noise models for EList objects
##################################

normal.model.to.t <- function(model) list(
    get.dist = function(elist, i, param) {
        normal.dist <- model$get.dist(elist, i, param[-1])
        mvt(normal.dist$mean, normal.dist$covar, param[1])
    },
    
    initial = function(elist) c(1.0, model$initial(elist)),
    
    describe = function(elist, param) {
        normal.desc <- model$describe(elist, param[-1])
        sprintf('prior df = %g\n%s', param[1], normal.desc)
    }
)

model.normal.standard <- list(
    # elist$weights if present is taken as proportional to 1.0/variance

    get.dist = function(elist, i, param) {
        if (is.null(elist$weights))
            var <- rep(param[1], ncol(elist))
        else
            var <- param[1] / elist$weights[i,]
        
        mvnormal(rep(0,ncol(elist)), as.diagonal.matrix(var))
    },

    initial = function(elist) c(1.0),
    
    describe = function(elist, param) {
        if (is.null(elist$weights))
            sprintf('variance = %g', param[1])
        else
            sprintf('variance = %g / weight', param[1])
    }
)

model.t.standard <- normal.model.to.t(model.normal.standard)


model.t.independent <- list(
    get.dist = function(elist, i, param) {
        if (is.null(elist$weights))
            var <- rep(1.0, ncol(elist))
        else
            var <- 1.0 / elist$weights[i,]
        mvt(rep(0,ncol(elist)), as.diagonal.matrix(var), 0.0)
    },
    
    initial = function(elist) { double() },
    
    describe = function(elist, param) { 'Independent t-tests.' }
)


model.normal.patseq <- list(
    # elist$E is average read tail length
    # elist$other$counts is polya read count

    get.dist = function(elist, i, param) {
        var <- param[1] / elist$other$counts[i,] + (elist$E[i,]**2) * param[2]
        mvnormal(rep(0,ncol(elist)), as.diagonal.matrix(var))
    },
    
    initial = function(elist) { c(250.0, 0.001) },
    
    describe = function(elist, param) {
        sprintf('variance = %.1f / polya_read_count + %.6f * tail_length^2',param[1],param[2])
    }
)

model.t.patseq <- normal.model.to.t(model.normal.patseq)


model.normal.quadratic <- list(
    # variance as a quadratic function of elist$other$x
    # Note: setting elist$other$x = elist$E will result in overfitting,
    #       in particular the constant term will be too small
    
    get.dist = function(elist, i, param) {
        expr <- elist$other$x[i,]
        var <- pmax(1e-10, expr*expr*param[1] + expr*param[2] + param[3] )
        mvnormal(rep(0,ncol(elist)), as.diagonal.matrix(var))
    },
    
    initial = function(elist) { c(1.0,1.0,1.0) },
    
    describe = function(elist, param) {
        sprintf('variance = %g x^2 + %g x + %g', param[1],param[2],param[3])
    }
)

model.t.quadratic <- normal.model.to.t(model.normal.quadratic)


##################################
# Convenience functions for EList objects
##################################


fit.elist <- function(
        elist, design, model=model.t.standard, 
        noise.design=NULL,
        cores=1) {
    result <- fit(
        data = elist$E,
        design = design,
        get.dist = function(i,param) model$get.dist(elist,i,param),
        initial = model$initial(elist),
        noise.design = noise.design,
        cores = cores
    )    
    result$elist <- elist    
    result$noise.description <- sprintf(
        "%s\n%s", 
        model$describe(elist, result$param), 
        result$noise.description
        )
    
    result
}


#
# Input:
#
# Specify either coefs or contrasts or neither.
# If neither is given, this tests the fit of the noise!
#
# Output:
#
# Intended to match the output of limma's topTableF
#
test.fit <- function(fit, coefs=NULL, contrasts=NULL, sort=TRUE) {
    stopifnot(is.null(coefs) + is.null(contrasts) <= 1)
    
    if (!is.null(coefs))
        p.values <- p.values.coefs(fit, coefs)
    else if (!is.null(contrasts))
        p.values <- p.values.contrasts(fit, contrasts)
    else
        p.values <- fit$noise.p.values
    
    fdrs <- p.adjust(p.values, method='fdr')
    
    # Handle weighting correctly
    average.fit <- fit.coef(fit, cbind(rep(1,ncol(fit$data))))


    table <- data.frame(matrix(nrow=length(p.values),ncol=0))
    
    for(i in seq_len(ncol(fit$coef)))
        table[,colnames(fit$coef)[i]] <- fit$coef[,i]
    
    table[,'AveExpr'] <- average.fit$coef[,1]

    table[,'P.Value'] <- p.values
    table[,'adj.P.Val'] <- fdrs
    
    table[,'noise.P.Value'] <- fit$noise.p.values
    
    if (!is.null(fit$elist) && !is.null(fit$elist$genes)) {
        genes <- fit$elist$genes
        if (!is.null(rownames(genes)))
            rownames(table) <- rownames(fit$elist$genes)
        
        for(i in seq_len(ncol(genes)))
            table[,colnames(genes)[i]] <- genes[,i]
    }
    
    if (sort)
        table <- table[order(p.values),]
    
    table
}





