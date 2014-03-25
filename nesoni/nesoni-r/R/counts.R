
# Functions for dealing with table created by "nesoni count:"

read.counts <- function(filename, min.total=0, min.max=0, keep=NULL, norm.file=NULL, use.tmm=TRUE, quiet=FALSE) {    
    #data <- read.delim(filename, check.names=FALSE)
    #rownames(data) <- data$Feature
    #n_samples <- grep('^RPKM', colnames(data))[1] - 2
    #counts <- as.matrix( data[,2:(n_samples+1), drop=FALSE] )
    #gene <-data[, (n_samples*2+2):ncol(data)]
    
    data <- read.grouped.table(filename, require=c('Count'), default.group='Count')
    counts <- as.matrix( data$Count )
    
    gene <- data.frame(row.names=rownames(counts), locus_tag=rownames(counts))
    
    if (!is.null(data$Annotation)) { 
        gene <- cbind(gene, data$Annotation)
    }
    if (!is.null(data$Alignment)) {
        gene <- cbind(gene, 'Total reads'=rowSums(counts))
        gene <- cbind(gene, data$Alignment)
    }
    
    n_samples <- ncol(counts)
    
    if (!quiet)
        cat(sprintf("%d genes\n", nrow(counts)))
    
    have.norm <- !is.null(norm.file)
    if (have.norm) {
        #norm <- read.delim(norm.file, check.names=FALSE)
        #rownames(norm) <- norm$Sample
        norm <- read.grouped.table(norm.file, require=c('All'))$All
        
        # Make counts have same order as norm file, throw out columns not included
        counts <- counts[,rownames(norm),drop=FALSE]
        n_samples <- ncol(counts)
    }
    
    if (!is.null(keep)) {
        counts <- counts[,keep,drop=FALSE]
        if (have.norm)
            norm <- norm[keep,,drop=FALSE]
        n_samples <- ncol(counts)
    }
    
    totals <- mapply(function(i) sum(counts[i,]), 1:nrow(counts))
    maximums <- mapply(function(i) max(counts[i,]), 1:nrow(counts))    
    good <- totals >= min.total & maximums >= min.max

    if (!quiet)
        cat(sprintf("%d genes after filtering\n", sum(good)))
    
    result <- DGEList(counts=counts[good,], gene=gene[good,])

    mean.lib.size <- exp(mean(log(result$samples$lib.size)))

    if (!have.norm) {
        if (use.tmm) {
            result <- calcNormFactors(result)
        } else {
            result$samples$norm.factors <- rep(1, n_samples)
        }
        
        effective.sizes <- result$samples$lib.size * result$samples$norm.factors
        result$samples$normalizing.multiplier <- mean.lib.size / effective.sizes 
    } else {
        result$samples$normalizing.multiplier <- norm$Normalizing.multiplier
        effective.sizes <- mean.lib.size / result$samples$normalizing.multiplier
        result$samples$norm.factors <- effective.sizes / result$samples$lib.size
    }
    
    result$original.number.of.genes <- nrow(counts)    

    if (!quiet) {
        cat('                        Sample    Library size   Further         Normalizing\n')
        cat('                                                 norm factor     multiplier\n')
        for(i in 1:n_samples) {
            cat(sprintf("%30s %15d %10.2f %15.2f\n", colnames(result$counts)[i], result$samples$lib.size[i], result$samples$norm.factors[i], result$samples$normalizing.multiplier[i]))
        }
        cat('\n\n')
    }
    
    result
}

# Normalize to reads per million (+/- a few TMM tweaks)
reads.per.million <- function(dgelist) {
    effective.sizes <- dgelist$samples$lib.size * dgelist$samples$norm.factors
    t( t(dgelist$counts)/effective.sizes * 1e+06 )
}

# Normalize roughly preserving the overall count magnitude
normalized.counts <- function(dgelist) {
    t( t(dgelist$counts) * dgelist$samples$normalizing.multiplier )
}

# Similar to a log transformation, however glog2(0, m) = log2(m) rather than negative infinity.
#
# Stabilizes variance of values with a constant noise component plus a componet that scales with the data
#
glog2 <- function(x,m=1) log2((x+sqrt(x*x+4*m*m))/2)

# glog2 transformation, resulting in moderated log2 Reads Per Million values
# 
# B.P. Durbin, J.S. Hardin, D.M. Hawkins and D.M. Rocke (2002)
# A variance-stabilizing transformation for gene-expression microarray data.
# Bioinformatics (2002) 18 (suppl 1): S105-S110. 
#
glog2.rpm.counts <- function(dgelist, moderation=5.0) {
    lib.size <- dgelist$samples$lib.size * dgelist$samples$norm.factors
    mean.size <- mean(lib.size)

    #e <- glog2( t(t(dgelist$counts)*(1e6/lib.size)), moderation*1e6/mean.size )
    #e <- log2( t(t(dgelist$counts)*(1e6/lib.size)) + moderation*1e6/mean.size ) 
    e <- matrix(nrow=nrow(dgelist$counts),ncol=ncol(dgelist$counts))
    rownames(e) <- rownames(dgelist$counts)
    colnames(e) <- colnames(dgelist$counts)
    for(i in basic.seq(ncol(e))) {
        e[,i] <- glog2( dgelist$counts[,i] * (1e6/lib.size[i]), moderation*(1e6/mean.size))
        #e[,i] <- glog2(dgelist$counts[,i]*1.0, moderation) + log2(1e6/lib.size[i])
    }
    
    out <- list()
    out$E <- e
    out$lib.size <- lib.size
    out$genes <- dgelist$genes
    out$targets <- dgelist$samples
    new("EList", out)
}





row.apply <- function(data, func) {
    mapply(function(i) func(data[i,]), 1:nrow(data))
}



