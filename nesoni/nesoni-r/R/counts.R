
# Functions for dealing with table created by "nesoni count:"

read.counts <- function(filename, min.total=0, min.max=0, keep=NULL, norm.file=NULL, use.tmm=TRUE) {    
    #data <- read.delim(filename, check.names=FALSE)
    #rownames(data) <- data$Feature
    #n_samples <- grep('^RPKM', colnames(data))[1] - 2
    #counts <- as.matrix( data[,2:(n_samples+1), drop=FALSE] )
    #gene <-data[, (n_samples*2+2):ncol(data)]
    
    data <- read.grouped.table(filename, require=c('Count','Annotation'))
    counts <- as.matrix( data$Count )
    gene <- data$Annotation
    n_samples <- ncol(counts)
    
    
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

    cat('                        Sample    Library size   Further         Normalizing\n')
    cat('                                                 norm factor     multiplier\n')
    for(i in 1:n_samples) {
        cat(sprintf("%30s %15d %10.2f %15.2f\n", colnames(result$counts)[i], result$samples$lib.size[i], result$samples$norm.factors[i], result$samples$normalizing.multiplier[i]))
    }
    cat('\n\n')
    
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

# log transformation with various moderation options, resulting in log2 RPM values
# moderation can be
# - a number 
# - 'vst' to use DESeq's Variance Stabilizing Transformation
log.rpm.counts <- function(dgelist, moderation=5.0) {
    lib.size <- dgelist$samples$lib.size * dgelist$samples$norm.factors
    mean.size <- mean(lib.size)

    if (moderation == 'vst') {
        cds <- newCountDataSet(dgelist$counts, basic.seq(ncol(dgelist$counts)))
        pData(cds)$sizeFactor <- (dgelist$samples$lib.size * dgelist$samples$norm.factors) / 1e6
        cdsBlind <- estimateDispersions(cds, method='blind', fitType='local')
        e <- getVarianceStabilizedData(cdsBlind)
        
    } else {
        moderation <- as.numeric(moderation)
        is.na(moderation) && stop("expected log moderation to be either a number or 'vst'")        
        e <- t( log2((t(dgelist$counts)/lib.size + moderation/mean.size) * 1e+06) )
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



