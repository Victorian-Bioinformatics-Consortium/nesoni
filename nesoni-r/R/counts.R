
# Functions for dealing with table created by "nesoni count:"

read.counts <- function(filename, min.total=0, min.max=0, keep=NULL, norm.file=NULL, use.tmm=TRUE) {    
    data <- read.delim(filename, check.names=FALSE)
    rownames(data) <- data$Feature
    n_samples <- grep('^RPKM', colnames(data))[1] - 2
    counts <- as.matrix( data[,2:(n_samples+1), drop=FALSE] )
    gene <-data[, (n_samples*2+2):ncol(data)]
    
    cat(sprintf("%d genes\n", nrow(counts)))
    
    have.norm <- !is.null(norm.file)
    if (have.norm) {
        norm <- read.delim(norm.file, check.names=FALSE)
        rownames(norm) <- norm$Sample
        
        # Make counts have same order as norm file, throw out columns not included
        counts <- counts[,rownames(norm),drop=FALSE]
        n_samples <- ncol(counts)
    }
    
    if (!is.null(keep)) {
        counts <- counts[,keep,drop=FALSE]
        if (have.norm)
            norm <- norm[keep,]
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



row.apply <- function(data, func) {
    mapply(function(i) func(data[i,]), 1:nrow(data))
}

dendrogram.paths <- function(dend) {
   if (is.leaf(dend)) {
       ''
   } else {
       result <- c()
       for(i in 1:length(dend)) {
           result <- c(result, paste(sprintf('%d',length(dend)-i), dendrogram.paths(dend[[i]])))
       }
       result
   }
}



nesoni.heatmap <- function(x, reorder.columns=FALSE,dist.row=NA,dist.col=NA, ...) {    
    if(!is.matrix(x)) x <- as.matrix(x)

    if (all(is.na(dist.row))) 
        dist.row <- dist( t(scale(t(x))) )
    
    if (all(is.na(dist.col)) && reorder.columns) 
        dist.col <- dist( scale(t(x)) )
    
    method <- 'OLO'

    if (reorder.columns) {   
        control <- list(hclust = hclust(dist.col))
        dend.col <- as.dendrogram(
                seriation::seriate(dist.col, 
                method = method, control = control)[[1]])
        dend <- 'both'
    } else {    
        dend.col <- NA
        dend <- 'row'
    }
    
    control <- list(hclust = hclust(dist.row))
    dend.row <- as.dendrogram(
            seriation::seriate(dist.row, 
            method = method, control = control)[[1]])

    result <- gplots::heatmap.2(x, Colv = dend.col, Rowv = dend.row, dendrogram = dend, 
                      scale = "none", trace='none', density.info='none', 
                      lwid=c(2,par("din")[1]-2), lhei=c(1.5,par("din")[2]-1.5),
                      cexCol=1.0,
                      ...)
    
    invisible(result)
}    

unsigned.col <- hsv(h=seq(0.95,1.15, length.out=256)%%1.0, v=seq(0,1, length.out=256)**0.5,s=seq(1,0,length.out=256)**0.5)
signed.col <- hsv(h=(sign(seq(-1.0,1.0, length.out=256))*0.2+0.8)%%1.0, v=1,s=abs(seq(-1,1,length.out=256)))

hmap.elist <- function(filename.prefix, elist, 
                       min.sd=0.0, min.span=0.0, min.svd=0.0, svd.rank=NULL,
                       annotation=c('gene', 'product'), 
                       res=150, row.labels=NA, margins=c(20,20), 
                       main='log2 expression\ndifference from\nrow average', ...) {
    keep <- rep(TRUE, nrow(elist$E))

    if (min.sd > 0.0) {
        sd <- sqrt(row.apply(elist$E, var))
        keep <- (keep & sd >= min.sd)
    }
    
    if (min.span > 0.0) {    
        span <- row.apply(elist$E, max) - row.apply(elist$E, min)
        keep <- (keep & span >= min.span)
    }
    
    if (min.svd > 0.0) {
        if (is.null(svd.rank))
            svd.rank <- ncol(elist$E)-1
        s <- svd(t(scale(t(elist$E), center=TRUE,scale=FALSE)), nu=svd.rank,nv=svd.rank)
        cat('SVD d diagonal:\n')
        print(s$d[1:svd.rank])
        mag <- sqrt( rowSums(s$u*s$u) * nrow(s$u) / ncol(s$u) )
        keep <- (keep & mag >= min.svd)
    }
    
    elist <- elist[keep,]
    
    if (is.na(row.labels))
        row.labels <- (nrow(elist) <= 500)

    data <- t(scale(t(elist$E), center=TRUE,scale=FALSE))
    
    for(colname in annotation)
        if (!all(is.na(elist$gene[,colname])))
            rownames(data) <- paste(rownames(data), elist$gene[,colname])
    
    height <- if(row.labels) (16*nrow(data)+400)*res/150 else 2500*res/150    
    png(sprintf('%s.png',filename.prefix), width=1500*res/150, height=height, res=res)
    
    heatmap <- nesoni.heatmap(data, col=signed.col, symkey=TRUE, labRow=(if(row.labels) NULL else NA), margins=margins, main=main, ...)

    dev.off()

    shuffled.elist <- elist[rev(heatmap$rowInd),]
    
    table.filename <- sprintf('%s.csv', filename.prefix)
    
    sink(table.filename)
    cat('# Heatmap data\n')
    cat('#\n')
    cat('# Values given are log2 reads per million\n')
    cat('#\n')
    cat(sprintf('# %d genes shown\n', nrow(data)))
    cat('#\n')

    frame <- data.frame(name=rownames(shuffled.elist$E), row.names=rownames(shuffled.elist$E), check.names=FALSE) 
    for(colname in annotation) {
        frame[,colname] <- shuffled.elist$gene[,colname]
    }     
    frame[,'cluster hierarchy'] <- rev(dendrogram.paths(heatmap$rowDendrogram))
    frame <- data.frame(frame, shuffled.elist$E, check.names=FALSE)
    
    write.csv(frame, row.names=FALSE)
    sink()

    invisible(list(heatmap=heatmap, frame=frame))
}




