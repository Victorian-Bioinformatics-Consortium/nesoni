
# NMF on RNAseq expression data
# The weighting scheme is inspired by voom, but is more primitive: 
# it is weighted per-gene rather than per count


nmf.reordering <- function(fac.basis, fac.coef, other.basis, other.coef) {
    common.genes <- intersect(rownames(fac.basis), rownames(other.basis))
    stopifnot(length(common.genes) > 0)

    common.samples <- intersect(colnames(fac.coef), colnames(other.coef))
    stopifnot(length(common.samples) > 0)
    
    similarity <- function(i,j) {
       sum( fac.basis[common.genes,i] * other.basis[common.genes,j] ) *
       sum( fac.coef[i,common.samples] * other.coef[j,common.samples] )       
    }

    n <- ncol(fac.basis)
    n.other <- ncol(other.basis)

    num <- rep(0, n)
    den <- rep(0, n)
    for(i in 1:n) {
        for(j in 1:n.other) {
            s <- similarity(i,j)
            num[i] <- num[i] + s*j
            den[i] <- den[i] + s
        }
    }
    
    order(num / den)
}


# order hint is a list containing $fac.basis and $fac.coef
# (i.e. output of a previous invocation of nmf.counts)
# factorization will attempt to mimic order of this if given
#
nmf.dgelist <- function(dgelist, n.classes, scaling.iters=5, nmf.runs=20, order.hint=NULL) {
    library(NMF) #NMF needs to run some initialization
    
    counts <- normalized.counts(dgelist)

    row.average <- rowMeans(counts)

    # Initial guess is Poisson errors
    #scaler <- sqrt(rowSums(counts))
    scaler <- rowMeans(counts) + 0.5
    
    # Refine scaler    
    for(i in 1:scaling.iters) {
        massaged <- counts / scaler
        fac <- NMF::nmf(massaged, n.classes, method='lee', nrun=nmf.runs)
        fit <- NMF::fitted(fac) * scaler
        
        error <- sqrt(rowSums( (counts - fit)^2 ))
        loess.curve <- loess(log(error) ~ log(row.average), degree=2)
        new.scaler <- exp(loess.curve$fitted)
        
        plot(row.average, error/row.average, log='x')
        points(row.average, scaler/row.average, col='blue')
        points(row.average, new.scaler/row.average, col='red')

        scaler <- new.scaler        
    }

    # Calculate final factorization
    massaged <- counts / scaler
    fac <- NMF::nmf(massaged, n.classes, method='lee', nrun=nmf.runs)
    fac.basis <- basis(fac)
    fac.coef <- coef(fac)

    # Reorder classes
    if (!is.null(order.hint)) {
        reordering <- nmf.reordering(fac.basis, fac.coef, order.hint$fac.basis, order.hint$fac.coef)
    } else {
        reordering <- 1:n.classes
    }
    fac.basis <- fac.basis[, reordering]
    fac.coef <- fac.coef[reordering, ]

    # Remove scaling and
    # normalize factorization such that
    # gene.weight columns sum to 1
    gene.weights <- fac.basis * scaler
    temp <- colSums(gene.weights)
    sample.weights <- t( fac.coef * temp )
    gene.weights <- t(t(gene.weights) / temp)
    
    colnames(sample.weights) <- paste(1:n.classes,'of',n.classes,sep='')
    colnames(gene.weights) <- colnames(sample.weights)

    list(
        dgelist = dgelist,
        counts = counts,
        scaler = scaler,
        fac.object = fac,
        fac.reordering = reordering,
        fac.basis = fac.basis,
        fac.coef = fac.coef,
        gene.weights = gene.weights,
        sample.weights = sample.weights
    )
}


clean.nas <- function(vec, alt) { vec <- as.character(vec); ifelse(is.na(vec),alt,vec) } 

nmf.report <- function(prefix, item, glyph=NULL, annotations=c('gene','product')) {
    saveRDS(item, sprintf('%s.rds', prefix))
    
    n <- ncol(item$gene.weights)
    n.samples <- ncol(item$counts)
    n.genes <- nrow(item$counts)
      
    png(sprintf('%s-heatmap.png',prefix), width=600,height=800, pointsize=18)
    heatmap(scale(item$sample.weights[n.samples:1,],center=FALSE), Colv=NA,Rowv=NA,scale='none',col=unsigned.col,margins=c(3,12))
    dev.off()
    
    p.values <- matrix(1, nrow=n.genes, ncol=n)
    for(i in 1:n.genes) {
        non.zero <- (item$gene.weights[i,] > 1e-10)  #NMF doesn't produce exact zeros
        fit <- lm(item$counts[i,] ~ 0+I(item$sample.weights[,non.zero]))
        
        #Should be equal:
        #print(as.vector( item$gene.weights[i,non.zero] ))
        #print(as.vector( fit$coefficients ))
        
        p.values[i,non.zero] <- as.vector( summary(fit)$coefficients[,'Pr(>|t|)'] ) * 0.5
        #Single sided p values
    }
      
    
    for(j in 1:n) {
        p <- p.values[,j]
        fdr <- p.adjust(p)
    
        result <- data.frame(
            Feature = rownames(item$counts),
            check.names = FALSE
        )
        for(annotation in annotations)
            result[,annotation] <- item$dgelist$genes[,annotation]
            
        result$FDR <- fdr  
        result$Weight <- item$gene.weights[,j]
        result <- data.frame(result, item$gene.weights, check.names=FALSE)

        result <- result[ order(result$FDR), ]
        
        significant <- result$FDR < 0.05
    
        png(sprintf("%s-%s.png",prefix,colnames(item$gene.weights)[j]), width=800, height=800, pointsize=18)
        show <- min(50,nrow(result)):1
        temp <- item$gene.weights[as.character(result$Feature[show]),]
        for(annotation in annotations)
            rownames(temp) <- paste(rownames(temp), clean.nas(result[show, annotation],''))
        rownames(temp) <- paste(ifelse(significant[show], '', '[insignificant] '), rownames(temp), sep='')
        rownames(temp) <- paste(sprintf('%.6f', result$Weight[show]), rownames(temp))
        heatmap(temp / rowSums(temp), Colv=NA,Rowv=NA,scale='none',col=unsigned.col,margins=c(6,36))
        dev.off()
    
        result <- result[significant, ]    
        #write.csv(result, file=sprintf("%s-%s.csv", prefix, colnames(item$gene.weights)[j]), na='', row.names=FALSE)    
        write.table(result, file=sprintf("%s-%s.txt", prefix, colnames(item$gene.weights)[j]), na='', sep='\t', quote=FALSE, row.names=FALSE)    
    }
    
    
    
    base.prefix <- strsplit(prefix,'/')[[1]]
    base.prefix <- base.prefix[ length(base.prefix) ]
    
    sink( sprintf('%s.html', prefix) )
    cat('<html><head><title>', base.prefix, '</title></head><body>\n')
    cat('<h1>', base.prefix, '</h1>\n')
    cat(sprintf('<p><img src="%s-heatmap.png">', base.prefix))
    
    total.sample <- sum(item$sample.weights)
    
    cat('<table cellpadding="10">\n')
    for(j in 1:n) {
        class.name <- colnames(item$gene.weights)[j]
        cat('<tr>\n')
        cat(sprintf('<td valign="top"><b>%s</b><p>%.1f%% of total<p><a href="%s-%s.txt">[table]</a></td>\n', 
            class.name, sum(item$sample.weights[,j]) * 100.0 / total.sample, base.prefix, class.name))

        if (!is.null(glyph)) {
            glyph.name <- sprintf('%s-%s-glyph.png', prefix, class.name)
            glyph(glyph.name, item$sample.weights[,j])
            cat(sprintf('<td valign="top"><img src="%s-%s-glyph.png"/></td>', base.prefix, class.name))
        }

        cat(sprintf('<td valign="top"><div style="height:150px; overflow: hidden; border:1px solid black;"><a href="%s-%s.png"><img style="border: none" src="%s-%s.png"/></a></div></td>\n',
            base.prefix, class.name, base.prefix, class.name)) 
        cat('</tr>\n')
    }
    cat('</table>')
    
    cat('</body></html>\n')
    sink()
}





