
aspect.ratio <- function() { par()$fin[1] / par()$fin[2] }

put.plot <- function(px1,px2,py1,py2) {
    width <- par()$fin[1]
    height <- par()$fin[2]
    
    par(
        mai=c(0,0,0,0),
        plt=c(px1,px2,py1,py2),
        new=TRUE
    )
}

basic.image <- function(x,y,z, col, breaks) {
    image(x=x, y=y, z=z,col=col,breaks=breaks,xaxt="n",yaxt="n",xlab='',ylab='')
}

color.legend <- function(col, breaks, title='') {
    basic.image(x=breaks,y=c(0),z=matrix(breaks,ncol=1), col=col,breaks=breaks)
    axis(1, las=2)
    mtext(title)
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

basic.seq <- function(n) {
    # Count from 1 to n inclusive without being clever
    
    if (n < 1)
        numeric(0)
    else
        1:n
}

do.dendrogram <- function(mat) {
    if (nrow(mat) < 3) {
        list(
            dendrogram = NA,
            order = basic.seq(nrow(mat)),
            paths = rep('',nrow(mat))
        )
    } else {
    
        dist.mat <- dist(mat)
        control <- list(hclust = hclust(dist.mat))
        dend.mat <- as.dendrogram(
                seriation::seriate(dist.mat, 
                method = 'OLO', control = control)[[1]])
        
        list(dendrogram = dend.mat, 
             order = order.dendrogram(dend.mat),
             paths = dendrogram.paths(dend.mat)
        )
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
    
    heatmap <- nesoni.heatmap(data, col=signed.col, symkey=TRUE,symbreaks=TRUE, labRow=(if(row.labels) NULL else NA), margins=margins, main=main, ...)

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

