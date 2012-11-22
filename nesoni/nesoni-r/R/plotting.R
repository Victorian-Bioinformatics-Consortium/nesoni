
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

color.legend <- function(col, breaks, title='') {
    if (min(breaks) != max(breaks)) {
        basic.image(x=breaks,y=c(0),z=matrix(breaks,ncol=1), col=col,breaks=breaks)
        axis(1, las=2, line=0.5)
        mtext(title, line=0.5)
    }
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

do.dendrogram <- function(mat, enable=TRUE) {
    # Note: paths are given ordered by order

    if (nrow(mat) < 3 || !enable) {
        list(
            dendrogram = NULL,
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

trim.labels <- function(labels) {
    labels <- as.character(labels)
    n <- max(10, ceiling( mean(nchar(labels)) * 3.0 ))
    
    for(i in basic.seq(length(labels)))
        if (nchar(labels[i]) > n) {
            m <- n
            for(j in 1:n)
                if (any(substr(labels[i],j,j) == c(' ',',',';')))
                   m <- j
            labels[i] <- sprintf("%s...",substr(labels[i],1,m))
        }
    
    labels
}


multiplot <- function(plots, labels) {
    n.plot <- length(plots)
    
    height.inches <- par()$fin[2]

    dend.col.y1 <- 1.0 - 1.5/height.inches
    dend.col.y2 <- 1.0 - 0.1/height.inches
        
    y1 <- min(4 / height.inches, 0.4)
    y2 <- max(dend.col.y1 - 0.1/height.inches, 0.5)

    legend.y2 <- y1*0.75
    legend.y1 <- legend.y2 - 0.03*aspect.ratio()
    
    weights <- as.vector( mapply(function(i)i$weight, plots) )
    cumweights <- c(0,cumsum(weights))
    
    print(cumweights)

    annotation.x <- 0.5
    x1 <- numeric(n.plot)
    x2 <- numeric(n.plot)
    for(i in basic.seq(n.plot)) {
        margin <- 0.05 * weights[i] / cumweights[n.plot+1]
        x1[i] <- annotation.x * (cumweights[i] / cumweights[n.plot+1] + margin)
        x2[i] <- annotation.x * (cumweights[i+1] / cumweights[n.plot+1] - margin)
    }
    
    plot.new()

# = Do plots =========================================================================

    legend.x1 <- annotation.x + 0.05

    for(i in basic.seq(n.plot)) {
        plot <- plots[[i]]
        if (plot$type == 'dendrogram') {
            dendrogram <- plot$data$dendrogram
            if (!is.null(dendrogram)) {
                put.plot(x1[i],x2[i],y1,y2)
                plot(dendrogram, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
            }            
        
        } else if (plot$type == 'heatmap') {
            dendrogram <- plot$dendrogram$dendrogram
            if (!is.null(dendrogram)) {
                put.plot(x1[i],x2[i],dend.col.y1,dend.col.y2)
                plot(dendrogram, horiz=FALSE, axes=FALSE, yaxs="i", leaflab="none")                
            }

            if (plot$signed) {
                col <- signed.col
                extreme <- max(0.0,abs(plot$data),na.rm=TRUE)
                breaks <- seq(-extreme,extreme, length=length(col)+1)
            } else {
                col <- unsigned.col
                extreme <- max(0,plot$data,na.rm=TRUE)
                breaks <- seq(0,extreme, length=length(col)+1)
            }
            
            if (!is.null(plot$legend)) {
                legend.x2 <- legend.x1 + 0.1
                put.plot(legend.x1,legend.x2,legend.y1,legend.y2)
                color.legend(col, breaks, plot$legend)
                legend.x1 <- legend.x2 + 0.05       
            }
        
            put.plot(x1[i],x2[i],y1,y2)
            data <- plot$data
            if (nrow(data) == 1) data <- data[c(1,1),]
            if (ncol(data) == 1) data <- data[,c(1,1)]
            basic.image(basic.seq(ncol(data))-0.5, basic.seq(nrow(data))-0.5, t(data), col=col, breaks=breaks)            
            axis(1, at=basic.seq(ncol(data))-0.5, labels=colnames(data), las=2, tick=FALSE)        
            mtext(plot$title, adj=0.5, line=0.5)
                        
        } else if (plot$type == 'bar') {
            if (length(plot$data)) {
                put.plot(x1[i],x2[i],y1,y2)
                
                barplot(as.vector(plot$data), space=0, horiz=TRUE, col='black', yaxs='i',xaxt='n',lab=c(3,3,7),labels=NULL)
                axis(1,line=0.5,las=2)
                mtext(plot$title, adj=0.5, line=0.5)
            }
        } else if (plot$type == 'scatter') {
            if (length(plot$data)) {
                put.plot(x1[i],x2[i],y1,y2)
                
                plot(as.vector(plot$data),basic.seq(length(plot$data))-0.5, ylim=c(0,length(plot$data)),yaxt='n',ylab=NA,yaxs='i',xaxt='n',xlab=NA,pch=18,bty='n',lab=c(3,3,7))
                axis(1,line=0.5,las=2)
                mtext(plot$title, adj=0.5, line=0.5)
            }
        }
    }    

# = Do annotations ===========================================================================
#
#   Depends on there being at least one plot, as this will be an axis on the
#   last plot drawn.

    line <- 0
    for(i in basic.seq(length(labels))) {
        if (i < length(labels))
            l <- trim.labels(labels[[i]])
        else
            l <- as.character(labels[[i]])
        axis(4, at=basic.seq(length(l))-0.5, labels=l, las=2, tick=FALSE, line=line)
        line <- line + 0.5 * (1+max(0,nchar(l)))
    }
}



nesoni.heatmap <- function(mat, 
                           labels=list(), 
                           reorder.columns=FALSE, 
                           sort.mat=NULL, 
                           signed=TRUE,
                           legend='log2 Reads Per Million\ndifference from row mean',
                           levels=NULL) {    
    n.rows <- nrow(mat)
    n.cols <- ncol(mat)
    
    if (is.null(sort.mat))
       sort.mat <- t(scale(t(mat)))
    
    dend.row <- do.dendrogram( sort.mat )
    dend.col <- do.dendrogram( t(sort.mat), enable=reorder.columns )
    
    plots <- list(
            list(
                weight=1,
                type='dendrogram',
                data=dend.row
            ),
            list(
                weight=2,
                type='heatmap',
                data=mat[dend.row$order,dend.col$order,drop=FALSE],
                signed=signed,
                dendrogram=dend.col,
                legend=legend
            )
    )
    
    if (!is.null(levels)) {
        plots[[3]] <- list(
            weight=0.25,
            type='bar',
            data=levels[dend.row$order],
            title='row\nmean'
        )
    }
    
    multiplot(
        plots,
        lapply(labels, function(item) item[dend.row$order])
    )

    list(
        dend.row = dend.row,
        dend.col = dend.col
    )   
    
    #if (n.rows < 2 || n.cols < 2) {
    #    plot.new()
    #    title(sprintf('Can\'t plot %d x %d heatmap', n.rows, n.cols))
    #} else {
    #    dend.row.x1 <- 1/90
    #    dend.row.x2 <- 1/9
    #    
    #    dend.col.y1 <- 1.0 - dend.row.x2 * aspect.ratio()
    #    dend.col.y2 <- 1.0 - (1.0-dend.col.y1)*0.1
    #    
    #    y1 <- 4 / par()$fin[2]
    #    y2 <- dend.col.y1
    #    
    #    x1 <- dend.row.x2
    #    x2 <- 1/3
    #
    #    legend.y2 <- y1*0.75
    #    legend.y1 <- legend.y2 - 0.03*aspect.ratio()
    #    legend.x1 <- 15/30
    #    legend.x2 <- 18/30
    #
    #    if (signed) {
    #        col <- signed.col
    #        extreme <- max(0.0,abs(mat),na.rm=TRUE)
    #        breaks <- seq(-extreme,extreme, length=length(col)+1)
    #    } else {
    #        col <- unsigned.col
    #        extreme <- max(0,mat,na.rm=TRUE)
    #        breaks <- seq(0,extreme, length=length(col)+1)
    #    }
    #            
    #    plot.new()
    #    
    #    if (!is.null(dend.col$dendrogram)) {
    #        put.plot(x1,x2, dend.col.y1,dend.col.y2)
    #        plot(dend.col$dendrogram, horiz=FALSE, axes=FALSE, yaxs="i", leaflab="none")
    #    }
    #    
    #    if (!is.null(dend.row$dendrogram)) {
    #        put.plot(dend.row.x1,dend.row.x2, y1,y2)
    #        plot(dend.row$dendrogram, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
    #    }
    #    
    #    put.plot(x1,x2, y1,y2)
    #    basic.image(1:n.cols,1:n.rows, t(mat[dend.row$order,dend.col$order,drop=FALSE]), col=col, breaks=breaks)
    #    axis(1, at=1:n.cols, labels=colnames(mat), las=2, tick=FALSE)
    #
    #    line <- 0
    #    for(i in basic.seq(length(labels))) {
    #        if (i < length(labels))
    #            l <- trim.labels(labels[[i]])
    #        else
    #            l <- as.character(labels[[i]])
    #        axis(4, at=1:n.rows, labels=l[dend.row$order], las=2, tick=FALSE, line=line)
    #        line <- line + 0.45 * max(0,nchar(l))
    #    }
    #    
    #    put.plot(legend.x1,legend.x2, legend.y1,legend.y2)
    #    color.legend(col, breaks, legend)        
    #}
    #
    #list(
    #    dend.row = dend.row,
    #    dend.col = dend.col
    #)   
}

#nesoni.heatmap <- function(x, reorder.columns=FALSE,dist.row=NA,dist.col=NA, ...) {    
    #if(!is.matrix(x)) x <- as.matrix(x)
    #
    #if (all(is.na(dist.row))) 
    #    dist.row <- dist( t(scale(t(x))) )
    #
    #if (all(is.na(dist.col)) && reorder.columns) 
    #    dist.col <- dist( scale(t(x)) )
    #
    #method <- 'OLO'
    #
    #if (reorder.columns) {   
    #    control <- list(hclust = hclust(dist.col))
    #    dend.col <- as.dendrogram(
    #            seriation::seriate(dist.col, 
    #            method = method, control = control)[[1]])
    #    dend <- 'both'
    #} else {    
    #    dend.col <- NA
    #    dend <- 'row'
    #}
    #
    #control <- list(hclust = hclust(dist.row))
    #dend.row <- as.dendrogram(
    #        seriation::seriate(dist.row, 
    #        method = method, control = control)[[1]])
    #
    #result <- gplots::heatmap.2(x, Colv = dend.col, Rowv = dend.row, dendrogram = dend, 
    #                  scale = "none", trace='none', density.info='none', 
    #                  lwid=c(2,par("din")[1]-2), lhei=c(1.5,par("din")[2]-1.5),
    #                  cexCol=1.0,
    #                  ...)
    #
    #invisible(result)
#}    

unsigned.col <- hsv(h=seq(0.95,1.15, length.out=256)%%1.0, v=seq(0,1, length.out=256)**0.5,s=seq(1,0,length.out=256)**0.5)
signed.col <- hsv(h=(sign(seq(-1.0,1.0, length.out=256))*0.2+0.8)%%1.0, v=1,s=abs(seq(-1,1,length.out=256)))


svd.gene.picker <- function(mat, svd.rank=NULL, min.svd=2.0) {
    if (is.null(svd.rank))
        svd.rank <- ncol(mat)-1
    s <- svd(t(scale(t(mat), center=TRUE,scale=FALSE)), nu=svd.rank,nv=svd.rank)
    cat('SVD d diagonal:\n')
    print(s$d[basic.seq(svd.rank)])
    mag <- sqrt( rowSums(s$u*s$u) * nrow(s$u) / ncol(s$u) )
    mag >= min.svd
}


hmap.elist <- function(filename.prefix, elist, 
                       min.sd=0.0, min.span=0.0, min.svd=0.0, svd.rank=NULL,
                       annotation=c('gene', 'product'), 
                       res=150, row.labels=NA,
                       reorder.columns = FALSE) {
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
        #if (is.null(svd.rank))
        #    svd.rank <- ncol(elist$E)-1
        #s <- svd(t(scale(t(elist$E), center=TRUE,scale=FALSE)), nu=svd.rank,nv=svd.rank)
        #cat('SVD d diagonal:\n')
        #print(s$d[basic.seq(svd.rank)])
        #mag <- sqrt( rowSums(s$u*s$u) * nrow(s$u) / ncol(s$u) )
        #keep <- (keep & mag >= min.svd)
        
        keep <- keep & svd.gene.picker(elist$E, svd.rank, min.svd)
    }
    
    elist <- elist[keep,]
    
    averages <- rowMeans(elist$E)
    
    if (is.na(row.labels))
        row.labels <- (nrow(elist) <= 300)

    data <- t(scale(t(elist$E), center=TRUE,scale=FALSE))

    labels <- list(rownames(data))
    
    for(colname in annotation)
        if (!all(is.na(elist$gene[,colname])))
            labels[[ length(labels)+1 ]] <- elist$gene[,colname]
        
    height <- if(row.labels) (25*nrow(data)+800)*res/150 else 2500*res/150    
    png(sprintf('%s.png',filename.prefix), width=2000*res/150, height=height, res=res)
    
    #heatmap <- nesoni.heatmap(data, col=signed.col, symkey=TRUE,symbreaks=TRUE, labRow=(if(row.labels) NULL else NA), margins=margins, main=main, ...)
    heatmap <- nesoni.heatmap(data, labels=if(row.labels) labels else list(), reorder.columns=reorder.columns, levels=averages)

    dev.off()

    shuffled.elist <- elist[rev(heatmap$dend.row$order),]
    
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
    #frame[,'cluster hierarchy'] <- rev(dendrogram.paths(heatmap$rowDendrogram))
    frame[,'cluster hierarchy'] <- rev(heatmap$dend.row$paths)
    frame <- data.frame(frame, shuffled.elist$E, check.names=FALSE)
    
    write.csv(frame, row.names=FALSE)
    sink()

    invisible(list(heatmap=heatmap, frame=frame))
}

