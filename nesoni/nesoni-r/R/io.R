
# Read a table
# If it contains a line starting with #Groups, group columns as indicated
# Returns a list group.name -> data.frame

read.grouped.table <- function(filename, require=c(), default.group='All') {
    groups <- c()
    tab.separated <- FALSE
    
    skip <- 0
    
    f <- file(filename,'r')
    repeat {
        line <- readLines(f,1)
        if (length(line) == 0) break;
        if (substr(line,1,1) != '#') {
            tab.separated <- (length(grep('\t',line)) > 0)
            break;
        }
        
        skip <- skip + 1
        
        parts <- strsplit(line,',')[[1]]
        if (parts[1] == '#Groups') {
            groups <- parts[basic.seq(length(parts)-1)+1]
        }
    }
    close(f)

    if (tab.separated)
        data <- read.delim(filename, skip=skip, check.names=FALSE)    
    else    
        data <- read.csv(filename, skip=skip, check.names=FALSE)
    
    rows <- data[,1]
    cols <- colnames(data)[basic.seq(ncol(data)-1)+1]    
    data <- data[ ,basic.seq(ncol(data)-1)+1, drop=FALSE]
    rownames(data) <- rows
    colnames(data) <- cols

    # === Fallbacks if groups not given ===
        
    if (!length(groups)) {
        rpkms <- grep('^RPKM', colnames(data))
        if (length(rpkms)) {
            # === Legacy count file ===
            n_samples <- rpkms[1] - 1
            groups <- c(
                rep('Count', n_samples),
                rep('RPKM', n_samples),
                rep('Annotation', ncol(data)-n_samples*2)
            )
        }
    }
    if (!length(groups)) {
        groups <- c(default.group)
    }

    
    i <- 2
    while(i <= ncol(data)) {
        if (is.null(groups[i]) || is.na(groups[i]) || groups[i] == '')
            groups[i] <- groups[i-1]
        i <- i + 1
    }
        
    groups <- factor(groups)
    
    result <- list()
    for(name in levels(groups)) {
        result[[name]] <- data[,groups == name,drop=FALSE]
    }
    
    for(item in require)
        if (is.null( result[[item]] ))
            stop('Table ',filename,' has no group ',item,'. Is it in the right format?')
    
    result
}


# groups should be a list of data frames
write.grouped.table <- function(groups, filename, comments=c(), rowname.title='Names') {
    # === Sanity check ===
        
    n.groups <- length(groups)
    n.groups >= 1 || stop('Need at least one group of columns to write')
    
    row.names <- rownames(groups[[1]])
    
    for(i in basic.seq(n.groups)) {
        is.data.frame(groups[[i]])                || stop('Not a data frame: ',i,' ',names(groups)[[i]])
        nrow(groups[[i]]) == nrow(groups[[1]])    || stop('Wrong number of rows: ',i,' ',names(groups)[[i]])
        all( rownames(groups[[i]]) == row.names ) || stop('Row names don\'t match: ',i,' ',names(groups)[[i]])
    }

    
    # === Construct data frame ===
    
    column.groups <- character()
    column.names <- character()
    frame <- data.frame(row.names = row.names, check.names=FALSE, check.rows=FALSE)
    
    col <- 1
    column.groups[col] <- ''
    column.names[col] <- rowname.title
    frame[,col] <- row.names
    col <- col+1

    for(i in basic.seq(n.groups)) {
        for(j in basic.seq(ncol(groups[[i]]))) {
            column.groups[col] <- names(groups)[[i]]
            column.names[col] <- colnames(groups[[i]])[j]
            frame[,col] <- groups[[i]][,j]
            col <- col+1
        }
    }
    
    colnames(frame) <- column.names

    # === Don't repeat column group names ===
    
    i <- length(column.groups)
    while(i >= 2) {
        if (column.groups[i] == column.groups[i-1])
            column.groups[i] <- ''
        i <- i - 1
    }


    # === Output ===
   
    sink(filename)
    
    for(comment in comments) {
        cat('#',comment,'\n',sep='')
    }
    
    cat('#Groups',paste(column.groups,collapse=','),'\n',sep='')
    
    write.csv(frame, row.names=FALSE)
    sink()
}





