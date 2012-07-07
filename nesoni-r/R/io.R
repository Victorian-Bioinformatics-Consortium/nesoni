
# Read a table
# If it contains a line starting with #Groups, group columns as indicated
# Returns a list group.name -> data.frame

read.grouped.table <- function(filename, require=c()) {
    groups <- c()
    tab.separated <- FALSE
    
    f <- file(filename,'r')
    repeat {
        line <- readLines(f,1)
        if (length(line) == 0) break;
        if (substr(line,1,1) != '#') {
            tab.separated <- (length(grep('\t',line)) > 0)
            break;
        }
        
        parts <- strsplit(line,',')[[1]]
        if (parts[1] == '#Groups') {
            groups <- parts[basic.seq(length(parts)-1)+1]
        }
    }
    close(f)

    if (tab.separated)
        data <- read.delim(filename, comment.char='#', check.names=FALSE)    
    else    
        data <- read.csv(filename, comment.char='#', check.names=FALSE)
    
    rownames(data) <- data[,1]
    data <- data[ ,2:ncol(data), drop=FALSE]


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
        groups <- c('All')
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
