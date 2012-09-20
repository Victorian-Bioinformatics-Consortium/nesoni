
basic.seq <- function(n) {
    # Count from 1 to n inclusive without being clever
    
    if (n < 1)
        numeric(0)
    else
        1:n
}

basic.subset <- function(object, rows) {
    # Get rows by boolean vector
        
    stopifnot(is.logical(rows))
    
    # Don't tile selection
    stopifnot( length(rows) == nrow(object) )
    
    # Don't drop down to vector if single row selected
    object[rows,,drop=FALSE]
}

basic.take <- function(object, rows) {
    #Get rows by index

    stopifnot(is.numeric(rows))    
    stopifnot(all(rows >= 1))
    stopifnot(all(rows <= nrow(object)))
    
    object[rows,,drop=FALSE]
}

basic.image <- function(x,y,z, col, breaks) {
    if (length(x) > 0 && length(y) > 0) {
        image(x=x, y=y, z=z,col=col,breaks=breaks,xaxt="n",yaxt="n",xlab='',ylab='')
    }
}
