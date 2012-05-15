
# Compare outputs from "test-counts:" (can be read using read.delim)
# or any data frame with columns "Feature" and "FDR", sorted by significance

compare.tests <- function(filename.prefix, a, b, n=50, title='', a.title='', b.title='') {
    png(sprintf('%s.png', filename.prefix), width=800,height=20*n)
    plot.new()
    plot.window(c(0,30), c(0,n+4), mar=c(0,0,0,0))
    
    text(15,n+3, title,adj=0.5)
    text(10,n+1.5, a.title,adj=1)
    text(20,n+1.5, b.title,adj=0)
    text(1,n+1.5,'FDR',adj=0)
    text(25,n+1.5,'FDR',adj=0)
    
    for(i in 1:n) {
        text(10,n+1-i, a$Feature[i], adj=1)
        text(1,n+1-i, sprintf('%.2g',a$FDR[i]), adj=0)
        
        text(20,n+1-i, b$Feature[i], adj=0)
        text(25,n+1-i, sprintf('%.2g',b$FDR[i]), adj=0)
    }
    
    for(i in 1:n)
      for(j in 1:n)
        if (as.character(a$Feature[i]) == as.character(b$Feature[j])) 
          lines(c(10,11,19,20),c(n+1-i,n+1-i,n+1-j,n+1-j))
    
    dev.off()
}
