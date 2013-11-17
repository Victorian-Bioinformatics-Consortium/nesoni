
import sys, os, re, subprocess

from nesoni import grace, config, legion, io, workspace, selection

def R_literal(item):
    if item is None:
        return 'NULL'
    elif isinstance(item, str):
        return "'" + item.replace('\\','\\\\').replace("'","\\'") + "'"
    elif isinstance(item, bool):
        return 'TRUE' if item else 'FALSE'
    elif isinstance(item, float):
        return repr(item)
    elif isinstance(item, int):
        return '%d' % item
    elif isinstance(item, list) or isinstance(item, tuple):
        if not len(item) or not (isinstance(item[0], list) or isinstance(item[0], tuple)):
            return 'c(' + ','.join( R_literal(subitem) for subitem in item ) + ')'
        
        #Matrix
        ncol = len(item[0])
        for item2 in item: assert len(item2) == ncol
        
        all = [ ]
        for item2 in item: all.extend(item2)
        return 'matrix(%s, nrow=%d, byrow=TRUE)' % (
            R_literal(all),
            len(item)
        )
         
    else:
        assert False, "Can't encode %s" % repr(item)

def run_script(script, only_tell=False, silent=False, **kwargs):
    script = (
        '{\n'
        "options(warn=1,error=function() {\n"
        " dump.frames(); cat(names(last.dump),sep='\n'); q('no',1);\n"
        "});\n"
        '\n'
        "if (capabilities('cairo')) {\n"
        "  options(bitmapType='cairo');\n"
        "}\n"
        '\n' +
        ''.join([
            key + ' <- ' + R_literal(kwargs[key]) + '\n'
            for key in kwargs
        ]) + 
        '\n' +
        script +
        '\ninvisible();\n'   #Rscript prints the return value of top level expressions. No thanks.
        '}\n'        
    )
       
    if only_tell:
        print script
        return
        
    if silent:
        stdout = subprocess.PIPE
    else:
        stdout = None
 
    process = subprocess.Popen(
        ['Rscript', '-'],
        bufsize=1<<24,
        stdin=subprocess.PIPE,        
        stdout=stdout,
        stderr=stdout,
        close_fds=True,
    )
    process.stdin.write(script)
    process.stdin.close()
    
    assert 0 == process.wait(), 'R failed'

@config.Bool_flag('tell', 'Show R+ code instead of executing it.')
class R_action(config.Action):
    tell = False

    def run(self):
        args = { }
        for parameter in self.parameters:
            if parameter.name == 'tell': continue
            args[parameter.name] = parameter.get(self)
        run_script(self.script, self.tell, **args)




@config.help("""\
Plot a grid of sample-sample scatter plots.
""")
@config.String_flag('spike_in', 'Comma separated list of spike-in control "genes".')
@config.Positional('counts', 'Output from "count:"')
class Plot_counts(config.Action_with_prefix, R_action):
    counts = None
    spike_in = ''

    script = r"""
    library(nesoni)
    dgelist <- read.counts(counts)    
    counts <- dgelist$counts
    n <- ncol(counts)
    
    pngname <- sprintf('%s-count.png', prefix)
    png(pngname, width=n*300, height=n*300 )
    
    if (nchar(spike_in) == 0) {
        spikes = c()
    } else {
        spikes = strsplit(spike_in,',')[[1]]
    }
    
    not_first <- FALSE
    for(i in 1:n) {
        for(j in 1:n) {
            if (i < j) {
                par(fig=c((i-1)/(n-1),i/(n-1),(n-j)/(n-1),(n-j+1)/(n-1)), new=not_first)
                plot(counts[,i], counts[,j], 
                     log='xy', pch=19, pty='s',
                     xlab=colnames(counts)[i], ylab=colnames(counts)[j])
                
                if (length(spikes)) {
                    points(counts[spikes,i], counts[spikes,j], 
                           pch=19, pty='s', col='red')
                }
                
                not_first <- TRUE
            }
        }
    }
    
    dev.off()
    """


TEST_COUNTS_HELP = """\

Usage:

    nesoni test-counts: [options] output_prefix counts.txt \\
        term [term...] \\
        [contrast: N [N...]] \\
        [with: term [term...]] \\
        [use: regex [regex...]]

Find significant differential expression from the output of "nesoni count:".

term is a regular expression on the sample names, or several 
terms joined by ^, indicating an interaction term (ie the XOR of the 
component terms). "with:" terms are included in the linear model, but 
not the significance test. The model will also include a constant term.

Examples: A, A^B

A term can be given a more intuitive name in the output by appending "=name", eg

Examples: A=fold-change, A^B=interaction-between-A-and-B

Trended common dispersion is used as the dispersion estimate by default.

Options:
   
   --mode MODE      - Within group variation estimation method. See below.
                      Default: voom
   
   --quantile-norm yes/no
                    - Use limma's normalizeQuantiles on
                      log((count+1)/gene_length). The normalization is 
                      then applied via the offset parameter in edgeR.                    
                      Use this if you have banana-shaped scatter plots.                      
                      Default: no
   
   --min-count NNN  - Discard features with less than this total count.
                      Default: 10

   --fdr N.NN       - False Discovery Rate cutoff for statistics and plots.
                      Default: 0.01
   --output-all yes/no
                    - List all genes in output, not just significant ones
                      Default: yes
   
   --constant-term yes/no 
                    - Include a constant term.
                      Default: yes
   
   --tell yes/no    - Output R code instead of executing it
                      Default: no

   term [term...]   - Test these terms
   
   contrast: N N N  - Instead of doing an ANOVA on multiple terms,
                      perform a contrast with the given weights.

   with: term [term...]
                    - Also include these terms in model

   use: regex [regex...] 
                    - Use only samples matching these regexes

Examples:

Test for effect of A

    nesoni test-counts: output counts.txt A

Test for effect of A, including B in the model 
so that within group variance isn't inflated by its effect

    nesoni test-counts: output counts.txt A with: B

Test for effect of A or B or both

    nesoni test-counts: output counts.txt A B

Test for an interaction between A and B

    nesoni test-counts: output counts.txt A^B with: A B

"""


MODE_HELP = {
'voom' : """   
Uses BioConductor package limma, with voom, to perform moderated t-tests.

The empirical prior on within group variance is a function of overall expression.

A plot ...-voom.png is produced showing the relation of within-group 
variance to expression level.

"voomed" columns in the output are log2 reads per million,
normalized using EdgeR's Trimmed Mean normalization.


References:

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing 
differential expression in microarray experiments. Statistical Applications 
in Genetics and Molecular Biology 3, No. 1, Article 3.

Law, CW, Chen, Y, Shi, W, and Smyth, GK (2011). Voom! Variance modelling 
powers an empirical Bayes linear modelling pipeline for RNA-Seq data. Submitted.
""",  
'nullvoom' : """    
Uses BioConductor package limma, with voom, to fit the null model 
to the data. A scaled chi-square distribution is then fitted to 
the residual variance, using maximum likelihood. Genes where the 
residual is significantly larger than would be expected from
the fitted chi-square distribution are reported as significant.

Use this when your biologist (bless them) has not done any 
biological replicates.

This mode is unable to distinguish differentially expressed genes from 
genes with high variability between samples.

This mode gives quite conservative p-values. Also note that this is
not a method endorsed by WEHI, it's just something I made up.

A plot ...-voom.png is produced showing the relation of variance 
to expression level.

A plot ...-qq.png is produced showing how well the chi-square 
distribution fits the residual variance.
""",
'glog' : """
This is a simpler alternative to voom. Counts are converted to 
log2 Reads Per Million using a variance stabilised transformation.

Let the generalised logarithm with moderation m be

   glog(x,m) = log2((x+sqrt(x*x+4*m*m))/2)

then the transformed values will be

   glog( count/library_size*1e6, log_moderation/mean_library_size*1e6 )

where log_moderation is a parameter.

The log2 RPM values are then tested using limma.

References:

B.P. Durbin, J.S. Hardin, D.M. Hawkins and D.M. Rocke (2002)
A variance-stabilizing transformation for gene-expression microarray data.
Bioinformatics (2002) 18 (suppl 1): S105-S110. 

Smyth, G. K. (2004). 
Linear models and empirical Bayes methods for assessing 
differential expression in microarray experiments. 
Statistical Applications in Genetics and Molecular Biology 3, No. 1, Article 3.
""",
'nullglog' : """
Uses BioConductor package limma, after tranformation as in the "glog" method, 
to fit the null model to the data. A scaled chi-square distribution is then 
fitted to the residual variance, using maximum likelihood. Genes where the 
residual is significantly larger than would be expected from the fitted 
chi-square distribution are reported as significant.

Use this when your biologist (bless them) has not done any 
biological replicates.

This mode is unable to distinguish differentially expressed genes from 
genes with high variability between samples.

This mode gives quite conservative p-values. Also note that this is
not a method endorsed by WEHI, it's just something I made up.

A plot ...-qq.png is produced showing how well the chi-square 
distribution fits the residual variance.
""",
'poisson' : """     
Uses BioConductor package EdgeR to look for differential expression 
with a Poisson noise model.

This will certainly give misleading results, and should not be used.


References:

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene
expression data. Bioinformatics 26, 139-140
""",
'common' : """      
Uses BioConductor package EdgeR to look for differential expression
using a negative binomial model, with the assumption that overdispersion 
is the same for all tags.

This mode is unable to distinguish differentially expressed gene from 
genes with high variability between samples.


References:

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene
expression data. Bioinformatics 26, 139-140
""",
'trend' : """       
Uses BioConductor package EdgeR to look for differential expression
using a negative binomial model, with the assumption that overdispersion 
is a smoothly varying function of overall intensity.

This mode is unable to distinguish differentially expressed gene from 
genes with high variability between samples.


References:

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene
expression data. Bioinformatics 26, 139-140
"""
}



COMMON = r"""

describe <- function(vec) {
    vec <- vec[!is.nan(vec)]
    
    sprintf(
         "mean %.2f  sd %.2f  min %.2f  median %.2f  max %.2f",
         mean(vec),
         sd(vec),
         min(vec),
         median(vec),
         max(vec)
    );
}


do.heatmap <- function(heatmap.features, heatmap.data, heatmap.labels, heatmap.values) {
    png(sprintf('%s-heatmap.png',OUTPUT_PLOT), res=150, width=2000, height=min(5000,25*length(heatmap.features)+800))

    if (ncol(heatmap.values) == 1) {
        ord <- order(heatmap.values[,1])
    } else {        
        ord <- do.dendrogram( t(scale(t(heatmap.data), scale=FALSE)) )$order
    }
    
    cols <- list(
            list(
                weight=1,
                type='heatmap',
                data=t(scale(t(heatmap.data),center=TRUE,scale=FALSE)) [ord,,drop=FALSE],
                signed=TRUE,
                legend='log2 Reads Per Million\ndifference from row mean'
            ),
            list(
                weight=0.25,
                type='bar',
                data=rowMeans(heatmap.data) [ord],
                title='mean\nlog2 RPM'
            )
        )
        
    for(i in basic.seq(ncol(heatmap.values)))
        cols[[length(cols)+1]] <- list(
                weight=0.25,
                type='bar',
                data=heatmap.values[,i][ord],
                title=colnames(heatmap.values)[i]
            )
    
    heatmap.labels.ord <- lapply(heatmap.labels,function(item) item[ord])
    
    multiplot(
        cols,
        heatmap.labels.ord
    )           
    dev.off()
}           


get.significant <- function(result) {
    significant <- rep(TRUE, nrow(result))
    if (LOG_FOLD_CUTOFF > 0.0) {
        for(i in basic.seq(N_TO_TEST)) {
            significant <- significant & (abs(result[,DESIGN_COLUMNS[i]]) >= LOG_FOLD_CUTOFF)
        }
    }
    
    significant <- significant & (result$FDR <= FDR_CUTOFF)
    significant
}

describe.significant <- function() {
    desc <- sprintf('a False Discovery Rate of %g', FDR_CUTOFF)
    if (LOG_FOLD_CUTOFF > 0.0) {
        desc <- sprintf('%s and an absolute log2 fold change of at least %g', desc, LOG_FOLD_CUTOFF)
    }
    desc
}


library(nesoni)
dgelist <- read.counts(FILENAME, min.total=MIN_COUNT, keep=KEEP, norm.file=NORM_FILE)

design <- DESIGN
colnames(design) <- DESIGN_COLUMNS

"""


EDGER = COMMON + r"""

library(edgeR)
library(limma)

stopifnot(!QUANTILE_NORM)

#Hmm
if (nrow(dgelist$counts) == 0) {
    sink(LOG_FILENAME, split=TRUE, append=TRUE)
    cat(nrow(deglist$counts), 'tags with total count at least', MIN_COUNT, '\n')
    cat('Insufficient data for a meaningful answer.\n')
} else {
    d <- dgelist
    
    if (MODE == 'trend') {
        cat('estimateGLMTrendedDisp\n')
        d <- estimateGLMTrendedDisp(d, design=design)
        #dispersion <- d$trended.dispersion
        cat('estimateGLMTagwiseDisp\n')
        d <- estimateGLMTagwiseDisp(d, design=design)
        dispersion <- d$tagwise.dispersion
        summary(dispersion)
    } else if (MODE == 'common') {
        cat('estimateGLMCommonDisp\n')
        d <- estimateGLMCommonDisp(d, design=design)
        #dispersion <- d$common.dispersion
        print(d$common.dispersion)
        cat('estimateGLMTagwiseDisp\n')
        d <- estimateGLMTagwiseDisp(d, design=design)
        dispersion <- d$tagwise.dispersion
        summary(dispersion)
    } else {
        cat('Poisson dispersion\n')
        dispersion <- 1e-6    #See edgeR Poisson example
    }
    
    cat('glmFit\n')
    fit <- glmFit(d, design=design, dispersion=dispersion)
    
    cat('glmLRT\n')
    #if (USE_CONTRAST) {
    #    contrast <- CONTRAST_WEIGHTS
    #    print(contrast)
    #    lrt <- glmLRT(d, fit, contrast=contrast)
    #} else {
    
    #API changed in edgeR
    # lrt <- glmLRT(d, fit, coef=1:N_TO_TEST)
    lrt <- glmLRT(fit, coef=1:N_TO_TEST)
    
    #}
    
    # Construct result frame
    result <- data.frame(
        Feature = rownames(dgelist$counts),
        row.names = rownames(dgelist$counts)
    )
    
    # 10/9/2012: logConc changed to logCPM (concentration per million?)
    result$"log2 average per million" <- lrt$table$logCPM

    #if (USE_CONTRAST) {
    #    result$"log2 contrast" <- lrt$table$logFC
    #}
    
    for(i in 1:N_TO_TEST) {
        result[, colnames(design)[i]] <- lrt$coefficients[,i] 
        #10/9/2012: Coefficients full now appears to be in log2
        #26/9/2013: Renamed coefficients.full -> coefficients
    }
                    
    #result$p <- lrt$table$p.value
    #result$FDR <- p.adjust(lrt$table$p.value, method='BH')
    #10/9/2012: p.value appears to have been renamed
    result$p <- lrt$table$PValue
    result$FDR <- p.adjust(lrt$table$PValue, method='BH')
    
    #for(i in (N_ALL_SAMPLES*2+2):ncol(data)) {
    #    result[ , colnames(data)[i]] = data[, i]
    #}
    #for(i in 2:(N_ALL_SAMPLES*2+1)) {
    #    result[ , colnames(data)[i]] = data[, i]
    #}
    
    for(i in 1:ncol(dgelist$genes)) {
        result[, colnames(dgelist$genes)[i]] = dgelist$genes[i]
    }
    
    reordering <- order(result$p)
    result <- result[ reordering,, drop=FALSE]
    rcounts <- dgelist$counts[ reordering,, drop=FALSE]
    
    #significant <- (result$FDR <= FDR_CUTOFF)
    significant <- get.significant(result)
        
    results_to_output <- result
    if (OUTPUT_COUNTS) {
        # Add raw counts
        results_to_output <- cbind(results_to_output, dgelist$counts[rownames(results_to_output),])
    }

    blank <- mapply(function(i){ all(is.na(results_to_output[,i])) }, 1:ncol(results_to_output))

    sink(OUTPUT_FILENAME_ALL)
    write.csv(results_to_output[,!blank,drop=FALSE], na='', row.names=FALSE)
    sink()

    sink(OUTPUT_FILENAME_SIGNIFICANT)
    write.csv(results_to_output[significant,!blank,drop=FALSE], na='', row.names=FALSE)
    sink()
    
    #if (USE_CONTRAST) {
    #    pngname = sprintf('%s.png', OUTPUT_PLOT)
    #    png(pngname, width=800, height=800 )
    #
    #    sane <- abs(result[,'log2 contrast']) < 1000.0
    #
    #    plot(result[sane,"log2 average per million"], 
    #         result[sane,"log2 contrast"],
    #         xlab="log2 average per million",
    #         ylab="log2 contrast",
    #         main='Contrast',
    #         pch=19, cex=0.25) 
    #    points(result[significant & sane,"log2 average per million"],
    #           result[significant & sane,"log2 contrast"],
    #           pch=19, cex=1.0,
    #           col="red")
    #    dev.off()
    #} else {

    if (N_TO_TEST == 1) # <-- Plots don't seem useful for ANOVA
        for(i in 1:N_TO_TEST) {
            if (N_TO_TEST == 1) {
                pngname = sprintf('%s.png', OUTPUT_PLOT)
            } else {
                pngname = sprintf('%s-%s.png', OUTPUT_PLOT, TERM_NAMES[i]) 
            }
            
            sane <- abs(result[,colnames(design)[i]]) < 1000.0
        
            png(pngname, width=800, height=800 )
            plot(result[sane,"log2 average per million"], 
                 result[sane,colnames(design)[i]],
                 xlab="log2 average per million",
                 ylab="log2 fold change",
                 main=TERM_NAMES[i],
                 pch=19, cex=0.25)
                 
            points(result[significant & sane,"log2 average per million"],
                   result[significant & sane,colnames(design)[i]],
                   pch=19, cex=1.0,
                   col="red")
        
            #points((result$"log2 average per million"),
            #       sqrt( dispersion[reordering] ),
            #       pch=19,
            #       col="blue")
            
            dev.off()
        }

    #}
    
    heatmap.features <- rev(as.character(result$Feature[significant]))
    heatmap.data <- glog2.rpm.counts(dgelist, GLOG_MODERATION)$E[heatmap.features,,drop=FALSE]
       #voom(dgelist)$E[heatmap.features,,drop=FALSE]
    heatmap.labels <- list(heatmap.features)
    heatmap.values <- result[heatmap.features,colnames(design)[1:N_TO_TEST],drop=FALSE]    
    
    #if (nrow(heatmap.values) > 0) { #Hnnggrrrakkkaahnaaarrrrg flrup
    #    for(i in basic.seq(ncol(heatmap.values))) {
    #        heatmap.values[ heatmap.values[,i] < -10.0, i ] <- -10.0
    #        heatmap.values[ heatmap.values[,i] > 10.0, i ] <- 10.0
    #    }
    #}
    
    for(i in basic.seq(ncol(dgelist$genes)-1)+1) {
        if ( colnames(dgelist$genes)[i] != "Length" &&
             colnames(dgelist$genes)[i] != "Total reads" &&
             colnames(dgelist$genes)[i] != "On same fragment" &&
             colnames(dgelist$genes)[i] != "Ambiguous alignment" )
            heatmap.labels[[length(heatmap.labels)+1]] <- dgelist$genes[heatmap.features,i]
    }
    do.heatmap(heatmap.features,heatmap.data,heatmap.labels,heatmap.values)
    
    
    sink(LOG_FILENAME, split=TRUE, append=TRUE)
    
    cat(dgelist$original.number.of.genes, 'features\n')
    cat('Discarded', dgelist$original.number.of.genes-nrow(dgelist$counts), 'features with total count less than', MIN_COUNT, '\n')
    cat('Kept', nrow(dgelist$counts), 'features\n')
    
    cat('\n')
    if (! QUANTILE_NORM) {
        print(d$samples)
    }
    
    cat('\nEstimated dispersion:\n')
    if (length(dispersion) > 1) {
        cat(describe(dispersion), '\n')
    } else {
        cat(dispersion, '\n')
    }
    
    cat('\n\nWith', describe.significant(), '\n')
    cat('\n',sum(significant), 'genes called as differentially expressed\n')
    
    #if (N_TO_TEST == 1) {
    #    cat(sum(result[significant,colnames(design)[1]] > 0.0), 'up\n')
    #    cat(sum(result[significant,colnames(design)[1]] < 0.0), 'down\n')
    #}
    
    cat("\nLimma's convest function estimates that the true number of tags", 
        "\nwhich are differentially expressed is around", (1.0-convest(result$p)) * length(result$p), "\n")
    
    cat('\nTotal counts:\n')
    cat(describe(rowSums(rcounts)), '\n')
    
    if (sum(significant)) {
        cat('\nTotal counts within significant tags:\n')
        cat(describe(rowSums(rcounts[significant,,drop=FALSE])), '\n')
    
        for(i in 1:N_TO_TEST) {
            cat('\nAbsolute log2', TERM_NAMES[i],'of significant tags:\n')
            cat(describe(abs(result[,colnames(design)[i]])), '\n')
        }
    }
}
"""



LIMMA = COMMON + r"""

library(limma)

null.design <- design[, (N_TO_TEST+1):ncol(design), drop=FALSE]

if (MODE == 'glog' || MODE == 'nullglog') {
    y <- glog2.rpm.counts(dgelist, GLOG_MODERATION)

    if (QUANTILE_NORM) {
        y$E <- normalizeBetweenArrays(y$E, method = 'quantile')
    }

} else {
    if (MODE == 'nullvoom' || MODE == 'nullglog') {
        voom.design <- null.design
    } else {
        voom.design <- design
    }
    
    if (QUANTILE_NORM) {
        norm <- 'quantile'
    } else {
        norm <- 'none'
    }
    
    png(sprintf('%s-voom.png', OUTPUT_PLOT), width=800, height=800 )
    y <- voom(dgelist, voom.design, plot=TRUE, normalize.method=norm)
    dev.off()
}

fit <- lmFit(y, design)

#if (USE_CONTRAST) {
#    fit <- contrasts.fit(fit, CONTRAST_WEIGHTS)
#    coef <- data.frame('log2 contrast' = fit$coefficients[,1], check.names=FALSE)
#} else {

coef <- fit$coefficients[,1:N_TO_TEST, drop=FALSE]

#}

if (MODE != 'nullvoom' && MODE != 'nullglog') {
    fit <- eBayes(fit)
    
    #if (USE_CONTRAST) {
    #    top <- topTable(fit, coef=1, sort.by='none', confint=TRUE, number=Inf)
    #} else {    
    
    top <- topTable(fit, coef=1:N_TO_TEST, sort.by='none', confint=N_TO_TEST==1, number=Inf)
    
    #}
     
    ordering <- order(top$P.Value)
    p <- top$P.Value

} else {
    #stopifnot(!USE_CONTRAST) #Not supported.

    fit <- lmFit(y, null.design)
    fit <- eBayes(fit)
    
    df <- ncol(dgelist$counts) + fit$df.prior - ncol(null.design)
    s2 <- fit$s2.post  #fit$sigma ^ 2
    
    #Trimmed chisq (untrimmed turns out to work just as well)
    #cut <- quantile(s2, 0.99)[[1]]
    #s2.good <- s2[ s2 < cut ]
    #n.good <- length(s2.good)        
    #f <- function(scale) {
    #    n.good*(log(pchisq(cut*scale,df))-log(scale)) -sum(dchisq(s2.good*scale,df,log=TRUE))
    #}
     
    #Untrimmed chisq, however there's an algebraic way to do this   
    #f <- function(scale) {
    #    -length(s2)*log(scale)-sum(dchisq(s2*scale,df,log=TRUE))
    #}    
    #search.max <- 2.0*df/median(s2)
    #scale <- optimize( f, interval=c(0.0, search.max) )$minimum
    # Was scale within expected bounds?
    #stopifnot(
    #    scale > search.max*0.01,
    #    scale < search.max*0.99
    #)    
    
    #Maximum likelihood scale of chi-square (special case of gamma distribution)
    #scale <- df*length(s2)/sum(s2)
    
    # fit$s2.prior allows that a proportion (1%) of genes are differentially expressed
    # Similar idea to trimmed chisq, above.
    scale <- df / fit$s2.prior

    png(sprintf('%s-qq.png', OUTPUT_PLOT), width=800, height=800 )
    qqplot(qchisq( ((1:length(s2))-0.5)/length(s2) , df)/scale, s2,
        xlab='Fitted chi-square quantiles',
        ylab='sigma^2',
        main='QQ plot of sigma^2 unexplained by null model vs chi-squared distribution' 
    ) 
    abline(0,1)
    dev.off()
    
    ordering <- order(s2, decreasing=TRUE) 
    p <- pchisq(s2*scale, df, lower.tail=FALSE)
}

result <- data.frame(Feature = rownames(y$genes), row.names=rownames(y$genes))

result$'average expression (log2 reads-per-million)' <- fit$Amean

for(i in basic.seq(ncol(coef))) {
    result[,colnames(coef)[i]] <- coef[,i]
}

if ((MODE == 'glog' || MODE == 'voom') && ncol(coef) == 1) {
    result[,'+/- at 95% confidence'] <- (top$CI.975-top$CI.025)*0.5
}

result$p <- p

result$FDR <- p.adjust(p, method='BH')

for(i in basic.seq(ncol(y$genes))) {
    result[,colnames(y$genes)[i]] <- y$genes[,i]
}

#------------------------------------------------------
result <- result[ordering,]

#significant <- (result$FDR <= FDR_CUTOFF)
significant <- get.significant(result)

results_to_output <- result

if (OUTPUT_COUNTS) {
    # Add raw counts
    results_to_output <- cbind(results_to_output, dgelist$counts[rownames(results_to_output),])
}

blank <- mapply(function(i){ all(is.na(results_to_output[,i])) }, 1:ncol(results_to_output))

sink(OUTPUT_FILENAME_ALL)
write.csv(results_to_output[,!blank,drop=FALSE], na='', row.names=FALSE)
sink()

sink(OUTPUT_FILENAME_SIGNIFICANT)
write.csv(results_to_output[significant,!blank,drop=FALSE], na='', row.names=FALSE)
sink()


if (ncol(coef) == 1) { # <-- Plots don't seem useful for ANOVA
    i <- 1    
    colname <- colnames(coef)[i]    
    
    #Strip "log2 "
    name <- substr(colname, 6,nchar(colname))
    
    pngname = sprintf('%s.png', OUTPUT_PLOT)

    png(pngname, width=800, height=800 )

    plot(result[,"average expression (log2 reads-per-million)"], 
         result[,colname],
         xlab="log2 concentration",
         ylab="log2 fold change",
         main=name,
         pch=19, cex=0.25)
    
    #if (MODE == 'glog' || MODE == 'voom') {
    #    title(xlab='Error bars indicate 95% confidence interval',adj=1)
    #    
    #    use <- abs(result[,colname]) >= sort(abs(result[,colname]))[floor(1+nrow(result) * 0.99)]
    #    arrows(result[,"average expression (log2 reads-per-million)"][use],
    #           y0=top[ordering,'CI.025'][use],
    #           y1=top[ordering,'CI.975'][use],
    #           angle=90,code=3,length=0.03,col=c('#8888ff'))  
    #    points(result[,"average expression (log2 reads-per-million)"], 
    #         result[,colname],
    #         pch=19, cex=0.25)
    #}
        
    points(result[,"average expression (log2 reads-per-million)"][significant],
           result[,colname][significant],
           pch=19, cex=1.0,
           col="red")
    dev.off()
}


heatmap.features <- rev(as.character(result$Feature[significant]))
heatmap.data <- y$E[heatmap.features,,drop=FALSE]
heatmap.labels <- list(heatmap.features)
heatmap.values <- result[heatmap.features,colnames(coef),drop=FALSE]
for(i in basic.seq(ncol(y$genes)-1)+1) {
    if ( colnames(y$genes)[i] != "Length" &&
         colnames(y$genes)[i] != "Total reads" &&
         colnames(y$genes)[i] != "On same fragment" &&
         colnames(y$genes)[i] != "Ambiguous alignment" )
        heatmap.labels[[length(heatmap.labels)+1]] <- y$genes[heatmap.features,i]
}
do.heatmap(heatmap.features,heatmap.data,heatmap.labels,heatmap.values)

sink(LOG_FILENAME, split=TRUE, append=TRUE)

cat('                   Sample    Library size   Further Trimmed-Mean normalization\n')
for(i in 1:ncol(dgelist$counts)) {
    cat(sprintf("%25s %15d       %f\n", rownames(dgelist$samples)[i], dgelist$samples$lib.size[i], dgelist$samples$norm.factors[i]))
}
cat('\n\n')

if (MODE == 'glog' || MODE == 'nullglog') {
    cat('Using glog_moderation =',GLOG_MODERATION,'\n\n')
}

cat(dgelist$original.number.of.genes, 'features\n')
cat('Discarded', dgelist$original.number.of.genes-nrow(dgelist$counts), 'features with total count less than', MIN_COUNT, '\n')
cat('Kept', nrow(dgelist$counts), 'features\n')

#if (MODE == 'voom' || MODE == 'glog') {
cat(fit$df.prior, 'prior df (the prior is like this many extra samples)\n')
#}

cat('\nWith', describe.significant(), '\n\n')
cat(' ', sum(significant), 'genes called as differentially expressed\n\n')

cat("Limma's convest function estimates that the true number of tags\n")
cat("which are differentially expressed is around", (1.0-convest(p)) * length(p), "\n\n")

"""



GOSEQ = r"""

library(nesoni)
library(goseq)

report <- function(result,col, significant,direction,categories) {
    good <- result
    good$fdr <- p.adjust(good[,col], method='BH')    
    good <- good[order(good[,col]),]
    
    good <- good[good$fdr < 0.5,,drop=FALSE]
    n <- 0
    while(n < nrow(good) && 
          ((n < 10 && good$fdr[n+1] < 1.0) || good$fdr[n+1] <= FDR_CUTOFF))
      n <- n + 1
    good <- good[basic.seq(n),,drop=FALSE]
    
    significants <- names(significant)[significant]
    
    all.genes <- list()
    for(i in basic.seq(nrow(good))) {
        all.genes[[i]] <- sort(unique(categories[categories[,2] == good[i,'category'],1]))
    }
    
    seen <- rep(FALSE, nrow(good))    
    
    for(i in basic.seq(nrow(good))) {
        if (!seen[i]) {
            for(j in (i:nrow(good))) {
                if (identical(all.genes[[i]],all.genes[[j]])) {
                    cat(if (good[j,'fdr'] > FDR_CUTOFF) '(' else ' ')
                    cat(sprintf('FDR %-10.2g  %s',good[j,'fdr'],good[j,'category']))
                    cat(cat(if (good[j,'fdr'] > FDR_CUTOFF) ')\n' else '\n'))
                    seen[j] <- TRUE
                }
            }
            
            genes <- all.genes[[i]]
            siggenes <- genes[genes %in% significants]
            cat(sprintf(' %-16s',sprintf('%d of %d',length(siggenes),length(genes))))        
            for(gene in siggenes) {
                    if (direction[gene] < 0) cat('-')
                    if (direction[gene] > 0) cat('+')
                    cat(sprintf('%s ',gene))
            }
            cat('\n\n')
        }
    }
    
    if (nrow(good) == 0) cat('none\n\n')
}

analyse <- function(significant,direction,length,categories) {
    pwf <- tryCatch({
            nullp(significant, bias.data=length, plot.fit=FALSE)
        }, error=function(err) { 
            NULL 
        })
    
    if (is.null(pwf)) {
        cat('\n\ngoseq function "nullp" failed.\n\n')
    } else {
        result <- goseq(pwf,gene2cat=categories)
        cat('\nOver-represented categories:\n\n')
        report(result,'over_represented_pvalue', significant,direction,categories)
        cat('\nUnder-represented categories:\n\n')
        report(result,'under_represented_pvalue', significant,direction,categories)    
    }
}


data <- read.grouped.table(INPUT_FILENAME_ALL)$All
data.significant <- read.grouped.table(INPUT_FILENAME_SIGNIFICANT)$All
categories <- read.csv(CATEGORY_FILENAME)

#significant <- data$FDR <= FDR_CUTOFF
significant <- rownames(data) %in% rownames(data.significant)

if (substr(colnames(data)[2],1,5) == 'log2 ' && substr(colnames(data)[3],1,5) != 'log2 ') {
    direction <- sign(data[,2])
} else {
    direction <- rep(0,nrow(data))
}
names(direction) <- rownames(data)

names(significant) <- rownames(data)
length <- data$Length

sink(LOG_FILENAME, split=TRUE, append=TRUE)
cat('\n\n=============================\n')
cat(NAME)
cat('\n\n')

#cat(sprintf('%d of %d significantly DE genes\n',sum(significant),length(significant)))

analyse(significant,direction,length,categories)

"""



def term_specification(term):
    if '=' not in term: return term
    return term.split('=',1)[0]

def term_name(term):
    if '=' not in term: return term
    return term.split('=',1)[1]


@config.help("""\
Find significant differential expression from the output of "count:".

Terms are selection expressions on the sample names and any tags they were \
given with "tag:". See the main nesoni help page for selection expression syntax. \
"with:" terms are included in the linear model, but not the significance test. \
The model will also include a constant term.

A term can be given a more intuitive name in the output by appending "=name". \
Examples: A=fold-change A^B=interaction-between-A-and-B


Example 1: Say we have two conditions, "control" and "experimental". Here's are some possible tests:

experimental
- What is the fold change from control to experimental? \
Hypothesis 0 is that all samples have the same log expression level. \
Hypothesis 1 is that the control samples have a baseline log expression level, \
and the experimental samples have that baseline plus a further amount.


Example 2: Say we have two strains, "strain1" and "strain2" and two time points "time1" and "time2". \
Here are some possible tests:

time2
- What is the fold change from time1 to time2? \
But note that any difference between strains will be seen as within group variation!

time2 --select strain1
- What is the fold change from time1 to time2, looking at strain1 samples only?

time2 with: strain2
- What is there a change from time1 to time2? \
Any difference between strains is accounted for.

strain2^time2 with: strain2 time2
- Is there an interaction effect between strain and time?


Example 3: Say we have three time points, "time1", "time2" and "time3". Here are some possible tests:

time3 --select time2/time3
- Is there a change between time2 and time3?

time2 time3
- Is there any change over time? \
This is an ANOVA style test that is potentially more sensitive than comparing each pair of time points individually.


""" +
''.join(
        '\n\n'+
        '='*20 +
        '\n--mode %s \n' % mode +
        MODE_HELP[mode]        
        for mode in ('voom', 'glog', 'nullvoom', 'nullglog', 'poisson', 'common', 'trend')
))
@config.String_flag('mode', 'Analysis mode. See below.')
@config.Float_flag('glog_moderation', 
    'For "--mode glog" only.'
    ' Amount of moderation used in log transformation.'
    )
@config.String_flag('select', 'Selection expression. Which samples to use.')
@config.Int_flag('min_count', 'Discard features with less than this total count.')
@config.Float_flag('fdr', 'False Discovery Rate cutoff.')
@config.Float_flag('log_fold', 
     'Absolute log2 fold change cutoff. '
     'If there are several test terms, at least one must exceed the given fold change.')
@config.Bool_flag('output_counts', 'Include raw read counts in output.')
@config.Bool_flag('constant_term', 'Include a constant term in the model.')
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.Bool_flag('tell', 'Output R code instead of executing it.')
@config.Positional('counts_file', 'The table of counts produced by "count:".')
@config.Main_section('test', 'Terms to test.')
#@config.Float_section('contrast', 'Instead of doing an ANOVA on multiple terms, perform a contrast with the given weights.')
@config.Section('with_', 'Also include these terms in the model.')
#@config.Section('use', 'Only use samples matching these regexes.')
@config.Section('goseq', 
     'One or more files containing gene categories. '
     'goseq will be used to find over-represented categories in the significantly differentially expressed genes.\n\n'
     'A category file should be a CSV file. The first line should contain the headings:\n\n'
     'Gene,Category\n\n'
     'Subsequent lines should give pairings of genes with categories.\n\n'
     'The filename may be given as "filename.csv=heading", in which case the heading will be used in place of the filename in the output.'
     )
class Test_counts(config.Action_with_prefix):
    select = 'all'
    min_count = 10
    constant_term = True
    mode = 'voom'
    glog_moderation = 5.0
    quantile_norm = False
    fdr = 0.01
    log_fold = 0.0
    output_counts = False
    tell = False
    counts_file = None
    norm_file = None
    
    test = [ ]
    #contrast = [ ]
    with_ = [ ]
    use = [ ] 
    goseq = [ ]

    def run(self):
        min_count = self.min_count
        constant_term = self.constant_term
        mode = self.mode
        glog_moderation = self.glog_moderation
        quantile_norm = self.quantile_norm
        fdr = self.fdr
        log_fold = self.log_fold
        only_tell = self.tell
        output_prefix = self.prefix
        filename = self.counts_file
        test_terms = self.test
        with_terms = self.with_
        #contrast_weights = self.contrast
        select = self.select
        norm_file = self.norm_file
        goseq = self.goseq
        output_counts = self.output_counts
    
        log = grace.Log()
    
        assert mode in ('voom', 'glog', 'nullvoom', 'nullglog', 'poisson', 'common', 'trend')
        
        assert test_terms, 'No terms to test'
        
        #use_contrast = bool(contrast_weights)
        #assert not use_contrast or len(contrast_weights) == len(test_terms), 'Need one contrast weight per term'
        
        all_terms = test_terms + with_terms 
        if constant_term: all_terms += [ 'all' ]
        n_to_test = len(test_terms)
        
        #if use_contrast:
        #    contrast_weights += [ 0.0 ] * (len(all_terms)-len(test_terms))
        
        reader = io.Table_reader(filename, 'Count')
        reader.close()
        all_samples = [ item for i, item in enumerate(reader.headings) if reader.groups[i] == 'Count' ]
        tags = { }
        for item in all_samples:
            tags[item] = [ item ]
        for line in reader.comments:
            if line.startswith('#sampleTags='):
                parts = line[len('#sampleTags='):].split(',')
                tags[parts[0]] = parts
        
        n_all_samples = len(all_samples)
            
        keep = [ selection.matches(select, tags[item]) for item in all_samples ]
        samples = [ all_samples[i] for i in xrange(n_all_samples) if keep[i] ]
        n_samples = len(samples)
        
        assert n_samples, 'No samples selected'
        assert n_samples > 1, 'Only one sample selected'
        
        model = [ ]
        for term in all_terms:        
            spec = term_specification(term)
            model.append([ 1 if selection.matches(spec, tags[item]) else 0 for item in samples ])
        model = zip(*model) #Transpose
            
        log.log('Model matrix (terms being tested for significance shown in { }):\n')
        for i in xrange(n_samples):
            log.log('{ ')
            j = 0
            while True:
                if j == n_to_test:
                    log.log('} ')
                if j >= len(all_terms): 
                    break
                log.log('%d ' % model[i][j])
                j += 1
            log.log(samples[i] + '\n')
        log.log('\n')    
        
        
        assert len(all_terms) <= n_all_samples, 'Can\'t have more linear model terms than samples.'
        if mode not in ('poisson','nullvoom','nullglog'):
            assert len(all_terms) < n_all_samples, (
                'Can\'t have as many linear model terms as samples, within group variation can\'t be estimated. '
                'I wouldn\'t recommend using --mode poisson or --mode nullvoom, '
                'you should go and sequence some more samples.'
            )
        
        
        log.log('Analysis mode: %s\n%s\n'%(mode,MODE_HELP[mode]))
        
        log_filename = output_prefix + '-info.txt'
        design_columns = [ 'log2 '+term_name(item) for item in all_terms ]
        term_names = [ term_name(item) for item in all_terms ]
        
        log.attach(open(log_filename,'wb'))
        log.close()    
        
        if mode in ('voom', 'glog', 'nullvoom', 'nullglog'):
            script = LIMMA
        else:
            script = EDGER
        
        if os.path.exists(output_prefix+'.png'):
            os.unlink(output_prefix+'.png')
        
        run_script(
            script,
            
            FILENAME = filename,
            N_SAMPLES = n_samples,
            N_ALL_SAMPLES = n_all_samples,
            N_TO_TEST = n_to_test,
            FDR_CUTOFF = fdr,
            LOG_FOLD_CUTOFF = log_fold,
            OUTPUT_COUNTS = output_counts,
            MIN_COUNT = min_count,
            MODE = mode,
            GLOG_MODERATION = glog_moderation,
            QUANTILE_NORM = quantile_norm,
            KEEP = keep,
            OUTPUT_FILENAME_ALL = output_prefix + '-all.csv',
            OUTPUT_FILENAME_SIGNIFICANT = output_prefix + '.csv',
            LOG_FILENAME = log_filename,
            OUTPUT_PLOT = output_prefix,
            #USE_CONTRAST = use_contrast,
            #CONTRAST_WEIGHTS = contrast_weights,
            DESIGN = model,
            DESIGN_COLUMNS = design_columns,
            TERM_NAMES = term_names,
            NORM_FILE = norm_file,
            
            only_tell=only_tell
            )
        
        for item in goseq:
            run_script(
                GOSEQ,
                
                INPUT_FILENAME_ALL = output_prefix + '-all.csv',
                INPUT_FILENAME_SIGNIFICANT = output_prefix + '.csv',
                LOG_FILENAME = log_filename,
                NAME = term_name(item),
                CATEGORY_FILENAME = term_specification(item),
                FDR_CUTOFF = fdr,
            
                only_tell=only_tell
                )

    def report(self, reporter):
        prefix = os.path.basename(self.prefix)
        image = reporter.get(self.prefix + '-heatmap.png', image=True, title='')
        sig_csv = reporter.get(self.prefix + '.csv', title='[DE genes table]')
        all_csv = reporter.get(self.prefix + '-all.csv', title='[All genes table]')
        maybe_maplot = (
            (' &sdot; ' + reporter.get(self.prefix+'.png', image=True,title='[MA-plot]'))
            if os.path.exists(self.prefix+'.png') else ''
            )
        info = reporter.get(self.prefix + '-info.txt', title='[Info]')        
        
        text = (
            '<table><tr>\n'
            '<td valign="top">%(image)s</td>\n'
            '<td valign="top"><b>%(prefix)s</b>\n'
            '<br/>%(sig_csv)s &sdot; %(all_csv)s &sdot; %(info)s %(maybe_maplot)s\n'
            '</tr></table>'
            ) % locals()
        
        reporter.write(text)
            


#TEST_POWER_HELP = """\
#
#Usage:
#
#    nesoni test-power: [options] output_prefix [of: [test-count options]] 
#
#Test the statistical power of "nesoni test-counts".
#
#Options:
#
#    --m NNN            - Number of differentially expressed tags
#    
#    --n NNN            - Total number of tags
#    
#    --reps N           - Within group replicates
#    
#    --count NNNN       - Average count
#    
#    --dispersion N.NN  - Dispersion
#    
#    --log-fold N.NN    - log2 Fold change of differentially expressed tags
#
#"""
#
#POWER_TEMPLATE = r"""
#
#library(MASS)
#
#m <- %(m)d
#n <- %(n)d
#reps <- %(reps)d
#avg_count <- %(count)d
#dispersion <- %(dispersion)f
#log_fold <- %(log_fold)f
#
#mynegbin <- function(u, d) {
#    rnegbin(1, u, 1.0/d)
#}
#
#counts <- matrix(0, nrow=n, ncol=(reps*2) )
#names <- rep('', n)
#
#for(i in 1:n) {
#    a <- 1.0
#    b <- 1.0
#    if (i <= m) {
#        b <- b * (2 ** sample(c(log_fold,-log_fold),1))
#        names[i] <- paste('DifferentialTag',i,sep='')
#    } else {
#        s <- 0.0
#        names[i] <- paste('Tag',i,sep='')
#    }
#    scale <- avg_count / (0.5*(a+b))
#    a <- a * scale
#    b <- b * scale
#    
#    for(j in 1:reps) {
#        counts[i,j]      <- mynegbin(a, dispersion)
#        counts[i,reps+j] <- mynegbin(b, dispersion)
#    }
#}
#
#data <- data.frame(
#    Feature=names,
#    counts,
#    counts,
#    Product=names
#)
#
#for(i in 1:reps) {
#   colnames(data)[i+1] <- paste('Control',i,sep='')
#   colnames(data)[reps+i+1] <- paste('Experimental',i,sep='')
#}
#for(i in 1:(reps*2)) {
#   colnames(data)[i+reps*2+1] <- paste('RPKM', colnames(data)[i+1])
#}
#
#write.table(data, %(filename_literal)s, sep='\t', quote=FALSE, row.names=FALSE)
#
#"""
#
#POWER_REPORT_TEMPLATE = """
#
#claimed_fdr <- %(claimed_fdr)f
#
#data <- read.delim(%(output_filename_literal)s, check.names=FALSE)
#
#false_positives <- 0
#true_positives <- 0
#
#for(i in 1:nrow(data)) {
#    sig <- data$FDR[i] <= claimed_fdr
#    true <- regexpr('Differential', data$Feature[i])[[1]] != -1
#    
#    if (sig && !true) { false_positives <- false_positives + 1 }
#    if (sig && true) { true_positives <- true_positives + 1 }
#}
#
#sink(%(log_filename_literal)s, append=TRUE, split=TRUE)
#
#cat('\n\nFalse positives:', false_positives, '\n')
#cat(    'False negatives:', %(m)d - true_positives, '\n')
#cat(    'True positives:', true_positives, '\n')
#
#"""
#
#def test_power_main(args):
#    m, args = grace.get_option_value(args, '--m', int, 10)
#    n, args = grace.get_option_value(args, '--n', int, 1000)
#    reps, args = grace.get_option_value(args, '--reps', int, 2)
#    count, args = grace.get_option_value(args, '--count', int, 100)
#    dispersion, args = grace.get_option_value(args, '--dispersion', float, 0.1)
#    log_fold, args = grace.get_option_value(args, '--log-fold', float, 1.0)
#
#    if len(args) < 1:
#        print >> sys.stderr, TEST_POWER_HELP
#        raise grace.Help_shown()
#    
#    output_prefix, args = args[0], args[1:]
#
#    options = [ ]
#    def of(args):
#        options.extend(args)    
#    grace.execute(args, {'of': of})
#    
#    filename = output_prefix + '-input.txt'
#    filename_literal = R_literal(filename)
#    log_filename_literal = R_literal(output_prefix + '-info.txt')
#    
#    run_script(POWER_TEMPLATE % locals())
#
#    claimed_fdr = test_counts_main([ output_prefix, filename, 'Experimental' ] + options)
#
#    output_filename_literal = R_literal(output_prefix + '.txt')
#
#    run_script(POWER_REPORT_TEMPLATE % locals())

POWER_SCRIPT = """
library(nesoni)

library(MASS)
mynegbin <- function(u, d) {
    rnegbin(1, u, 1.0/max(1e-30,d))
}

q <- c(runif(DE, 1.0-DE_TOP, 1.0), runif(GENES-DE, 0.0, 1.0))
elevels <- qnorm(q, 0.0, GENE_SD)
dispersions <- rnorm(GENES, DISPERSION, DISPERSION_SD)

e <- matrix(0, nrow=GENES, ncol=REPS*2)
rnames <- rep('', GENES)
cnames <- rep('', REPS*2)

for(i in basic.seq(REPS)) {
    cnames[i] <- sprintf('control%d', i)
    cnames[i+REPS] <- sprintf('experimental%d',i)
}

for(i in basic.seq(REPS*2)) {
    e[,i] <- elevels
}

for(i in basic.seq(GENES)) {
    if (i < DE) {
        rnames[i] <- sprintf('Diff%d',i)
        diff <- DE_AMOUNT * sample(c(-1,1),1)
        for(j in basic.seq(REPS)+REPS)
            e[i,j] <- e[i,j] + diff
    } else {
        rnames[i] <- sprintf('Same%d',i)
    }
}

rownames(e) <- rnames
colnames(e) <- cnames

counts <- 2**e 
for(i in basic.seq(REPS*2)) {
    scale <- READS / sum(counts[,i])
    for(j in basic.seq(GENES)) {
        counts[j,i] <- mynegbin(counts[j,i]*scale, dispersions[j])
    }
}

annotation = data.frame(gene=rnames,row.names=rnames)

write.grouped.table(list(Count=data.frame(counts),Annotation=annotation), FILENAME)
"""

@config.Int_flag('genes', 'Number of genes')
@config.Int_flag('reads', 'Number of reads per sample')
@config.Int_flag('reps', 'Replicates in each group')
@config.Float_flag('gene_sd', 'Standard deviation of log2 average expression of each gene.')
@config.Float_flag('dispersion', 'Mean dispersion (biological/experimental noise level)')
@config.Float_flag('dispersion_sd', 'Standard deviation of dispersion.')
@config.Int_flag('de', 'Number of differentially expressed genes.')
@config.Float_flag('de_amount', 'Log2 fold change of differentially expressed genes.')
@config.Float_flag('de_top', 'Differentially expressed genes have (measured) expression levels in the top this proportion of genes.')
@config.Int_flag('repeats', 'How many trials of "test-counts:" to perform.')
@config.Configurable_section('test_counts', 'Options for "test-counts:"', presets=[
    ('default', lambda obj: Test_counts(), '')
    ])
class Test_power(config.Action_with_prefix):
    genes = 3000
    reads = 1000000
    reps = 3
    gene_sd = 3.0    
    dispersion = 0.1
    dispersion_sd = 0.0
    de = 100
    de_amount = 2.0
    de_top = 1.0
    repeats = 10    
    #test_counts=

    def run(self):
        assert self.repeats >= 1
        assert self.reps >= 1
        assert self.de <= self.genes

        futures = [ legion.future(self._run_test) for i in xrange(self.repeats) ]
        good = sum(item()[0] for item in futures)
        bad = sum(item()[1] for item in futures)

        self.log.log('\n')
        self.log.log(' True positives: %4d\n' % good)
        self.log.log('False positives: %4d (expected %.1f)\n' % (bad, self.test_counts.fdr*(good+bad)))
        self.log.log('False negatives: %4d\n' % (self.de*self.repeats-good))
    
    def _run_test(self):
        good = 0
        bad = 0
        with workspace.tempspace() as space:
            run_script(POWER_SCRIPT,
                FILENAME=space/'counts.csv',
                GENES=self.genes, READS=self.reads, REPS=self.reps,
                GENE_SD=self.gene_sd, DISPERSION=self.dispersion, DISPERSION_SD=self.dispersion_sd,
                DE=self.de, DE_AMOUNT=self.de_amount, DE_TOP=self.de_top,
                )
            self.test_counts(
                space/'test',
                space/'counts.csv',
                test=['experimental'],
                output_all=False
                ).run()
            
            with open(space/'test.txt', 'rb') as f:
                f.readline()
                for line in f:
                    if line.startswith('Diff'):
                        good += 1
                    else:
                        bad += 1
        return good, bad


HEATMAP_SCRIPT = """

library(nesoni)

dgelist <- read.counts(COUNTS, min.total=MIN_TOTAL, min.max=MIN_MAX, norm.file=NORM_FILE)

if (length(ORDER)) {
    dgelist <- dgelist[, ORDER]
}

#elist <- voom(dgelist, normalize.method=if(QUANTILE) 'quantile' else 'none')
elist <- glog2.rpm.counts(dgelist, GLOG_MODERATION)

hmap.elist(PREFIX, elist, min.sd=MIN_SD, min.span=MIN_SPAN, min.svd=MIN_SVD, svd.rank=SVD_RANK, reorder.columns=REORDER_COLUMNS)

"""


@config.help("""\
Produce a heatmap (and table) of log expression levels.

The parameter --glog-moderation controls how counts close to zero are treated. \
Counts close to zero are inherently noisier than larger counts. \
--glog-moderation controls the degree to which this variation is squashed down.

Hierachical clustering and ordering of rows is performed using the "seriation" package in R+.

You will need to use at least one of --min-sd, --min-span, or --min-svd.
""")
@config.Float_flag('glog_moderation', 
    'Amount of moderation used in log transformation.'
    ' See "glog" mode in "test-counts:".'
    )
@config.Int_flag('min_total', 'Expression level filter:\nExclude genes with less than this total number of reads')
@config.Int_flag('min_max', 'Expression level filter:\nExclude genes with no sample having at least this many reads')
@config.Float_flag('min_sd', 'Fold change filter:\nExclude genes with less than this standard deviation of log2 expression levels')
@config.Float_flag('min_span', 'Fold change filter:\nExclude genes where there is no pair of samples that differ in log2 expression by this much')
@config.Float_flag('min_svd', 'Fold change filter:\nUsing the SVD of log2 expression levels, '
                              'exclude genes where the sum of squares of the U matrix row for the gene is less than this many standard deviations.')
@config.Int_flag('svd_rank', 'Only use the top this many columns of the SVD U matrix when applying --min-svd.')
@config.Bool_flag('reorder_columns', 'Cluster and optimally order columns as well as rows.')
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
#@config.Bool_flag('quantile', 'Use quantile normalization.')
@config.Positional('counts', 'File containing output from "nesoni count:"')
@config.String_flag('select', 'Selection expression (see main help text). Which samples to use.')
@config.String_flag('sort', 'A sort expression to sort samples (see main help text).')
#@config.Section('order', 'Optionally, specify an order to show the columns in.')
class Heatmap(config.Action_with_prefix):
    glog_moderation = 5.0
    counts = None
    min_total = 0
    min_max = 0
    min_sd = 0.0
    min_span = 0.0
    min_svd = 0.0
    svd_rank = None
    reorder_columns = False
    norm_file = None
    #quantile = False
    select = 'all'
    sort = ''
    #order = [ ]

    def run(self):
        reader = io.Table_reader(self.counts, 'Count')
        reader.close()
        all_samples = [ item for i, item in enumerate(reader.headings) if reader.groups[i] == 'Count' ]
        tags = { }
        for item in all_samples:
            tags[item] = [ item ]
        for line in reader.comments:
            if line.startswith('#sampleTags='):
                parts = line[len('#sampleTags='):].split(',')
                tags[parts[0]] = parts
            
        samples = selection.select_and_sort(
            self.select, self.sort, all_samples, lambda sample: tags[sample])

        run_script(HEATMAP_SCRIPT,
            PREFIX=self.prefix,
            GLOG_MODERATION=self.glog_moderation,
            COUNTS=self.counts,
            MIN_TOTAL=self.min_total,
            MIN_MAX=self.min_max,
            MIN_SD=self.min_sd,
            MIN_SPAN=self.min_span,
            MIN_SVD=self.min_svd,
            SVD_RANK=self.svd_rank,
            REORDER_COLUMNS=self.reorder_columns,
            NORM_FILE=self.norm_file,
            #QUANTILE=self.quantile,
            ORDER=samples,
        )




SIMILARITY_SCRIPT = """

library(nesoni)

dgelist <- read.counts(COUNTS, norm.file=NORM_FILE)

if (length(ORDER)) {
    dgelist <- dgelist[, ORDER]
}

elist <- glog2.rpm.counts(dgelist, GLOG_MODERATION)
mat <- elist$E

sink(sprintf("%s.nex", PREFIX))
cat("#NEXUS\n")
cat("begin taxa;\n")
cat(sprintf("dimensions ntax=%d;\n",ncol(mat)))
cat("taxlabels\n")
for(i in basic.seq(ncol(mat))) {
    cat(sprintf("%s\n", colnames(mat)[i]));
}
cat(";\n")
cat("end;\n")
cat("begin distances;\n")
cat("format triangle=lower diagonal nolabels;\n")
cat("matrix\n")
for(i in basic.seq(ncol(mat))) {
    for(j in basic.seq(i)) {
        diff <- mat[,j]-mat[,i]
        dist <- sqrt(mean(diff*diff))
        cat(sprintf("%f ",dist));
    }
    cat("\n")
}
cat(";\n")
cat("end;\n")
sink()

library(limma)
png(sprintf("%s-plotMDS.png",PREFIX),width=800,height=800)
plotMDS(elist)
dev.off()

"""


@config.help("""\
The similarity/difference between samples is visualized in several ways:

- The plotMDS function from limma.

- A SplitsTree4 generated NeighborNet based on root mean square difference in glog2 gene expression level between each pair of samples.

These plots may be helpful in identifying contaminated or mislabeled samples.
""")
@config.Float_flag('glog_moderation', 
    'Amount of moderation used in log transformation.'
    ' See "glog" mode in "test-counts:".'
    )
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.Positional('counts', 'File containing output from "nesoni count:"')
@config.String_flag('select', 'Selection expression (see main help text). Which samples to use.')
#@config.String_flag('sort', 'A sort expression to sort samples (see main help text).')
class Similarity(config.Action_with_prefix):
    glog_moderation = 5.0
    counts = None
    norm_file = None
    select = 'all'
    sort = ''

    def run(self):
        reader = io.Table_reader(self.counts, 'Count')
        reader.close()
        all_samples = [ item for i, item in enumerate(reader.headings) if reader.groups[i] == 'Count' ]
        tags = { }
        for item in all_samples:
            tags[item] = [ item ]
        for line in reader.comments:
            if line.startswith('#sampleTags='):
                parts = line[len('#sampleTags='):].split(',')
                tags[parts[0]] = parts
            
        samples = selection.select_and_sort(
            self.select, self.sort, all_samples, lambda sample: tags[sample])

        run_script(SIMILARITY_SCRIPT,
            PREFIX=self.prefix,
            GLOG_MODERATION=self.glog_moderation,
            COUNTS=self.counts,
            NORM_FILE=self.norm_file,
            ORDER=samples,
        )

        io.execute(
            'SplitsTree +g -i INPUT -x COMMAND',
            INPUT=self.prefix + '.nex',
            COMMAND='UPDATE; '
                    'SAVE FILE=\'%s.nex\' REPLACE=yes; '
                    'EXPORTGRAPHICS format=svg file=\'%s.svg\' REPLACE=yes TITLE=\'NeighborNet of expression levels\'; ' 
                    'QUIT' 
                    % (self.prefix, self.prefix),
            )

    def report(self, reporter):    
        reporter.heading('Sample similarity')
        
        reporter.p(
            'The following plots attempt to summarize the similarity/differences in expression patterns between samples, '
            'based on the glog2-transformed normalized read counts. '
            'Samples from the same experimental group should cluster together.'
            )
        
        reporter.p(
            reporter.get(self.prefix + '-plotMDS.png',
                title = 'limma\'s "plotMDS" Multi-Dimensional Scaling plot of sample similarity',
                image = True
                )
            )
        
        reporter.p(
            reporter.get(self.prefix + '.svg',
                title = 'Split Network visualization of sample similarity.',
                image = True
                ) +
            '<br>(Visualization of euclidean distances as a split network. '
            'Note: This is <i>not</i> a phylogenetic network.)'
            )



COMPARE_TESTS_SCRIPT = """

library(nesoni)

compare.tests(
    PREFIX,
    read.csv(A),
    read.csv(B),
    n=N,
    title=TITLE,
    a.title=a_TITLE,
    b.title=b_TITLE 
)

"""

@config.help("""\
Compare two tables produced by "test-counts:".

Produces an image comparing the order of the top differentially expressed genes.
""")
@config.Int_flag('n', 'Show this many genes')
@config.String_flag('title', 'Image title')
@config.String_flag('a_title', 'File A title')
@config.String_flag('b_title', 'File B title')
@config.Positional('a', 'Input file A')
@config.Positional('b', 'Input file B')
class Compare_tests(config.Action_with_prefix):
    n = 50
    title = ''
    a_title = ''
    b_title = ''
    a = None
    b = None
    
    def run(self):
        run_script(COMPARE_TESTS_SCRIPT,
            PREFIX=self.prefix,
            A=self.a,
            B=self.b,
            N=self.n,
            TITLE=self.title,
            a_TITLE=self.a_title,
            b_TITLE=self.b_title,
        )


NORMALIZATION_SCRIPT = """

library(nesoni)

dgelist <- read.counts(COUNTS_FILENAME, use.tmm=USE_TMM)

result <- data.frame(
    Normalizing.multiplier = dgelist$samples$normalizing.multiplier,
    row.names = rownames(dgelist$samples)
)

write.grouped.table(list(All=result), sprintf("%s.csv",PREFIX),comments=c('#Normalization'))

library(limma)

n <- nrow(result)

png(sprintf("%s-raw.png",PREFIX), width=500, height=100+20*n)
par(mar=c(5,10,4,2))
boxplot( voom(dgelist,lib.size=1e6)$E[,n:1], horizontal=TRUE, las=1, 
    main='Unnormalized',
    xlab='log2 reads (voom)')
dev.off()

png(sprintf("%s-libsize.png",PREFIX), width=500, height=100+20*n)
par(mar=c(5,10,4,2))
boxplot( voom(dgelist, lib.size=dgelist$samples$lib.size)$E[,n:1], horizontal=TRUE, las=1, 
    main='Normalization by library size\n(= total count of reads aligning to genes)',
    xlab='log2 reads-per-million (voom)')
dev.off()

if (USE_TMM) {
    png(sprintf("%s-norm.png",PREFIX), width=500, height=100+20*n)
    par(mar=c(5,10,4,2))
    boxplot( voom(dgelist)$E[,n:1], horizontal=TRUE, las=1, 
        main='Normalization by TMM effective library size',
        xlab='log2 reads-per-million (voom)')
    dev.off()
}

"""

@config.help("""\
Creates <prefix>.csv containing a table of \
normalizing multipliers as calculated by EdgeR using Trimmed Mean of M-values normalization.
""")
@config.Bool_flag('tmm', 'Use TMM normalization (in addition to normalizing by number of mapped reads).')
@config.Positional('counts_filename', 'CSV file containing counts.')
class Norm_from_counts(config.Action_with_prefix):
    tmm = True
    counts_filename = None

    def run(self):
        assert self.counts_filename, 'Counts filename required.'
        run_script(NORMALIZATION_SCRIPT,
            PREFIX=self.prefix,
            COUNTS_FILENAME=self.counts_filename,
            USE_TMM=self.tmm,
        )



GLOG_SCRIPT = r"""
library(nesoni)
dgelist <- read.counts(COUNTS, norm.file=NORM_FILE)
elist <- glog2.rpm.counts(dgelist, GLOG_MODERATION)

write.grouped.table(list(glog2=data.frame(elist$E), Annotation=elist$gene), sprintf("%s.csv",PREFIX),comments=c('#glog2 reads-per-million expression levels'))
"""

@config.Float_flag('glog_moderation', 
    'Amount of moderation used in log transformation.'
    ' See "glog" mode in "test-counts:".'
    )
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.Positional('counts', 'File containing output from "nesoni count:"')
class Glog(config.Action_with_prefix):
    prefix = None
    glog_moderation = 5.0
    norm_file = None
    counts = None

    def run(self):
        run_script(GLOG_SCRIPT,
            PREFIX=self.prefix,
            GLOG_MODERATION=self.glog_moderation,
            NORM_FILE=self.norm_file,
            COUNTS=self.counts,
            )



NMF_SCRIPT = """

library(nesoni)

if (is.null(ORDER_HINT)) {
    order.hint <- NULL
} else {
    order.hint <- readRDS(ORDER_HINT)
}

if (is.null(GLYPH)) {
    glyph <- NULL
} else {
    glyph <- eval(parse(text=GLYPH))
}

if (JUST_REPORT) {
    factorization <- readRDS(sprintf('%s.rds', PREFIX))
} else {

    dgelist <- read.counts(COUNTS, min.total=MIN_TOTAL, min.max=MIN_MAX, norm.file=NORM_FILE)
    factorization <- nmf.dgelist(dgelist, n.classes=RANK,scaling.iters=SCALING_ITERS,nmf.runs=NMF_RUNS,order.hint=order.hint )

}

nmf.report(PREFIX, factorization, glyph)

"""

@config.help("""\
Produce a Non-negative Matrix Factorization (NMF) from the output of "count:".

Non-negative Matrix Factorization is a form of fuzzy clustering. The rank r NMF of an n x m count matrix \
consists of two matricies, an n x r matrix and an r x m matrix, that multiply together to approximate the \
original count matrix.

The NMF algorithm used here is the Lee algorithm implemented in the R+ package NMF, \
which seeks to minimize the squared error. \
Using a method inspired by the "voom" function of limma, \
rows in the count matrix are weighted as a smooth function of the row average, \
such that each row contributes approximately equally to squared error. 
""")
@config.Int_flag('rank', 'Number of classes in factorization.')
@config.Int_flag('nmf_runs', 'Number of times to run NMF algorithm. (NMF algorithm does not always converge to the optimum.)')
@config.Int_flag('scaling_iters', 'Number of times to iterate to determine the row weighting.')
@config.String_flag('order_hint', '.rds file from a previous NMF. Classes will be ordered as similarly as possible to this.') 
@config.Int_flag('min_total', 'Exclude genes with less than this total number of reads')
@config.Int_flag('min_max', 'Exclude genes with no sample having at least this many reads')
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.String_flag('glyph', 'R code for a function(png.filename, sample.weights) to produce a glyph for each class.')
@config.Positional('counts', 'File containing output from "count:"')
@config.Hidden('just_report')
class NMF(config.Action_with_prefix):
    rank = 5
    nmf_runs = 10
    scaling_iters = 5
    order_hint = None
    min_total = 10
    min_max = 0
    counts = None
    norm_file = None
    glyph = None
    
    just_report = False

    def run(self):
        run_script(NMF_SCRIPT,
            PREFIX=self.prefix,
            RANK=self.rank,
            NMF_RUNS=self.nmf_runs,
            SCALING_ITERS=self.scaling_iters,
            ORDER_HINT=self.order_hint,
            MIN_TOTAL=self.min_total,
            MIN_MAX=self.min_max,
            COUNTS=self.counts,
            NORM_FILE=self.norm_file,
            JUST_REPORT=self.just_report,
            GLYPH=self.glyph
        )
            








