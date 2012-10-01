
import sys, re, subprocess

from nesoni import grace, config, legion, io

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
@config.Positional('counts', 'Output from "count:"')
class Plot_counts(config.Action_with_prefix, R_action):
    counts = None

    script = r"""
    library(nesoni)
    dgelist <- read.counts(counts)    
    counts <- dgelist$counts
    n <- ncol(counts)
    
    pngname <- sprintf('%s-count.png', prefix)
    png(pngname, width=n*300, height=n*300 )
    
    not_first <- FALSE
    for(i in 1:n) {
        for(j in 1:n) {
            if (i < j) {
                par(fig=c((i-1)/(n-1),i/(n-1),(n-j)/(n-1),(n-j+1)/(n-1)), new=not_first)
                plot(counts[,i], counts[,j], 
                     log='xy', pch=19, pty='s',
                     xlab=colnames(counts)[i], ylab=colnames(counts)[j])
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

"voomed" columns in the output are log2 reads per million,
normalized using EdgeR's Trimmed Mean normalization.  
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
         "mean %.2f  sd %.2f  min %.2f  max %.2f",
         mean(vec),
         sd(vec),
         min(vec),
         max(vec)
    );
}


#data <- read.delim(FILENAME, check.names=FALSE)
#rownames(data) <- data$Feature
#
#counts <- data[,2:(N_ALL_SAMPLES+1), drop=FALSE]
#
#counts <- counts[,KEEP, drop=FALSE]
#
#good   <- (rowSums(counts) >= MIN_COUNT)
#
#data   <- data[good,, drop=FALSE]
#counts <- counts[good,, drop=FALSE]

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
        dispersion <- d$trended.dispersion
        summary(dispersion)
    } else if (MODE == 'common') {
        cat('estimateGLMCommonDisp\n')
        d <- estimateGLMCommonDisp(d, design=design)
        dispersion <- d$common.dispersion
        print(dispersion)
    } else {
        cat('Poisson dispersion\n')
        dispersion <- 1e-6    #See edgeR Poisson example
    }
    
    cat('glmFit\n')
    fit <- glmFit(d, design=design, dispersion=dispersion)
    
    cat('glmLRT\n')
    if (USE_CONTRAST) {
        contrast <- CONTRAST_WEIGHTS
        print(contrast)
        lrt <- glmLRT(d, fit, contrast=contrast)
    } else {
        lrt <- glmLRT(d, fit, coef=1:N_TO_TEST)
    }
    
    # Construct result frame
    result <- data.frame(
        Feature = rownames(dgelist$counts)
    )
    
    # 10/9/2012: logConc changed to logCPM (concentration per million?)
    result$"log2 average per million" <- lrt$table$logCPM

    if (USE_CONTRAST) {
        result$"log2 contrast" <- lrt$table$logFC
    }
    
    for(i in 1:N_TO_TEST) {
        result[, colnames(design)[i]] <- lrt$coefficients.full[,i] 
        #10/9/2012: Coefficients full now appears to be in log2
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
    
    significant <- (result$FDR <= FDR_CUTOFF)
        
    if (OUTPUT_ALL) {
        results_to_output <- result
    } else {
        results_to_output <- result[significant, ]
    }
    blank <- mapply(function(i){ all(is.na(results_to_output[,i])) }, 1:ncol(result))
    write.table(results_to_output[,!blank], OUTPUT_FILENAME, sep='\t', na='', quote=FALSE, row.names=FALSE)
    
    if (USE_CONTRAST) {
        pngname = sprintf('%s.png', OUTPUT_PLOT)
        png(pngname, width=800, height=800 )

        sane <- abs(result[,'log2 contrast']) < 1000.0

        plot(result[sane,"log2 average per million"], 
             result[sane,"log2 contrast"],
             xlab="log2 average per million",
             ylab="log2 contrast",
             main='Contrast',
             pch=19, cex=0.25) 
        points(result[significant & sane,"log2 average per million"],
               result[significant & sane,"log2 contrast"],
               pch=19, cex=0.25,
               col="red")
        dev.off()
    } else {
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
                       pch=19, cex=0.25,
                       col="red")
            
                #points((result$"log2 average per million"),
                #       sqrt( dispersion[reordering] ),
                #       pch=19,
                #       col="blue")
                
                dev.off()
            }
    }
    
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
    
    cat('\n\nWith a False Discovery Rate of', FDR_CUTOFF, '\n')
    cat('\n',sum(significant), 'significantly differentially expressed tags\n')
    
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



VOOM = COMMON + r"""

library(limma)

if (MODE == 'nullvoom') {
    voom.design <- design[, (N_TO_TEST+1):ncol(design), drop=FALSE]
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

fit <- lmFit(y, design)
if (USE_CONTRAST) {
    fit <- contrasts.fit(fit, CONTRAST_WEIGHTS)
    coef <- data.frame('log2 contrast' = fit$coefficients[,1], check.names=FALSE)
} else {
    coef <- fit$coefficients[,1:N_TO_TEST, drop=FALSE]
}

if (MODE == 'voom') {
    fit <- eBayes(fit)
    
    if (USE_CONTRAST) {
        top <- topTable(fit, coef=1, sort.by='none', number=Inf)
    } else {    
        top <- topTable(fit, coef=1:N_TO_TEST, sort.by='none', number=Inf)
    }
     
    ordering <- order(top$P.Value)
    p <- top$P.Value

} else {
    stopifnot(MODE == 'nullvoom')
    stopifnot(!USE_CONTRAST) #Not supported.

    nullfit <- lmFit(y, voom.design)
    df <- ncol(dgelist$counts) - ncol(voom.design)
    s2 <- nullfit$sigma ^ 2
    
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
    scale <- df*length(s2)/sum(s2)

    png(sprintf('%s-qq.png', OUTPUT_PLOT), width=800, height=800 )
    qqplot(qchisq( ((1:length(s2))-0.5)/length(s2) , df)/scale, s2,
        xlab='Fitted chi-square quantiles',
        ylab='Sample quantiles',
        main='QQ plot of sigma^2 unexplained by null model vs chi-squared distribution' 
    ) 
    abline(0,1)
    dev.off()
    
    ordering <- order(s2, decreasing=TRUE) 
    p <- pchisq(s2*scale, df, lower.tail=FALSE)
}

result <- data.frame(Feature = rownames(y$genes))

result$'average expression (log2 reads-per-million)' <- fit$Amean

for(i in 1:ncol(coef)) {
    result[,colnames(coef)[i]] <- coef[,i]
}

result$p <- p

result$FDR <- p.adjust(p)

#for(i in (N_ALL_SAMPLES*2+2):ncol(data)) {
#    result[,colnames(data)[i]] = data[,i]
#}
for(i in 1:ncol(y$genes)) {
    result[,colnames(y$genes)[i]] <- y$genes[,i]
}

for(i in 1:ncol(y$E)) {
    result[,sprintf('%s voomed', colnames(y$E)[i])] <- y$E[,i]
}

#for(i in 2:(N_ALL_SAMPLES*2+1)) {
#    result[,colnames(data)[i]] = data[,i]
#}

#------------------------------------------------------
result <- result[ordering,]

significant <- (result$FDR <= FDR_CUTOFF)
    
if (OUTPUT_ALL) {
    results_to_output <- result
} else {
    results_to_output <- result[significant,]
}
blank <- mapply(function(i){ all(is.na(results_to_output[,i])) }, 1:ncol(result))
write.table(results_to_output[,!blank,drop=FALSE], OUTPUT_FILENAME, sep='\t', na='', quote=FALSE, row.names=FALSE)

if (ncol(coef) == 1) # <-- Plots don't seem useful for ANOVA
    for(i in 1:ncol(coef)) {
        colname <- colnames(coef)[i]    
        
        #Strip "log2 "
        name <- substr(colname, 6,nchar(colname))
        
        if (ncol(coef) == 1) {
            pngname = sprintf('%s.png', OUTPUT_PLOT)
        } else {
            pngname = sprintf('%s-%s.png', OUTPUT_PLOT, name) 
        }
    
        png(pngname, width=800, height=800 )
        plot(result[,"average expression (log2 reads-per-million)"], 
             result[,colname],
             xlab="log2 concentration",
             ylab="log2 fold change",
             main=name,
             pch=19, cex=0.25)
             
        points(result[,"average expression (log2 reads-per-million)"][significant],
               result[,colname][significant],
               pch=19, cex=0.25,
               col="red")
        dev.off()
    }


sink(LOG_FILENAME, split=TRUE, append=TRUE)

cat('                   Sample    Library size   Further Trimmed-Mean normalization\n')
for(i in 1:ncol(dgelist$counts)) {
    cat(sprintf("%25s %15d       %f\n", rownames(dgelist$samples)[i], dgelist$samples$lib.size[i], dgelist$samples$norm.factors[i]))
}
cat('\n\n')

cat(dgelist$original.number.of.genes, 'features\n')
cat('Discarded', dgelist$original.number.of.genes-nrow(dgelist$counts), 'features with total count less than', MIN_COUNT, '\n')
cat('Kept', nrow(dgelist$counts), 'features\n')

if (MODE == 'voom') {
    cat(fit$df.prior, 'df prior (the prior is like this many extra samples)\n\n')
}

cat('With a False Discovery Rate of', FDR_CUTOFF, 'there are\n')
cat(sum(significant), 'significantly differentially expressed tags\n\n')

cat("Limma's convest function estimates that the true number of tags\n")
cat("which are differentially expressed is around", (1.0-convest(p)) * length(p), "\n\n")

"""






#cat('estimateCRDisp\n')
#d <- estimateCRDisp(d, design=design, tagwise=%(tagwise_literal)s, prior.n=%(prior_n)f, trend=%(trend_literal)s)
#
#if (%(tagwise_literal)s) {
#    dispersion <- d$CR.tagwise.dispersion
#} else {
#    dispersion <- d$CR.common.dispersion
#}

def term_specification(term):
    if '=' not in term: return term
    return term.split('=',1)[0]

def term_name(term):
    if '=' not in term: return term
    return term.split('=',1)[1]


@config.help("""\
Find significant differential expression from the output of "count:".

Terms are a regular expression on the sample names, or several 
terms joined by ^, indicating an interaction term (ie the XOR of the 
component terms). "with:" terms are included in the linear model, but 
not the significance test. The model will also include a constant term.

Examples: A A^B

A term can be given a more intuitive name in the output by appending "=name", eg

Examples: A=fold-change A^B=interaction-between-A-and-B
""" +
''.join(
        '\n\n'+
        '='*20 +
        '\n--mode %s \n' % mode +
        MODE_HELP[mode]        
        for mode in ('voom', 'nullvoom', 'poisson', 'common', 'trend')
))
@config.String_flag('mode', 'Analysis mode. See below.')
@config.Int_flag('min_count', 'Discard features with less than this total count.')
@config.Float_flag('fdr', 'False Discovery Rate cutoff for statistics and plots.')
@config.Bool_flag('output_all', 'List all genes in output, not just significant ones.')
@config.Bool_flag('constant_term', 'Include a constant term in the model.')
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.Bool_flag('tell', 'Output R code instead of executing it.')
@config.Positional('counts_file', 'The table of counts produced by "count:".')
@config.Main_section('test', 'Terms to test.')
@config.Float_section('contrast', 'Instead of doing an ANOVA on multiple terms, perform a contrast with the given weights.')
@config.Section('with_', 'Also include these terms in the model.')
@config.Section('use', 'Only use samples matching these regexes.')
class Test_counts(config.Action_with_prefix):
    min_count = 10
    constant_term = True
    mode = 'voom'
    quantile_norm = False
    fdr = 0.01
    output_all = True
    tell = False
    counts_file = None
    norm_file = None
    
    test = [ ]
    contrast = [ ]
    with_ = [ ]
    use = [ ] 

    def run(self):
        test_counts_run(    
            min_count=self.min_count, constant_term=self.constant_term, mode=self.mode, 
            quantile_norm=self.quantile_norm, fdr=self.fdr, output_all=self.output_all, 
            only_tell=self.tell,
            output_prefix=self.prefix, filename=self.counts_file, 
            test_terms=self.test, with_terms=self.with_, contrast_weights=self.contrast, use_terms=self.use,
            norm_file=self.norm_file,
        )


#def test_counts_main(args):
#    #use_expr, args = grace.get_option_value(args, '--use', str, '')
#    min_count, args = grace.get_option_value(args, '--min-count', int, 10)
#    #tagwise, args = grace.get_option_value(args, '--tagwise', grace.as_bool, False)
#    #trend, args = grace.get_option_value(args, '--trend', grace.as_bool, True)
#    #prior_n, args = grace.get_option_value(args, '--prior', float, 0.0)
#    constant_term, args = grace.get_option_value(args, '--constant-term', grace.as_bool, True)
#    
#    mode, args = grace.get_option_value(args, '--mode', str, 'voom')
#    
#    quantile_norm, args = grace.get_option_value(args, '--quantile-norm', grace.as_bool, False)
#    
#    fdr, args = grace.get_option_value(args, '--fdr', float, 0.01)
#    output_all, args = grace.get_option_value(args, '--output-all', grace.as_bool, True)
#    only_tell, args = grace.get_option_value(args, '--tell', grace.as_bool, False)
#    
#    if len(args) < 2:
#        print >> sys.stderr, TEST_COUNTS_HELP
#        
#        for mode in ('voom', 'nullvoom', 'poisson', 'common', 'trend'):
#            print >> sys.stderr, '='*20
#            print >> sys.stderr, '--mode', mode
#            print >> sys.stderr, MODE_HELP[mode]
#        
#        raise grace.Help_shown()
#    
#    output_prefix, filename, args = args[0], args[1], args[2:]
#    
#    test_terms = [ ]
#    with_terms = [ ]
#    contrast_weights = [ ]
#    use_terms = [ ]
#    def default(args):
#        grace.expect_no_further_options(args)
#        test_terms.extend(args)
#    def with_(args):    
#        grace.expect_no_further_options(args)
#        with_terms.extend(args)
#    def contrast(args):
#        # -1 is not a flag... grace.expect_no_further_options(args)
#        contrast_weights.extend([ float(item) for item in args ])
#    def use(args):
#        grace.expect_no_further_options(args)
#        use_terms.extend(args)
#    
#    grace.execute(args, {'with': with_, 'contrast':contrast, 'use':use}, default)
#
#    test_counts_run(    
#        min_count, constant_term, mode, quantile_norm, fdr, output_all, only_tell,
#        output_prefix, filename, test_terms, with_terms, contrast_weights, use_terms,
#        None
#    )

def test_counts_run(
    min_count, constant_term, mode, quantile_norm, fdr, output_all, only_tell,
    output_prefix, filename, test_terms, with_terms, contrast_weights, use_terms,
    norm_file
):
    log = grace.Log()

    assert mode in ('voom', 'nullvoom', 'poisson', 'common', 'trend')
    
    if not use_terms:
        use_terms = [''] #Default to all
    
    assert test_terms, 'No terms to test'
    
    use_contrast = bool(contrast_weights)
    assert not use_contrast or len(contrast_weights) == len(test_terms), 'Need one contrast weight per term'
    
    all_terms = test_terms + with_terms 
    if constant_term: all_terms += [ '' ]
    n_to_test = len(test_terms)
    
    if use_contrast:
        contrast_weights += [ 0.0 ] * (len(all_terms)-len(test_terms))
    
    #f = open(filename,'rb')
    #header = f.readline()
    #f.close()
    #parts = header.rstrip('\n').split('\t')
    #
    #n_all_samples = parts.index('RPKM '+parts[1])-1    
    #all_samples = parts[1:n_all_samples+1]
    
    data = io.read_grouped_table(filename,{'Count':str,'Annotation':str})
    assert len(data['Count']), 'Count file is empty'
    
    all_samples = data['Count'].values()[0].keys()
    n_all_samples = len(all_samples)
    
    
    keep = [ any( re.search(use_expr, item) for use_expr in use_terms ) for item in all_samples ]
    samples = [ all_samples[i] for i in xrange(n_all_samples) if keep[i] ]
    n_samples = len(samples)
    
    model = [ ]
    for term in all_terms:
        column = [ False ] * n_samples
        for expression in term_specification(term).split('^'):
            for i in xrange(n_samples):
                if re.search(expression, samples[i]):
                    column[i] = not column[i]
        model.append([ int(item) for item in column ])
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
    if mode not in ('poisson','nullvoom'):
        assert len(all_terms) < n_all_samples, (
            'Can\'t have as many linear model terms as samples, within group variation can\'t be estimated. '
            'I wouldn\'t recommend using --mode poisson or --mode nullvoom, '
            'you should go and sequence some more samples.'
        )
    
    
    log.log('Analysis mode: %s\n%s\n'%(mode,MODE_HELP[mode]))
    
    log_filename = output_prefix + '-info.txt'
    design_columns = [ 'log2 '+term_name(item) for item in all_terms ]
    term_names = [ term_name(item) for item in all_terms ]

    #filename_literal = R_literal(filename)
    #mode_literal = R_literal(mode)
    #quantile_norm_literal = R_literal(quantile_norm)
    #keep_literal = R_literal(keep)
    #output_filename_literal = R_literal(output_prefix + '.txt')
    #log_filename_literal = R_literal(log_filename)
    #output_plot_literal = R_literal(output_prefix)
    #
    #use_contrast_literal = R_literal(use_contrast)
    #contrast_weights_literal = R_literal(contrast_weights)
    #
    ##tagwise_literal = R_literal(tagwise)
    ##trend_literal = R_literal(trend)
    #
    #design_literal = 'matrix(%s, nrow=%d, byrow=TRUE)' % (
    #    R_literal( sum(model, ()) ),
    #    n_samples
    #)
    #
    #design_columns_literal = R_literal(design_columns)
    #terms_literal = R_literal(term_names)
    
    log.attach(open(log_filename,'wb'))
    log.close()    
    
    if mode in ('voom', 'nullvoom'):
        script = VOOM
    else:
        script = EDGER
    
    run_script(
        script,
        
        FILENAME = filename,
        N_SAMPLES = n_samples,
        N_ALL_SAMPLES = n_all_samples,
        N_TO_TEST = n_to_test,
        FDR_CUTOFF = fdr,
        OUTPUT_ALL = output_all,
        MIN_COUNT = min_count,
        MODE = mode,
        QUANTILE_NORM = quantile_norm,
        KEEP = keep,
        OUTPUT_FILENAME = output_prefix + '.txt',
        LOG_FILENAME = log_filename,
        OUTPUT_PLOT = output_prefix,
        USE_CONTRAST = use_contrast,
        CONTRAST_WEIGHTS = contrast_weights,
        DESIGN = model,
        DESIGN_COLUMNS = design_columns,
        TERM_NAMES = term_names,
        NORM_FILE = norm_file,
        
        only_tell=only_tell
    )
        
    return fdr #test-power needs this


TEST_POWER_HELP = """\

Usage:

    nesoni test-power: [options] output_prefix [of: [test-count options]] 

Test the statistical power of "nesoni test-counts".

Options:

    --m NNN            - Number of differentially expressed tags
    
    --n NNN            - Total number of tags
    
    --reps N           - Within group replicates
    
    --count NNNN       - Average count
    
    --dispersion N.NN  - Dispersion
    
    --log-fold N.NN    - log2 Fold change of differentially expressed tags

"""

POWER_TEMPLATE = r"""

library(MASS)

m <- %(m)d
n <- %(n)d
reps <- %(reps)d
avg_count <- %(count)d
dispersion <- %(dispersion)f
log_fold <- %(log_fold)f

mynegbin <- function(u, d) {
    rnegbin(1, u, 1.0/d)
}

counts <- matrix(0, nrow=n, ncol=(reps*2) )
names <- rep('', n)

for(i in 1:n) {
    a <- 1.0
    b <- 1.0
    if (i <= m) {
        b <- b * (2 ** sample(c(log_fold,-log_fold),1))
        names[i] <- paste('DifferentialTag',i,sep='')
    } else {
        s <- 0.0
        names[i] <- paste('Tag',i,sep='')
    }
    scale <- avg_count / (0.5*(a+b))
    a <- a * scale
    b <- b * scale
    
    for(j in 1:reps) {
        counts[i,j]      <- mynegbin(a, dispersion)
        counts[i,reps+j] <- mynegbin(b, dispersion)
    }
}

data <- data.frame(
    Feature=names,
    counts,
    counts,
    Product=names
)

for(i in 1:reps) {
   colnames(data)[i+1] <- paste('Control',i,sep='')
   colnames(data)[reps+i+1] <- paste('Experimental',i,sep='')
}
for(i in 1:(reps*2)) {
   colnames(data)[i+reps*2+1] <- paste('RPKM', colnames(data)[i+1])
}

write.table(data, %(filename_literal)s, sep='\t', quote=FALSE, row.names=FALSE)

"""

POWER_REPORT_TEMPLATE = """

claimed_fdr <- %(claimed_fdr)f

data <- read.delim(%(output_filename_literal)s, check.names=FALSE)

false_positives <- 0
true_positives <- 0

for(i in 1:nrow(data)) {
    sig <- data$FDR[i] <= claimed_fdr
    true <- regexpr('Differential', data$Feature[i])[[1]] != -1
    
    if (sig && !true) { false_positives <- false_positives + 1 }
    if (sig && true) { true_positives <- true_positives + 1 }
}

sink(%(log_filename_literal)s, append=TRUE, split=TRUE)

cat('\n\nFalse positives:', false_positives, '\n')
cat(    'False negatives:', %(m)d - true_positives, '\n')
cat(    'True positives:', true_positives, '\n')

"""

def test_power_main(args):
    m, args = grace.get_option_value(args, '--m', int, 10)
    n, args = grace.get_option_value(args, '--n', int, 1000)
    reps, args = grace.get_option_value(args, '--reps', int, 2)
    count, args = grace.get_option_value(args, '--count', int, 100)
    dispersion, args = grace.get_option_value(args, '--dispersion', float, 0.1)
    log_fold, args = grace.get_option_value(args, '--log-fold', float, 1.0)

    if len(args) < 1:
        print >> sys.stderr, TEST_POWER_HELP
        raise grace.Help_shown()
    
    output_prefix, args = args[0], args[1:]

    options = [ ]
    def of(args):
        options.extend(args)    
    grace.execute(args, {'of': of})
    
    filename = output_prefix + '-input.txt'
    filename_literal = R_literal(filename)
    log_filename_literal = R_literal(output_prefix + '-info.txt')
    
    run_script(POWER_TEMPLATE % locals())

    claimed_fdr = test_counts_main([ output_prefix, filename, 'Experimental' ] + options)

    output_filename_literal = R_literal(output_prefix + '.txt')

    run_script(POWER_REPORT_TEMPLATE % locals())





HEATMAP_SCRIPT = """

library(nesoni)

dgelist <- read.counts(COUNTS, min.total=MIN_TOTAL, min.max=MIN_MAX, norm.file=NORM_FILE)

if (length(ORDER)) {
    dgelist <- dgelist[, ORDER]
}

elist <- voom(dgelist, normalize.method=if(QUANTILE) 'quantile' else 'none')
hmap.elist(PREFIX, elist, min.sd=MIN_SD, min.span=MIN_SPAN, min.svd=MIN_SVD, svd.rank=SVD_RANK, reorder.columns=REORDER_COLUMNS)

"""


@config.help("""\
Produce a heatmap (and table) of voomed expression levels (log2(count+0.5) normalized by effective library size).

Hierachical clustering and ordering of rows is performed using the "seriation" package in R+.

You will typically need to turn on both an "expression level filter" and a "fold change filter" \
in order to show a useful list of genes.
""")
@config.Int_flag('min_total', 'Expression level filter:\nExclude genes with less than this total number of reads')
@config.Int_flag('min_max', 'Expression level filter:\nExclude genes with no sample having at least this many reads')
@config.Float_flag('min_sd', 'Fold change filter:\nExclude genes with less than this standard deviation of log2 expression levels')
@config.Float_flag('min_span', 'Fold change filter:\nExclude genes where there is no pair of samples that differ in log2 expression by this much')
@config.Float_flag('min_svd', 'Fold change filter:\nUsing the SVD of log2 expression levels, '
                              'exclude genes where the sum of squares of the U matrix row for the gene is less than this many standard deviations.')
@config.Int_flag('svd_rank', 'Only use the top this many columns of the SVD U matrix when applying --min-svd.')
@config.Bool_flag('reorder_columns', 'Cluster and optimally order columns as well as rows.')
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.Bool_flag('quantile', 'Use quantile normalization.')
@config.Positional('counts', 'File containing output from "nesoni count:"')
@config.Section('order', 'Optionally, specify an order to show the columns in.')
class Heatmap(config.Action_with_prefix):
    counts = None
    min_total = 10
    min_max = 0
    min_sd = 0.0
    min_span = 0.0
    min_svd = 0.0
    svd_rank = None
    reorder_columns = False
    norm_file = None
    quantile = False
    order = [ ]

    def run(self):
        run_script(HEATMAP_SCRIPT,
            PREFIX=self.prefix,
            COUNTS=self.counts,
            MIN_TOTAL=self.min_total,
            MIN_MAX=self.min_max,
            MIN_SD=self.min_sd,
            MIN_SPAN=self.min_span,
            MIN_SVD=self.min_svd,
            SVD_RANK=self.svd_rank,
            REORDER_COLUMNS=self.reorder_columns,
            NORM_FILE=self.norm_file,
            QUANTILE=self.quantile,
            ORDER=self.order,
        )


COMPARE_TESTS_SCRIPT = """

library(nesoni)

compare.tests(
    PREFIX,
    read.delim(A),
    read.delim(B),
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

dgelist <- read.counts(sprintf("%s.txt",PREFIX), use.tmm=USE_TMM)

result <- data.frame(
    Sample = rownames(dgelist$samples),
    Normalizing.multiplier = dgelist$samples$normalizing.multiplier
)

write.table(result, sprintf("%s-norm.txt",PREFIX), sep='\t', quote=FALSE, row.names=FALSE)

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
From <prefix>.txt creates <prefix>-norm.txt containing a table of \
normalizing multipliers as calculated by EdgeR using Trimmed Mean of M-values normalization.
""")
@config.Bool_flag('tmm', 'Use TMM normalization (in addition to normalizing by number of mapped reads).')
class Norm_from_counts(config.Action_with_prefix):
    tmm = True

    def log_filename(self):
        return self.prefix + '-norm_log.txt'

    def run(self):
        run_script(NORMALIZATION_SCRIPT,
            PREFIX=self.prefix,
            USE_TMM=self.tmm,
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
            








