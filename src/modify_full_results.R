#!/usr/bin/env Rscript

## This script is a modified function from the script "sourcetracker_for_qiime.r" by Dan Knights.
## The function that creates the folder full_results is modified to fit our desires.
## Created by Moa Haegglund


helpstr <- c('-i OTU table, QIIME formatted in classic format.\n-t Utilized reference taxonomy (e.g. taxonomy file from GreenGenes).\n-r The results file of the SourceTracker analysis (results.RData)\n-o Name of the folder where the modified full results data will be stored.')

allowed.args <- list('-i' = NULL,'-t' = NULL, '-r' = NULL, '-o' = NULL)

"parse.args" <- function(allowed.args,helplist=NULL){
    argv <- commandArgs(trailingOnly=TRUE)
    # print help string if requested
    if(!is.null(helpstr) && sum(argv == '-h')>0){
        cat('',helpstr,'',sep='\n')
        q(runLast=FALSE)
    }
    argpos <- 4
    for(name in names(allowed.args)){
        argpos <- which(argv == name)
        if(length(argpos) > 0){
            # test for flag without argument
            if(argpos == length(argv) || substring(argv[argpos + 1],1,1) == '-')
                allowed.args[[name]] <- TRUE
            else {
                allowed.args[[name]] <- argv[argpos + 1]
            }
        }
    }
    return(allowed.args)
}

# parse arg list
arglist <- parse.args(allowed.args)

otu_table <- read.table(arglist[['-i']],sep='\t',skip=1,header=T,row.names=1,comment='')
taxonomy <- arglist[['-t']]
results <- arglist[['-r']]
output_name <- arglist[['-o']]


load(results)
taxa <- read.table(taxonomy,sep='\t',row.names=1,comment='')
colnames(taxa)<-'taxonomy'
res.mean <- apply(results$full.results,c(2,3,4),mean)
sample.sums <- apply(results$full.results[1,,,,drop=F],c(3,4),sum)

# create dir
subdir <- paste(output_name,sep='/')
dir.create(subdir,showWarnings=FALSE, recursive=TRUE)
# write each env separate file
for(i in 1:length(results$train.envs)){
    env.name <- results$train.envs[i]
    filename <- sprintf('%s/%s_contributions.txt', subdir, env.name)
    res.mean.i <- res.mean[i,,]
# handle the case where there is only one sink sample
    if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)
    env.mat <- res.mean.i/sample.sums
    colnames(env.mat)<-results$samplenames
    env.mat[is.na(env.mat)]<-0
    row.number= 1:dim(env.mat)[1]
    env.mat2 <- cbind(env.mat,row.number)
    rownames(env.mat2) <- rownames(otu_table)
    with.taxa <- merge(env.mat2,taxa,by='row.names')
    rownames(with.taxa)<-with.taxa[,1]
    with.taxa <- with.taxa[,-1]
    if (length(results$samplenames)==1)  {
        with.taxa2 <- with.taxa[with.taxa[,1]!=0,]
    } else {
        with.taxa2 <- with.taxa[rowSums(with.taxa[,1:length(results$samplenames)])!=0,]
    }
    if (dim(with.taxa2)[1]!=0) {
        sink(filename)
        cat('OTU.ID\t')
        write.table(with.taxa2,quote=F,sep='\t')
        sink(NULL)
    }
}
