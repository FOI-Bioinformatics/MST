#!/usr/bin/env Rscript

## This script is a modified function from the script "sourcetracker_for_qiime.r" by Dan Knights.
## The function that creates the folder full_results is modified to fit our desires.


helpstr <- c('-d Path to directory where the modified full result files are stored..\n-s The source to evaluate.\n--sample_ids Comma-separated list of sampleIDs to evaluate [optional].\n-p Path to store files\n-v Value for which probabilities higher than should be colored red (e. g. 0.7) [optional].\n-m Include OTUs with probability larger than.')

allowed.args <- list('-d' = NULL, '-s' = NULL, '--sample_ids' = NULL, '-p' = NULL, '-v' = NULL, '-m' = NULL)

"parse.args" <- function(allowed.args,helplist=NULL){
    argv <- commandArgs(trailingOnly=TRUE)
    # print help string if requested
    if(!is.null(helpstr) && sum(argv == '-h')>0){
        cat('',helpstr,'',sep='\n')
        q(runLast=FALSE)
    }
    argpos <-6
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

arglist <- parse.args(allowed.args)
source <- arglist[['-s']]
folder <- arglist[['-p']]
min_value <- as.numeric(arglist[['-m']])

file <- sprintf('%s%s_contributions.txt',arglist[['-d']],source)
datat = read.table(file,sep='\t',header=T, comment='')
datat2 = datat[,1:(dim(datat)[2]-2)]

if (is.null(arglist[['-v']])==FALSE){
value <- as.numeric(arglist[['-v']])}

if (is.null(arglist[['--sample_ids']])==FALSE){
samples <- as.list(strsplit(arglist[['--sample_ids']],",")[[1]])
} else {
all_columns <- as.list(colnames(datat2))
samples <- all_columns[-1]
}

for (i in samples){ls 
	tmp <- cbind(datat2[,1], datat2[,i])
	output_name <- sprintf('%s%s_%s.txt', folder, source, i)
	write.table('ring_internal_separator_thickness\t1\t0.5\nring_width\t1\t0.5', output_name, row.names=F, col.names=F, quote=F)
	output <- as.data.frame(tmp[which(tmp[,2]>=min_value),])
	output$height <- "ring_height"
	output$value <- 1
	output <- output[,c('V1','height', 'value','V2')]
	head(output)
	if (is.null(arglist[['-v']])==FALSE){	
		output2 <- as.data.frame(tmp[which(tmp[,2]>=min_value),])
		output2$height <- "ring_color"
		output2$value <- 1
		output2[which(output2$V2>=value),"V2"]<-'#e71b05'
		output2[which(output2$V2!='#e71b05'),"V2"]<-'#000000'	
		output = rbind(output, output2)}
	output <- output[,c('V1','height', 'value','V2')]
	write.table(output, output_name, sep='\t',row.names=F,col.names=F, append = TRUE, quote=F)
}
