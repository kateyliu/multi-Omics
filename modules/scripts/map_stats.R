#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

map_stats_f <- function(map_stats_in, map_stats_out) {
    data <- read.csv(map_stats_in, sep=",", header=TRUE, check.names=F )
    #print(data)
    #rownames(data) <- data[,1]
    #data[,1] <- NULL
    x <- data.frame(data)
    colnames(x) <- c('Sample', 'Total_Reads', 'Mapped_Reads')
    
    #print(x)
    x1 <- melt(x, id.var="Sample")
    #print(x1)
    
    #png(map_stats_out, width = 8, height = 8, unit="in",res=300)
    #CHECK for crowded plots! MAX 50 samples for default
    if (nrow(x) < 50) {
        png(map_stats_out, width = 2400, height = 2400, unit="px",res=300)
    } else {
        img_h = 50 * nrow(x) #each bar gets 50px height
        png(map_stats_out, width = 2400, height = img_h, unit="px",res=300)
    }

    upper_limit <- max(x$Total_Reads)
    limits <- seq( 0, upper_limit, length.out=10)

    cust_labels <- vector("character",length=length(limits))

    if( nchar(upper_limit) < 7 ) {
         cust_labels <- paste(round(limits/1000),"K",sep="") 
         limits <- round(limits/1000) * 1000
    } else {
         cust_labels <- paste(round(limits/1000000),"M",sep="") 
         limits <- round(limits/1000000) * 1000000
    }

    #NOTE: lightblue-ish #91b6d4
    #colors <- c(Total_Reads="Grey", Mapped_Reads="#91b6d4")
    colors <- c(Total_Reads="Grey", Mapped_Reads="steelblue")
    #KEY directive to make the barplot horizontal: coord_flip()
    #ggplot(x1, aes(x=value, y=Sample, fill=variable)) + geom_bar(stat="identity", position="identity") + coord_flip()
    ggplot(x1,aes(x=Sample, y=value, fill=variable)) +
        geom_bar( stat = "identity", position="identity") +
        scale_y_continuous("",limits=c(0,upper_limit), labels=cust_labels, breaks=limits) +
        scale_fill_manual(values=colors) +
        labs( title="Read Alignment Report\n\n", x = "", y="") +
        guides(fill=guide_legend(title=NULL)) +
        theme_bw() +
        theme(axis.text.y = element_text(angle=0, hjust = 1, vjust=0.5, size=10), legend.position="top") +
        coord_flip()

    #dev.off()

}
args <- commandArgs( trailingOnly = TRUE )
arg_in = args[1]
arg_out = args[2]
#...
map_stats_f(arg_in, arg_out)
