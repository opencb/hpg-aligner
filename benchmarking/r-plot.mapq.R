#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args) != 4) {
    print("Mismatch number of arguments.")
    print("Usage: r-plot.mapq.R bwa-mem-eval hpg-eval out-filename title-name")
    q()
}

bwamem_eval <- args[1]
hpg_eval <- args[2]
out_name <- args[3]
title <- args[4]

bwamem <- read.table(bwamem_eval);
hpg <- read.table(hpg_eval);

# 2,000,000 reads / 100 (to get percent ) => 20,000
bwamem$x <- bwamem$V3/bwamem$V2
bwamem$y <- bwamem$V2/20000

hpg$x <- hpg$V3/hpg$V2
hpg$y <- hpg$V2/20000

pdf(out_name)

plot(bwamem$x,bwamem$y,type="b",xlim=c(1e-7,1e-1),ylim=c(75,100), log="x",col="red",pch=".", xlab="#wrong mappings / #mapped", ylab="#mapped / total (%)",main=title)
lines(hpg$x, hpg$y,type="b",col="green", pch=".")

points(bwamem$x[41],bwamem$y[41],col="black",pch=19)
points(hpg$x[41],hpg$y[41],col="black",pch=19)

legend("bottomright",c("BWA MEM","HPG Aligner"),pch=c(".","."),lty=c(1,1), col=c("red","green"))

dev.off()

