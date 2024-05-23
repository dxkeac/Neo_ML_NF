#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)
args = commandArgs(TRUE)
library(sequenza)

seq_file=args[1]
sample_id=args[2]
out_dir=args[3]

#getwd()
#setwd('xxx')
#help(sequenza.extract)

ext <- sequenza.extract(seq_file, verbose = T, assembly = "hg38", chromosome.list=c(paste0("chr", c(1:22, "X", "Y"))))
CP <- sequenza.fit(ext)
sequenza.results(sequenza.extract = ext, cp.table = CP, sample.id = sample_id, out.dir = out_dir)

cint <- get.ci(CP)
cellularity <- cint$max.cellularity
#write.table(cellularity, paste(out_dir,paste(sample_id,"_purity.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
ploidy <- cint$max.ploidy
purity_ploidy <- paste(cellularity,ploidy,sep="\t")

write.table(purity_ploidy, paste(out_dir,paste(sample_id,"_purity_ploidy.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)

#write.table(ploidy, paste(out_dir,paste(sample_id,"_ploidy.txt",sep=""),sep='/'), col.names = FALSE, row.names = FALSE, sep = "\t", quote=FALSE)
