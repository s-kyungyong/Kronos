##################
# Subset results #
##################
#Load package on R3.6.1
library("edgeR")
library("dplyr")
library("tidyr")
library("data.table")


#############################
####  EXPLORE ABUNDANCE  ####
#############################
x21 = read.delim(snakemake@input[['phas21']],header=T)
x24 = read.delim(snakemake@input[['phas24']],header=T)

#Parse phasiRNA loci in pre-/mid-/post-meiotic groups
x21PRE=subset(x21, peakPRE_AABB=="TRUE")
x21NOT=subset(x21, peakPRE_AABB=="FALSE")
x24PRE=subset(x24, peakPRE_AABB=="TRUE")
x24NOT=subset(x24, peakPRE_AABB=="FALSE")


# Write tables
write.table(x21PRE, snakemake@output[['phas21PRE']], sep = "\t", row.names=F, quote=FALSE)
write.table(x21NOT, snakemake@output[['phas21NOT']], sep = "\t", row.names=F, quote=FALSE)

write.table(x21PRE$Name, snakemake@output[['list21PRE']], sep = "\t", row.names=F, quote=FALSE, col.names=F)
write.table(x21NOT$Name, snakemake@output[['list21NOT']], sep = "\t", row.names=F, quote=FALSE, col.names=F)

write.table(x24PRE, snakemake@output[['phas24PRE']], sep = "\t", row.names=F, quote=FALSE)
write.table(x24NOT, snakemake@output[['phas24NOT']], sep = "\t", row.names=F, quote=FALSE)

write.table(x24PRE$Name, snakemake@output[['list24PRE']], sep = "\t", row.names=F, quote=FALSE, col.names=F)
write.table(x24NOT$Name, snakemake@output[['list24NOT']], sep = "\t", row.names=F, quote=FALSE, col.names=F)
