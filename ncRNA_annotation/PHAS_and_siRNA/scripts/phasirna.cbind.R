################################
# Merge result and count files #
################################
#Load package on R3.6.1
library("edgeR")
library("dplyr")
library("tidyr")
library("data.table")


#############################
####  EXPLORE ABUNDANCE  ####
#############################
x1 = read.delim(snakemake@input[['result']],header=T)
x2 = read.delim(snakemake@input[['raw']],header=T)
x3 = read.delim(snakemake@input[['rpm']],header=T)

#Bind columns result and abundance files
#col=cbind(x1,x2,x3)
col1 <- x1 %>% full_join(x2, by = "Name")
col2 <- col1 %>% full_join(x3, by = "Name")

#Remove loci that are not Dicer products
srna=subset(col2, DicerCall != "N")

# Write tables
write.table(srna, snakemake@output[['srna']], sep = "\t", row.names=F, quote=FALSE)
