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
x = read.delim(snakemake@input[['all']],header=T)



##################################################
#Sum of normalized reads per developmental stages#
##################################################
x$sumPRE_AABB <- rowSums(x[ , c(76:78)], na.rm=TRUE)
x$sumMID_AABB <- rowSums(x[ , c(79:81)], na.rm=TRUE)
x$sumPOS_AABB <- rowSums(x[ , c(82:84)], na.rm=TRUE)

x$sumPRE_aabb <- rowSums(x[ , c(85:87)], na.rm=TRUE)
x$sumMID_aabb <- rowSums(x[ , c(88:90)], na.rm=TRUE)
x$sumPOS_aabb <- rowSums(x[ , c(91:93)], na.rm=TRUE)

x$sumPRE_aAbb <- rowSums(x[ , c(94:96)], na.rm=TRUE)
x$sumMID_aAbb <- rowSums(x[ , c(97:99)], na.rm=TRUE)
x$sumPOS_aAbb <- rowSums(x[ , c(100:102)], na.rm=TRUE)

x$sumPRE_aabB <- rowSums(x[ , c(103:105)], na.rm=TRUE)
x$sumMID_aabB <- rowSums(x[ , c(106:108)], na.rm=TRUE)
x$sumPOS_aabB <- rowSums(x[ , c(109:111)], na.rm=TRUE)

x$sumPRE_aabb_18 <- rowSums(x[ , c(112:114)], na.rm=TRUE)
x$sumMID_aabb_18 <- rowSums(x[ , c(115:117)], na.rm=TRUE)
x$sumPOS_aabb_18 <- rowSums(x[ , c(118:120)], na.rm=TRUE)

x$sumPRE_aabb_22 <- rowSums(x[ , c(121:123)], na.rm=TRUE)
x$sumMID_aabb_22 <- rowSums(x[ , c(124:126)], na.rm=TRUE)
x$sumPOS_aabb_22 <- rowSums(x[ , c(127:129)], na.rm=TRUE)




######################################################
#Average of normalized reads per developmental stages#
######################################################
x$avgPRE_AABB <- rowMeans(x[ , c(76:78)], na.rm=TRUE)
x$avgMID_AABB <- rowMeans(x[ , c(79:81)], na.rm=TRUE)
x$avgPOS_AABB <- rowMeans(x[ , c(82:84)], na.rm=TRUE)

x$avgPRE_aabb <- rowMeans(x[ , c(85:87)], na.rm=TRUE)
x$avgMID_aabb <- rowMeans(x[ , c(88:90)], na.rm=TRUE)
x$avgPOS_aabb <- rowMeans(x[ , c(91:93)], na.rm=TRUE)

x$avgPRE_aAbb <- rowMeans(x[ , c(94:96)], na.rm=TRUE)
x$avgMID_aAbb <- rowMeans(x[ , c(97:99)], na.rm=TRUE)
x$avgPOS_aAbb <- rowMeans(x[ , c(100:102)], na.rm=TRUE)

x$avgPRE_aabB <- rowMeans(x[ , c(103:105)], na.rm=TRUE)
x$avgMID_aabB <- rowMeans(x[ , c(106:108)], na.rm=TRUE)
x$avgPOS_aabB <- rowMeans(x[ , c(109:111)], na.rm=TRUE)

x$avgPRE_aabb_18 <- rowMeans(x[ , c(112:114)], na.rm=TRUE)
x$avgMID_aabb_18 <- rowMeans(x[ , c(115:117)], na.rm=TRUE)
x$avgPOS_aabb_18 <- rowMeans(x[ , c(118:120)], na.rm=TRUE)

x$avgPRE_aabb_22 <- rowMeans(x[ , c(121:123)], na.rm=TRUE)
x$avgMID_aabb_22 <- rowMeans(x[ , c(124:126)], na.rm=TRUE)
x$avgPOS_aabb_22 <- rowMeans(x[ , c(127:129)], na.rm=TRUE)





###########################################################
#Classify sRNA loci per category of pre-/mid-/post-meiotic#
###########################################################
x$peakPRE_AABB <- ifelse(x$avgPRE_AABB>x$avgMID_AABB & x$avgPRE_AABB>x$avgPOS_AABB, "TRUE", "FALSE")
x$peakMID_AABB <- ifelse(x$avgPRE_AABB<x$avgMID_AABB & x$avgMID_AABB>x$avgPOS_AABB, "TRUE", "FALSE")
x$peakPOS_AABB <- ifelse(x$avgPRE_AABB<x$avgPOS_AABB & x$avgMID_AABB<x$avgPOS_AABB, "TRUE", "FALSE")

x$peakPRE_aabb <- ifelse(x$avgPRE_aabb>x$avgMID_aabb & x$avgPRE_aabb>x$avgPOS_aabb, "TRUE", "FALSE")
x$peakMID_aabb <- ifelse(x$avgPRE_aabb<x$avgMID_aabb & x$avgMID_aabb>x$avgPOS_aabb, "TRUE", "FALSE")
x$peakPOS_aabb <- ifelse(x$avgPRE_aabb<x$avgPOS_aabb & x$avgMID_aabb<x$avgPOS_aabb, "TRUE", "FALSE")

x$peakPRE_aAbb <- ifelse(x$avgPRE_aAbb>x$avgMID_aAbb & x$avgPRE_aAbb>x$avgPOS_aAbb, "TRUE", "FALSE")
x$peakMID_aAbb <- ifelse(x$avgPRE_aAbb<x$avgMID_aAbb & x$avgMID_aAbb>x$avgPOS_aAbb, "TRUE", "FALSE")
x$peakPOS_aAbb <- ifelse(x$avgPRE_aAbb<x$avgPOS_aAbb & x$avgMID_aAbb<x$avgPOS_aAbb, "TRUE", "FALSE")

x$peakPRE_aabB <- ifelse(x$avgPRE_aabB>x$avgMID_aabB & x$avgPRE_aabB>x$avgPOS_aabB, "TRUE", "FALSE")
x$peakMID_aabB <- ifelse(x$avgPRE_aabB<x$avgMID_aabB & x$avgMID_aabB>x$avgPOS_aabB, "TRUE", "FALSE")
x$peakPOS_aabB <- ifelse(x$avgPRE_aabB<x$avgPOS_aabB & x$avgMID_aabB<x$avgPOS_aabB, "TRUE", "FALSE")

x$peakPRE_aabb_18 <- ifelse(x$avgPRE_aabb_18>x$avgMID_aabb_18 & x$avgPRE_aabb_18>x$avgPOS_aabb_18, "TRUE", "FALSE")
x$peakMID_aabb_18 <- ifelse(x$avgPRE_aabb_18<x$avgMID_aabb_18 & x$avgMID_aabb_18>x$avgPOS_aabb_18, "TRUE", "FALSE")
x$peakPOS_aabb_18 <- ifelse(x$avgPRE_aabb_18<x$avgPOS_aabb_18 & x$avgMID_aabb_18<x$avgPOS_aabb_18, "TRUE", "FALSE")

x$peakPRE_aabb_22 <- ifelse(x$avgPRE_aabb_22>x$avgMID_aabb_22 & x$avgPRE_aabb_22>x$avgPOS_aabb_22, "TRUE", "FALSE")
x$peakMID_aabb_22 <- ifelse(x$avgPRE_aabb_22<x$avgMID_aabb_22 & x$avgMID_aabb_22>x$avgPOS_aabb_22, "TRUE", "FALSE")
x$peakPOS_aabb_22 <- ifelse(x$avgPRE_aabb_22<x$avgPOS_aabb_22 & x$avgMID_aabb_22<x$avgPOS_aabb_22, "TRUE", "FALSE")




#REMOVE LOW ABUNDANT LOCI
u=subset(x, Reads>=500)

#SORT LOCI PER CATEGORY OF PRE-/MID-/POST-MEIOTIC GROUPS
z = u %>% arrange(desc(peakPRE_AABB), peakMID_AABB, peakPOS_AABB)


#CLASSIFY 24-nt LOCI BASED ON THE PHASE SCORE; SIGNIFICANT (PhaseScore >= 40)  OR AMBIGIOUS (PhaseScore >= 30 BUT < 40) OR NOT (PhaseScore < 30)
v=subset(z, PhaseScore>=40)
w=subset(z, PhaseScore<40)
wA=subset(w, PhaseScore>=30)
wB=subset(w, PhaseScore<30)


#ADD A BIOTYPE COLUMN TO ANNOTATE  "PHASIRNA", "HC-SIRNA" AND "AMBIGUOUS" LOCI
v$biotype=c("phasiRNA")
wA$biotype=c("ambiguous")


#PARSE LOCI IN CLASS OF 21-nt and 24-nt
#PHASIRNA, SIRNA AND AMBIGUOUS CATEGORIES
#Parse phasiRNA loci in classes of 21-nt and 24-nt
v21=subset(v, DicerCall==21)
v24=subset(v, DicerCall==24)
wA21=subset(wA, DicerCall==21)
wA24=subset(wA, DicerCall==24)
wB24=subset(wB, DicerCall==24)


#Add biotype column to annotate phasiRNA
v21$category=c("21PHAS")
v24$category=c("24PHAS")
wA21$category=c("21PHAS")
wA24$category=c("24PHAS")

wB24$biotype=c("hc-siRNA")
wB24$category=c("24nt hc-siRNA")


# Write tables
write.table(z, snakemake@output[['srna']], sep = "\t", row.names=F, quote=FALSE)
write.table(v, snakemake@output[['phas']], sep = "\t", row.names=F, quote=FALSE)
write.table(wA, snakemake@output[['ambiguous']], sep = "\t", row.names=F, quote=FALSE)

write.table(v21, snakemake@output[['phas21']], sep = "\t", row.names=F, quote=FALSE)
write.table(v24, snakemake@output[['phas24']], sep = "\t", row.names=F, quote=FALSE)
write.table(wA21, snakemake@output[['ambiguous21']], sep = "\t", row.names=F, quote=FALSE)
write.table(wA24, snakemake@output[['ambiguous24']], sep = "\t", row.names=F, quote=FALSE)
write.table(wB24, snakemake@output[['sirna']], sep = "\t", row.names=F, quote=FALSE)








#EXTRACT LINES OF PHAS LOCI NOT SIGNIFICANT AT PhaseScore >= 30
y= subset(x, PhaseScore <30)

#ADD A BIOTYPE COLUMN AND ANNOTATE "NOT PHASIRNA"
y$biotype=c("Not phasirna")

#CLASSIFY NOT PHASIRNA LOCI PER CATEGORY OF PRE-/MID-/POST-MEIOTIC GROUPS
z = y %>% arrange(desc(peakPRE_AABB), peakMID_AABB, peakPOS_AABB)

#PARSE LOCI IN CLASS OF SIRNA AND AMBIGUOUS CATEGORIES
zA=subset(z, PhaseScore <20)
zB=subset(z, PhaseScore>=20)

#ADD A CATEGORY COLUMN TO ANNOTATE SIRNA AND AMBIGUOUS GROUPS
zA$category=c("sirna")
zB$category=c("ambiguous")


#WRITE TABLE
write.table(zA, snakemake@output[['sirna']], sep = "\t", row.names=F, quote=FALSE)
write.table(zB, snakemake@output[['ambiguous']], sep = "\t", row.names=F, quote=FALSE)










#########################################
# Source of information for R functions #
#########################################
#Find the sum of rows
    #data$new <- rowSums(data, na.rm=TRUE)
    #data$row_sum <- rowSums(data[ , c(1,3)], na.rm=TRUE)
    #https://www.statology.org/sum-specific-columns-in-r/
#Find the average of rows:
    #data$row_mean <- rowMeans(data, na.rm=TRUE)
    #data$new <- rowMeans(data[ , c(1,2)], na.rm=TRUE)
    #https://www.statology.org/average-across-columns-in-r/
#Logical function:
    #df$rating <- ifelse(df$points>15 | df$assists>8, "good", "bad")
    #https://www.statology.org/r-if-statement-multiple-conditions/
#Sort the dataframe
    #res = mtcars %>% arrange(mpg, cyl)
    #https://koalatea.io/r-dplyr-arrange/
