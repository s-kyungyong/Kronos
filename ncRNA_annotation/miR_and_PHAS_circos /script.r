R
library(OmicCircos)

X <- read.delim("chrom_lengths.tab",sep = "\t",header=F,stringsAsFactors=FALSE)

miR <- read.delim("miRNAs.mapping.tab",sep="\t",header=F,stringsAsFactors=FALSE)

colnames(X) <- c("seg.name", "seg.Start", "seg.End", "name", "color")

PHAS21lab <- read.delim("21PHAS.mapping.tab",sep="\t",header=F,stringsAsFactors=FALSE)

PHAS24 <- read.delim("24PHAS.mapping.tab",sep="\t",header=F,stringsAsFactors=FALSE)

db <- segAnglePo(seg.dat=X, seg= X$seg.name) ;

colors <- as.vector(X$color)

mapX <- as.data.frame(cbind(X$seg.name,as.integer(X$seg.End),as.character(X$name)))

moreGray<-"#525252"


db <- as.data.frame(db, stringsAsFactors = FALSE)

# Convert relevant columns to numeric
db$seg.start <- as.numeric(db$seg.start)
db$seg.end   <- as.numeric(db$seg.end)
db$mid <- round((db$seg.start + db$seg.end) / 2)

# Create label mapping dataframe
mapX <- data.frame(
  V1 = db$seg.name,
  V2 = db$mid,
  V3 = db$seg.name  # Label text
)

par(mar=c(2,2,2,2),bg = "white")
plot(c(1,800),c(1,800),type="n", axes=FALSE,xlab="",ylab="",main="")

#chromossome names

#in here try to change the position of the labels
circos(R=400, cir=db, W=0, mapping= mapX, type="label2",side="out", col= moreGray, cex=1);
	
circos(R=390,cir=db,type="chr",col= colors,print.chr.lab=FALSE,W=20,scale=FALSE)

circos(R=332,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=330,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=327,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=322,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=317,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)


#MIR
circos(R=315, cir=db, W=20, mapping=miR , type="s2", col="#CC5EDEBF",B=F,cex =1);
	
#MIR LABLES
circos(R=330, cir=db, W=20, mapping= miR, type="label",side="out", col=moreGray , cex=0.25);
#21PHAS 

circos(R=262,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=257,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=252,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=247,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)


circos(R=243, cir=db, W=20, mapping= PHAS21lab , type="s2", col="#0095FEBF",B=F,cex = 1);
#21PHAS LABELS

#24PHAS
circos(R=192,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=187,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=182,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)
circos(R=177,cir=db,type="chr",col= moreGray,print.chr.lab=FALSE,W=20,scale=FALSE)

circos(R=174, cir=db, W=20, mapping= PHAS24 , type="s2", col="#DA7724BF",B=F,cex = 1);
