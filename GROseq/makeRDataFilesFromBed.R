# Nasun
NH00 <- read.table("NH.00Min.bed", sep="\t")
NH10 <- read.table("NH.10Min.bed", sep="\t")
NH40 <- read.table("NH.40Min.bed", sep="\t")
NH160 <-read.table("NH.160Min.bed", sep="\t")
# Leighton
LC00 <- read.table("LC.00Min.bed", sep="\t")
LC10 <- read.table("LC.10Min.bed", sep="\t")
LC40 <- read.table("LC.40Min.bed", sep="\t")
LC160 <-read.table("LC.160Min.bed", sep="\t")

# Remove spike-ins
NH00 <- NH00[grep("chr|rRNA", NH00[[1]]),]
NH10 <- NH10[grep("chr|rRNA", NH10[[1]]),]
NH40 <- NH40[grep("chr|rRNA", NH40[[1]]),]
NH160<-NH160[grep("chr|rRNA",NH160[[1]]),]
LC00 <- LC00[grep("chr|rRNA", LC00[[1]]),]
LC10 <- LC10[grep("chr|rRNA", LC10[[1]]),]
LC40 <- LC40[grep("chr|rRNA", LC40[[1]]),]
LC160<-LC160[grep("chr|rRNA",LC160[[1]]),]


# Order all probes.
NH00 <- NH00[order(NH00$V2),]
NH10 <- NH10[order(NH10$V2),]
NH40 <- NH40[order(NH40$V2),]
NH160 <-NH160[order(NH160$V2),]
LC00 <- LC00[order(LC00$V2),]
LC10 <- LC10[order(LC10$V2),]
LC40 <- LC40[order(LC40$V2),]
LC160 <-LC160[order(LC160$V2),]

# Write image
save.image("../SOAP.0-160.NH.LC.RData")

## Join
E20m <- data.frame(	chrom=		c(as.character(NH00[[1]]), as.character(LC00[[1]])), 
			chromStart=	c(NH00[[2]], LC00[[2]]),
			chromEnd=	c(NH00[[3]], LC00[[3]]),
			name=		c(NH00[[4]], LC00[[4]]),
			score=		c(NH00[[5]], LC00[[5]]),
			strand=		c(as.character(NH00[[6]]), as.character(LC00[[6]])) )

E210m <- data.frame(    chrom=          c(as.character(NH10[[1]]), as.character(LC10[[1]])),
                        chromStart=     c(NH10[[2]], LC10[[2]]),
                        chromEnd=       c(NH10[[3]], LC10[[3]]),
                        name=           c(NH10[[4]], LC10[[4]]),
                        score=          c(NH10[[5]], LC10[[5]]),
                        strand=         c(as.character(NH10[[6]]), as.character(LC10[[6]])) )

E240m <- data.frame(    chrom=          c(as.character(NH40[[1]]), as.character(LC40[[1]])),
                        chromStart=     c(NH40[[2]], LC40[[2]]),
                        chromEnd=       c(NH40[[3]], LC40[[3]]),
                        name=           c(NH40[[4]], LC40[[4]]),
                        score=          c(NH40[[5]], LC40[[5]]),
                        strand=         c(as.character(NH40[[6]]), as.character(LC40[[6]])) )

E2160m <- data.frame(    chrom=          c(as.character(NH160[[1]]), as.character(LC160[[1]])),
                        chromStart=     c(NH160[[2]], LC160[[2]]),
                        chromEnd=       c(NH160[[3]], LC160[[3]]),
                        name=           c(NH160[[4]], LC160[[4]]),
                        score=          c(NH160[[5]], LC160[[5]]),
                        strand=         c(as.character(NH160[[6]]), as.character(LC160[[6]])) )

# remove old.
remove(NH00)
remove(NH10)
remove(NH40)
remove(NH160)
remove(LC00)
remove(LC10)
remove(LC40)
remove(LC160)

# Must order again ...
E20m <- E20m[order(E20m[[2]]),]
E210m <- E210m[order(E210m[[2]]),]
E240m <- E240m[order(E240m[[2]]),]
E2160m <-E2160m[order(E2160m[[2]]),]

# Save
save.image("../SOAP.0-160.combined.RData")
remove(E20m)
remove(E210m)
remove(E240m)
remove(E2160m)

#####################################################################
###### Remove tRNA/rRNA
#####################################################################
load("../SOAP.0-160.NH.LC.RData")

tRNA <- read.table("hg18.rnaRepeats.tsv")
tRNA <- data.frame(tRNA[[1]], tRNA[[2]]-44, tRNA[[3]]+44, tRNA[[4]])

require(GROseq)
REMOVE <- AssociateWithInterval(p=data.frame(NH00[[1]], NH00[[2]], NH00[[3]], c(1:NROW(NH00))), f=tRNA)
NH00 <- NH00[is.na(REMOVE),]
NH00 <- NH00[!(as.character(NH00[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(NH10[[1]], NH10[[2]], NH10[[3]], c(1:NROW(NH10))), f=tRNA)
NH10 <- NH10[is.na(REMOVE),]
NH10 <- NH10[!(as.character(NH10[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(NH40[[1]], NH40[[2]], NH40[[3]], c(1:NROW(NH40))), f=tRNA)
NH40 <- NH40[is.na(REMOVE),]
NH40 <- NH40[!(as.character(NH40[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(NH160[[1]], NH160[[2]], NH160[[3]], c(1:NROW(NH160))), f=tRNA)
NH160 <- NH160[is.na(REMOVE),]
NH160 <- NH160[!(as.character(NH160[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(LC00[[1]], LC00[[2]], LC00[[3]], c(1:NROW(LC00))), f=tRNA)
LC00 <- LC00[is.na(REMOVE),]
LC00 <- LC00[!(as.character(LC00[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(LC10[[1]], LC10[[2]], LC10[[3]], c(1:NROW(LC10))), f=tRNA)
LC10 <- LC10[is.na(REMOVE),]
LC10 <- LC10[!(as.character(LC10[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(LC40[[1]], LC40[[2]], LC40[[3]], c(1:NROW(LC40))), f=tRNA)
LC40 <- LC40[is.na(REMOVE),]
LC40 <- LC40[!(as.character(LC40[[1]]) == "rRNA"),]

REMOVE <- AssociateWithInterval(p=data.frame(LC160[[1]], LC160[[2]], LC160[[3]], c(1:NROW(LC160))), f=tRNA)
LC160 <- LC160[is.na(REMOVE),]
LC160 <- LC160[!(as.character(LC160[[1]]) == "rRNA"),]

save.image("../SOAP.0-160.NH.LC.noRNAreps.RData")

## Join
E20m <- data.frame(	chrom=		c(as.character(NH00[[1]]), as.character(LC00[[1]])), 
			chromStart=	c(NH00[[2]], LC00[[2]]),
			chromEnd=	c(NH00[[3]], LC00[[3]]),
			name=		c(NH00[[4]], LC00[[4]]),
			score=		c(NH00[[5]], LC00[[5]]),
			strand=		c(as.character(NH00[[6]]), as.character(LC00[[6]])) )

E210m <- data.frame(    chrom=          c(as.character(NH10[[1]]), as.character(LC10[[1]])),
                        chromStart=     c(NH10[[2]], LC10[[2]]),
                        chromEnd=       c(NH10[[3]], LC10[[3]]),
                        name=           c(NH10[[4]], LC10[[4]]),
                        score=          c(NH10[[5]], LC10[[5]]),
                        strand=         c(as.character(NH10[[6]]), as.character(LC10[[6]])) )

E240m <- data.frame(    chrom=          c(as.character(NH40[[1]]), as.character(LC40[[1]])),
                        chromStart=     c(NH40[[2]], LC40[[2]]),
                        chromEnd=       c(NH40[[3]], LC40[[3]]),
                        name=           c(NH40[[4]], LC40[[4]]),
                        score=          c(NH40[[5]], LC40[[5]]),
                        strand=         c(as.character(NH40[[6]]), as.character(LC40[[6]])) )

E2160m <- data.frame(    chrom=          c(as.character(NH160[[1]]), as.character(LC160[[1]])),
                        chromStart=     c(NH160[[2]], LC160[[2]]),
                        chromEnd=       c(NH160[[3]], LC160[[3]]),
                        name=           c(NH160[[4]], LC160[[4]]),
                        score=          c(NH160[[5]], LC160[[5]]),
                        strand=         c(as.character(NH160[[6]]), as.character(LC160[[6]])) )

# remove old.
remove(NH00)
remove(NH10)
remove(NH40)
remove(NH160)
remove(LC00)
remove(LC10)
remove(LC40)
remove(LC160)

# Must order again ...
E20m <- E20m[order(E20m[[2]]),]
E210m <- E210m[order(E210m[[2]]),]
E240m <- E240m[order(E240m[[2]]),]
E2160m <-E2160m[order(E2160m[[2]]),]

# Save
save.image("../SOAP.0-160.combined.noRNAreps.RData")
