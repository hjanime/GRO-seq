###########################################################
##################################### Timecourse analysis.
###########################################################

# 0-10 minute analysis.
require(GROseq)
source("../../E2Changed/FitModel.R")

## Analyze read counts in the background of "AllTranscripts".
print("Reducing to +1kb to +13kb")
Gbed <- read.table("../../transcriptHMM/bedCG/MakeFinalVersion/FinalTranscripts.bed")
tID <- paste(Gbed[[1]],Gbed[[2]],Gbed[[6]], sep="")
Gbed <- data.frame(Chr=as.character(Gbed[[1]]), Start=as.integer(Gbed[[2]]), End=as.integer(Gbed[[3]]),
                        Str=as.character(Gbed[[6]]), RefSeqID=as.character(tID), MGI=as.character(tID))
print(head(Gbed))
G <- LimitToXkb(Gbed, SIZE=13000)

## Add in enhancer transcripts.
at <- read.table("../E2.Timecourse.Transcription.Annotated.DIVtoENHasENH.tsv", header=T)
et <- at[at$TYPE=="ENH",]

G <- rbind(data.frame(Chr=as.character(et[[1]]), Start=as.integer(et[[2]]), End=as.integer(et[[3]]),
                        Str=as.character(et[[4]]), ID=as.character(et[[5]]), TYPE=rep("ENH", NROW(et))), data.frame(G, TYPE=rep("OTHER")))
remove(at,et)

# Load data, and calculate numbers ...
load("../../SOAP.0-160.NH.LC.RData")
nrNH00 <- NROW(NH00)
nrNH10 <- NROW(NH10)
nrNH40 <- NROW(NH40)
nrNH160<- NROW(NH160)
nrLC00 <- NROW(LC00)
nrLC10 <- NROW(LC10)
nrLC40 <- NROW(LC40)
nrLC160<- NROW(LC160)
load("../../SOAP.0-160.NH.LC.noRNAreps.RData")

# Count 0 and 40 minutes.
RefSeqNH00 <- CountReadsInInterval(f= G, p= NH00[,c(1:3,6)])
RefSeqNH10 <- CountReadsInInterval(f= G, p= NH10[,c(1:3,6)])
RefSeqNH40 <- CountReadsInInterval(f= G, p= NH40[,c(1:3,6)])
RefSeqNH160 <- CountReadsInInterval(f= G, p= NH160[,c(1:3,6)])

RefSeqLC00 <- CountReadsInInterval(f= G, p= LC00[,c(1:3,6)])
RefSeqLC10 <- CountReadsInInterval(f= G, p= LC10[,c(1:3,6)])
RefSeqLC40 <- CountReadsInInterval(f= G, p= LC40[,c(1:3,6)])
RefSeqLC160 <- CountReadsInInterval(f= G, p= LC160[,c(1:3,6)])

## Calculate expected counts.
require(edgeR)

#nz <- which(RefSeqNH00 > 0 & RefSeqNH10 > 0 & RefSeqNH40 > 0 & RefSeqNH160 > 0 & 
#			RefSeqLC00 > 0 & RefSeqLC10 > 0 & RefSeqLC40 > 0 & RefSeqLC160>0)

## Create edgeR list object
d <- list()
d$counts <- as.matrix(data.frame(       NH00= RefSeqNH00,
					LC00= RefSeqLC00,  
					NH10= RefSeqNH10,
					LC10= RefSeqLC10,
			      		NH40= RefSeqNH40,
					LC40= RefSeqLC40, 
					NH160=RefSeqNH160,
					LC160=RefSeqLC160))

dim(d$counts) <- c(NROW(d$counts), NCOL(d$counts))
colnames(d$counts) <- c("NH00", "LC00", "NH_E2_10m", "LC_E2_10m", "NH_E2_40m", "LC_E2_40m", "NH_E2_160m", "LC_E2_160m")

d$samples$files <- c("VEH", "VEH", "E2_10m", "E2_10m", "E2_40m", "E2_40m", "E2_160m", "E2_160m")
d$samples$group <- as.factor(c("VEH", "VEH", "E2_10m", "E2_10m", "E2_40m", "E2_40m", "E2_160m", "E2_160m")) 
d$samples$description <- c("VEH", "VEH", "E2_10m", "E2_10m", "E2_40m", "E2_40m", "E2_160m", "E2_160m") 
d$samples$lib.size <- c(nrNH00, nrLC00, nrNH10, nrLC10, nrNH40, nrLC40, nrNH160, nrLC160)
d <- new("DGEList",d)
d <- estimateCommonDisp(d)

# Get the fraction of each over 00 minute.  -- do we really need to? 
exp <- data.frame(G, (d$conc$conc.group*d$common.lib.size)[,c(4,1,3,2)])
exp <- exp[exp$TYPE == "ENH",]

write.table(exp, "ET.E2.Timecourse.tsv", row.names=FALSE, quote=FALSE, sep="\t")

# Filter changed genes...
# Get a list of unique RefSeq IDs that have changes in at least 1 timepoint.
c10M <- read.table("ET.E2.10M.FDR001.tsv", header=TRUE)
c40M <- read.table("ET.E2.40M.FDR001.tsv", header=TRUE)
c160M <- read.table("ET.E2.160M.FDR001.tsv", header=TRUE)
E2RegRefSeq <- unique(c(as.character(c10M$ID), as.character(c40M$ID), as.character(c160M$ID)))
remove(c10M)
remove(c40M)
remove(c160M)

INDX <- match(E2RegRefSeq, exp$ID)
expE2Reg <- exp[INDX,]
print("Number of unique MGI IDs regulated during at least 1 timepoint:")
print(NROW(unique(expE2Reg$ID)))
write.table(expE2Reg, "ET.E2Reg.Timecourse-edgeR.tsv", row.names=FALSE, quote=FALSE, sep="\t")

## These have to be separate!
#source("F.MakeMatrix.R")
