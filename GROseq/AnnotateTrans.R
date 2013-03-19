## Read in main E2reg file.
E2 <- read.table("E2.Timecourse.tsv",header=T)
TYPE <- rep("NA", NROW(E2))
CLASS <- rep("NA", NROW(E2))

## Read in annotation files...
PATH <- "/usr/projects/GROseq/SOAP/transcriptHMM/AnnotatePredictions/TMP/"

## Refseq Genes
refGenes <- read.table(paste(PATH,"refGene.Overlap.bed", sep=""))
iRG <- match(E2[[5]], refGenes[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "REF"
CLASS[!is.na(iRG)] <- as.character(refGenes[iRG[!is.na(iRG)],10])
sum(TYPE == "REF")

## RNA Genes
rnaGenes <- read.table(paste(PATH,"rnaGene.Overlap.bed", sep=""))
iRG <- match(E2[[5]], rnaGenes[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "RNA"
CLASS[!is.na(iRG)] <- as.character(rnaGenes[iRG[!is.na(iRG)],10])
sum(TYPE == "RNA")

## Antisense Transc.
ansTrans <- read.table(paste(PATH,"as-refGene.Overlap.bed", sep=""))
iRG <- match(E2[[5]], ansTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "AS"
CLASS[!is.na(iRG)] <- as.character(ansTrans[iRG[!is.na(iRG)],10])
sum(TYPE == "AS")

## Divergent transcripts
divTrans <- read.table(paste(PATH,"divergent.Overlap.bed", sep=""))
iRG <- match(E2[[5]], divTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "DIV"
CLASS[!is.na(iRG)] <- as.character(divTrans[iRG[!is.na(iRG)],10])
sum(TYPE == "DIV")

## Enhancer Transc.
enhTrans <- read.table(paste(PATH,"E2PeaksOverlap.bed", sep=""))
iRG <- match(E2[[5]], enhTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "ENH"
CLASS[!is.na(iRG)] <- "ENH"
sum(TYPE == "ENH")

## Sense refgene bad match.
senseGeneBad <- read.table(paste(PATH, "refGeneBadMatch.Overlap.bed", sep=""))
iRG <- match(E2[[5]], senseGeneBad[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "GENE_BadMatch"
CLASS[!is.na(iRG)] <- as.character(senseGeneBad[iRG[!is.na(iRG)],10])
sum(TYPE == "GENE_BadMatch")

## Antisense refgene bad match.
asGeneBad    <- read.table(paste(PATH, "refGeneASBadMatch.Overlap.bed", sep=""))
iRG <- match(E2[[5]], asGeneBad[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "AS_BadMatch"
CLASS[!is.na(iRG)] <- as.character(asGeneBad[iRG[!is.na(iRG)],10])
sum(TYPE == "AS_BadMatch")

## Intergenic transcription
othTrans <- read.table(paste(PATH,"refASBadMatch.NoOverlap.bed", sep=""))
iRG <- match(E2[[5]], othTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "INTERGENIC"
CLASS[!is.na(iRG)] <- as.character(othTrans[iRG[!is.na(iRG)],4])
sum(TYPE == "INTERGENIC")

AN <- cbind(E2, TYPE= as.character(TYPE), CLASS= as.character(CLASS))
write.table(AN, "E2.Timecourse.Transcription.Annotated.tsv", sep="\t", row.names=FALSE, col.names=T, quote=FALSE)

###########
# Genes coated with rna...
##
## For these, I will have to make new columns.
Eout <- NULL
TYPE <- NULL
CLASS <- NULL

rnaGenes <- read.table(paste(PATH,"rnaGene.Overlap.bed", sep=""))
iRG <- match(rnaGenes[[4]], E2[[5]])

## Ensure that there are no duplicates...
rID <- paste(rnaGenes[[7]],rnaGenes[[8]],rnaGenes[[12]],sep="")
if(NROW(rnaGenes) != NROW(unique(rID))) {
  print("WARNING: POTENTIAL DUPLICATE RNAs IN DATASET: Distinct RNA Genes.")
}

for(i in c(1:NROW(iRG))) {
    Eout <- rbind(Eout, E2[iRG[i],])
    TYPE <- c(TYPE, "RNA")
    CLASS <- c(CLASS, as.character(rnaGenes[i,10]))
}

ArnaGenes <- read.table(paste(PATH,"rnaGene.Genic.Overlap.bed", sep=""))
iRG <- match(ArnaGenes[[4]], E2[[5]])

## Ensure that there are no duplicates...
rID <- paste(ArnaGenes[[7]],ArnaGenes[[8]],ArnaGenes[[12]],sep="")
if(NROW(ArnaGenes) != NROW(unique(rID))) {
  print("WARNING: POTENTIAL DUPLICATE RNAs IN DATASET: Genic Overlap.")
}

for(i in c(1:NROW(iRG))) {
    Eout <- rbind(Eout, E2[iRG[i],])
    TYPE <- c(TYPE, "RNA_withGene")
    CLASS <- c(CLASS, as.character(ArnaGenes[i,10]))
}

AN <- cbind(Eout, TYPE= as.character(TYPE), CLASS= as.character(CLASS))
write.table(AN, "E2.RNA_Gene.Timecourse.Transcription.Annotated.tsv", sep="\t", row.names=FALSE, col.names=T, quote=FALSE)

