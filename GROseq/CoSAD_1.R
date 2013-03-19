###
###  CorrelateSenseAntisenseDivergent.R
###
###  Correlates sense transcription with divergent and/or antisense transcription.
###

## Set up folders...
DIRap <- "~/tnandu1/SOAP/transcriptHMM/AnnotatePredictions/TMP/"

##############################################################################################
###
###  For divegent transcription, two methods are used:
###  (1) for each sense transcript that changes,
###     (1A) Identify the set of all associated "divergent" transcripts.
###	(1B) Correlate sense and divergent transcripts -- scatterplot/density plot + Pearson
###
##############################################################################################

## Get the set of all TranscriptIDs that are associated with divergent transcripts.
## Take this from an intermdiate step in the pipeline.
DT    <- read.table(paste(DIRap, "divergent.Overlap.bed", sep=""))

## Make a table of divergent, sense pairs.

	## NOTE: Limit the analysis to transcripts that are upstream of detected transcripts NOT annotated,
	## genes as those are not guaranteed to be tested for regulation (i.e. if they are internal).
	## Equilavent to limiting the analysis to the most 5' end of each transcript.  
	## Note that those which are both annotated and de novo are represented multiple times, so
	## we are guaranteed to get all transcripts.

dtID <- DT[grep("chr", DT[[10]]),c(4,10)]
remove(DT)
colnames(dtID) <- c("DIV", "PRI")

## Identify rows where the divergent OR primary transcript are changed in at least 1 timepoint...
E2  <- read.table(paste("E2Reg.Transcription.Annotated.tsv", sep=""), header=T)
dtID <- dtID[((dtID[[1]] %in% E2[[5]]) | (dtID[[2]] %in% E2[[5]])),] ## Limit to transcripts that are in the regulated list...
remove(E2)
print(paste("Number of Primary/Divergent transcript pairs: ",NROW(dtID)))

## Get the timecourse for all transcripts.
TC  <- read.table(paste("E2.Timecourse.tsv", sep=""), header=T)
divINDX <- match(dtID$DIV, TC[[5]])
priINDX <- match(dtID$PRI, TC[[5]])

## Check for NAs.
if( sum(is.na(c(divINDX, priINDX))) > 0 ) 
	print("WARNING: NAs in divINDX | priINDX")

## Build a vector appending each of these...
div <- NULL
pri <- NULL
num <- 0
cond<- 0
for(i in 1:NROW(dtID)) {

  ## Only use points at which there is expression information.
  if((TC[divINDX[i],7] >= 1) & (TC[priINDX[i],7] >= 1)) {
     num <- num+1
   ## Normalize each to the 0 minute timepoint... Thus, what we are plotting is fold-changes over 0 min.
   for(indx in c(8:10)) {
     if((TC[divINDX[i],indx] >= 1) & (TC[priINDX[i],indx] >= 1)) {
       cond <- cond+1
       div <- c(div, as.real(TC[divINDX[i], indx ])/TC[divINDX[i],7])
       pri <- c(pri, as.real(TC[priINDX[i], indx ])/TC[priINDX[i],7])
     }
   }

  }

}
print(paste("The number of non-zero Primary/Divergent pairs:", num, " Conditions:", cond))

## Get the correlation, and plot log-log.
print(cor.test(log(div,2), log(pri,2)))

require(MASS)
png("PrimaryDivergent.Transcript.png", height=750, width=750)
	par(font=3, font.lab=2, font.axis=2, mar=c(5.5, 5.5, 2, 2) + 0.3)
	plot(log(div,2), log(pri,2), col="black", cex=1, pch=16, lwd=3, cex.axis=2, cex.lab=2,
							xlab="Divergent Transcript Fold-Change (log2)", 
							ylab="Primary Transcript Fold-Change (log2)")
	contour(kde2d(log(div,2), log(pri,2),n=250), add=T, lwd=2, col="red")
dev.off()

###  The first one looks good ... this is not a priority!
###
###  (2) for each sense transcript that changes,
###     (2A) Using CountReadsInInterval, count the 500bp interval 5' of the transcript.
###	(2B) Correlate sense and divergent transcripts -- scatterplot/density plot + Pearson
###
###  All of this can (but does not have to be) repeated with RefSeq gene annotations...
###


##############################################################################################
###
###  Antisense transcription is defined according to coding genes.  Thus, this analysis
###  also takes an annotation-centric view.
###  (1) for each genic transcript that changes in response to estrogen,
###     (1A) Identify the set of all associated "antisense" transcripts.
###	(1B) Correlate sense and divergent transcripts -- scatterplot/density plot + Pearson
###
##############################################################################################

### Read annotation file.  Write out separate files for annotated genes and antisense transcription.
at <- read.table(paste("E2.Timecourse.Transcription.Annotated.tsv",sep=""), header=T)
at <- data.frame(at, rep(1,NROW(at)))
write.table(at[which(at$TYPE=="REF"),c(1:3,5,13,4)], "TMP.SenseGenes.bed", row.names=F, quote=F, col.names=F, sep="\t")
write.table(at[which(at$TYPE=="AS"),c(1:3,5,13,4)], "TMP.AntiSenseGenes.bed", row.names=F, quote=F, col.names=F, sep="\t") 
remove(at)

### Create a table of sense/antisense transcript pairs using Kent Tools' overlapSelect
### Note that this command places sense genes first, and corresponding antisense genes second.
RUN <- "overlapSelect -mergeOutput TMP.AntiSenseGenes.bed TMP.SenseGenes.bed  TMP.out"
system(RUN)
dtID <- read.table("TMP.out")[,c(4,10)]
colnames(dtID) <- c("SENSE", "ANTISENSE")

### Now its just like the divergent transcription example ...
### Identify rows where the divergent OR primary transcript are changed in at least 1 timepoint...
E2  <- read.table(paste("E2Reg.Transcription.Annotated.tsv", sep=""), header=T)
dtID <- dtID[((dtID[[1]] %in% E2[[5]]) | (dtID[[2]] %in% E2[[5]])),] ## Limit to transcripts that are in the regulated list...
remove(E2)
print(paste("Number of Primary/Divergent transcript pairs: ",NROW(dtID)))

## Get the timecourse for all transcripts.
TC  <- read.table(paste("E2.Timecourse.tsv", sep=""), header=T)
senseINDX     <- match(dtID$SENSE, TC[[5]])
antisenseINDX <- match(dtID$ANTISENSE, TC[[5]])

## Check for NAs.
if( sum(is.na(c(senseINDX, antisenseINDX))) > 0 ) 
	print("WARNING: NAs in senseINDX | antisenseINDX")

## Build a vector appending each of these...
sense <- NULL
antisense <- NULL
num <- 0
cond<- 0
for(i in 1:NROW(dtID)) {

  ## Only use points at which there is expression information.
  if((TC[senseINDX[i],7] >= 1) & (TC[antisenseINDX[i],7] >= 1)) {
     num <- num+1
   ## Normalize each to the 0 minute timepoint... Thus, what we are plotting is fold-changes over 0 min.
   for(indx in c(8:10)) {
     if((TC[senseINDX[i],indx] >= 1) & (TC[antisenseINDX[i],indx] >= 1)) {
       cond <- cond+1
       sense     <- c(sense,     as.real(TC[senseINDX[i],     indx ])/TC[senseINDX[i],    7])
       antisense <- c(antisense, as.real(TC[antisenseINDX[i], indx ])/TC[antisenseINDX[i],7])
     }
   }

  }

}
print(paste("The number of non-zero Sense/Antisense pairs:", num, " Conditions:", cond))

## Get the correlation, and plot log-log.
print(cor.test(log(sense,2), log(antisense,2)))

require(MASS)
png("PrimaryAntisense.Transcript.png", height=750, width=750)
	par(font=3, font.lab=2, font.axis=2, mar=c(5.5, 5.5, 2, 2) + 0.3)
	plot(log(antisense,2), log(sense,2), col="black", cex=1, pch=16, lwd=3, cex.axis=2, cex.lab=2,
							xlab="Anti-sense Transcript Fold-Change (log2)", 
							ylab="Sense Transcript Fold-Change (log2)")
	contour(kde2d(log(antisense,2), log(sense,2),n=500), add=T, lwd=2, col="red")
dev.off()

