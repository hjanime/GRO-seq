##
## Script combines transcripts that are within the same gene annotation, combining smaller transcripts
## for genes with low regulation into a single transcript representing the gene.
##

## Read input files.
T <- read.table("../Hmm-ADC.B-200.V5.bed")
G <- read.table("../fb.ensG.knG.refG.bed")

## Genes should already be ordered, but just in case.
#G <- G[order(G[[2]]),] ## Assume that genes are already ordered.

## Split only on genes >1000bp; Only want to split up pretty major transcrpitional units.
G <- G[which((G[[3]]-G[[2]]) > 1000),]

# Casting -- makes more human readable.
chrG    <- as.character(G[[1]])
chrT    <- as.character(T[[1]])
startG  <- as.integer(G[[2]])
startT  <- as.integer(T[[2]])
endG    <- as.integer(G[[3]])
endT    <- as.integer(T[[3]])
strandG <- as.character(G[[6]])
strandT <- as.character(T[[6]])

bT <- NULL
TH <- 0.8 ## Threshold for calling the gene part of the transcript
for(i in 1:NROW(G)) {
#  print(i)
  ## Identify the indies of genes that are inside the transcript.
  INDX <- which(chrT == chrG[i] & strandT == strandG[i] & endT > startG[i] & startT < endG[i]) 
	## Any transcripts that overlap at all.

  ## Foreach index, identify how much overlap.
  INDX2 <- NULL
  for(j in INDX) {
    OVERLAP <- (min(endG[i], endT[j]) - max(startG[i], startT[j]))/ (1+endT[j]-startT[j]) ## Overlap is in terms of the transcript -- .
#    print(paste(i, j, endG[i], endT[j], startG[i], startT[j], OVERLAP))
    if(OVERLAP >= TH)
      INDX2 <- c(INDX2, j)
  }

  if(NROW(INDX2) > 1) { ## Write out one transcript; from MIN(starts) to MAX(ends).
    START <- min(startT[INDX2])
    END   <- max(endT[INDX2])
    bT <- rbind(bT, data.frame(c= chrG[i], s= as.integer(START), e= as.integer(END),#as.integer(startG[INDX2[2]]-1), 
		n= paste(chrG[i], as.integer(START), strandG[i], sep=""), sc= 1, st= strandG[i]))
  }

}


## Paste together transcripts that are entirely inside of a single annotation??
print(head(bT))
print(paste("BEFORE=",NROW(T), "AFTER=",NROW(bT)))


# Break off 1 transcript at each step; repeat for multiple steps.
write.table(bT, "CombinedTranscripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
save.image("CombinedTranscripts.RData")


