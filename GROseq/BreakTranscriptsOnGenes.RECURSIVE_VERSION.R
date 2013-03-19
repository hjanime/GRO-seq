##
##  Break transcripts on genes.
##
##  The point of this is to condition transcript predictions on
##  known gene annotations.  To do this, transcript predictions
##  are used as is, except when they fully overlap multiple well
##  annotated genes.  In this case, predictions are split based
##  on gene annotations.
##


## Read input files.
T <- read.table("../Hmm-ADC.B-250.V10.bed")
G <- read.table("../fb.ensG.knG.refG.bed")

## Genes should already be ordered, but just in case.
#G <- G[order(G[[2]]),] ## Assume that genes are already ordered.

## Split only on genes >1000bp; Only want to split up pretty major transcrpitional units.
G <- G[which((G[[3]]-G[[2]]) > 1000),]

## Run this function until transcript numbers no longer change.
sepTranscripts <- function(G, T) {

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
for(i in 1:NROW(T)) {
#  print(i)
  ## Identify the indies of genes that are inside the transcript.
  INDX <- which(chrG == chrT[i] & strandG == strandT[i] & endG > startT[i] & startG < endT[i]) 
	## Any transcripts that overlap at all.

  ## Foreach index, identify how much overlap.
  INDX2 <- NULL
  for(j in INDX) {
    OVERLAP <- (min(endT[i], endG[j]) - max(startT[i], startG[j]))/ (endG[j]-startG[j]) ## Overlap is in terms of the gene.

    if(OVERLAP >= TH)
      INDX2 <- c(INDX2, j)
  }

  if(NROW(INDX2) < 2) {
    bT <- rbind(bT, data.frame(c= chrT[i], s= as.integer(startT[i]), e= as.integer(endT[i]), 
				n= paste(chrT[i], as.integer(startT[i]), strandT[i], sep=""), sc= 1, st= strandT[i]))
            ## Write out unaltered
  }
  else { ## Write out two.
    ## Write out first transcript
    if(strandT[i] == "+") {
       ThreePrimeG <- startG[INDX2[2]]-2
    }
    else {
	if(strandT[i] == "-") {
		ThreePrimeG <- endG[INDX2[1]]
	}
	else {
		print("WARNING!  STrand not specified!")
	}
    }
    bT <- rbind(bT, data.frame(c= chrT[i], s= as.integer(startT[i]), e= as.integer(ThreePrimeG),#as.integer(startG[INDX2[2]]-1), 
		n= paste(chrT[i], as.integer(startT[i]), strandT[i], sep=""), sc= 1, st= strandT[i]))

    ## Write out 2nd transcript transcript
     if(strandT[i] == "+") {
        STARTg <- as.integer(startG[INDX2[2]])
     }
     else {
        if(strandT[i] == "-") {
                STARTg <- as.integer(endG[INDX2[1]]+2)
        }
        else {
                print("WARNING!  STrand not specified!")
        }
    }

    bT <- rbind(bT, data.frame(c= chrT[i], s= as.integer(STARTg), e= as.integer(endT[i]), 
		n= paste(chrT[i], as.integer(STARTg), strandT[i], sep=""), sc= 1, st= strandT[i]))
  }

}

return(bT)
}

## Break off 1 transcript at each step; repeat for multiple steps.
nPrevRow <- NROW(T)
nCurrRow <- 1
i        <- 0
bT <- T
while(nPrevRow != nCurrRow) {
  nPrevRow <- nCurrRow
  bT <- sepTranscripts(G,bT)
  nCurrRow <- NROW(bT)

  print(paste("Done with iteration", i, "PrevNum=", nPrevRow, "CurrentNum=", nCurrRow))
  i <- i+1
}
print(head(bT))

## Paste together transcripts that are entirely inside of a single annotation??

write.table(bT, "BrokenTranscripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
save.image("BrokenTranscripts.RData")


