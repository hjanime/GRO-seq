require(GROseq)
sz <- 50

## Load data
#load("../SOAP.0-160.combined.RData")
load("../AllDataInOneVar.RData")
E2alldata <- E2alldata[grep("chr", E2alldata[[1]]),]
print(head(E2alldata))
print(NROW(E2alldata))

## Run the HMM.
DOIT <- function(E2alldata, B, V, PRE) {
	BothStrands <- DetectTranscriptsEM(E2alldata, LtProbB=B, UTS=V, thresh=1, debug=TRUE)

	## Write out hg18 transcripts.
	write.table(BothStrands[[4]], file=pipe(paste("gzip -c > bedCG/Hmm-", PRE, ".B", B, ".V", V, ".bed.gz", sep="")),
		quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

## Run the transition from a transcript --> no transcript over a grid.
for(V in c(10, 20, 15, 5)) {
#  for(tB in (c(-250, -300, -500))) {
  for(tB in (c(-200, -150, -100))) {
	print(paste("Doing all adata in one var; tB=",tB,"V=",V))
	DOIT(E2alldata, tB, V,"ADC")

  }
}

