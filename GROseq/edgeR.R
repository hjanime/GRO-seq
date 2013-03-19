##
## This script file is the SOAP analysis.
##
## Changes in expression are called using the edgeR package
## in Bioconductor.  
##
## edgeR uses an empirical Bayes rule to calculate a weight
## parameter, alpha, that describes the varience in the data.
##

## edgeRFindChangedGnees identifies genes that change in short-read
## data using the R package defined above, followed by fdr correction
## for multiple hypothesis testing.
##
edgeRFindChangedGenes <- function(RefSeqNH00,nrNH00,
					RefSeqLC00,nrLC00, 
					RefSeqNH10,nrNH10, 
					RefSeqLC10,nrLC10, FILENAME="ModelFit-10m.png", 
					pvalThresh=20.001, expCounts=15000000, G=Gbed, debug=TRUE) {
	## Making a DGE list object will kill those with 0 counts -- pre-emptive elimination.
	nZI <- !(RefSeqNH00 == 0 & RefSeqLC00 == 0 & RefSeqNH10 == 0 & RefSeqLC10 == 0) 
#	nZI <- rep(TRUE, NROW(RefSeqNH00))

	## Use edgeR to determine p-values
	require(edgeR)
	if(debug) print("Building DGEList Object.")

	## Create edgeR list object
	d <- list()
	d$counts <- as.matrix(data.frame(	NH00= RefSeqNH00[nZI], 
						LC00= RefSeqLC00[nZI],  
                                        	NHE2= RefSeqNH10[nZI], 
						LCE2= RefSeqLC10[nZI]))
	dim(d$counts) <- c(NROW(d$counts), NCOL(d$counts))
## Setting rownames to RefSeqIDs crashes b/c duplicate genes have same RefSeqID
#	rownames(d$counts) <- paste(G[[1]],":",G[[2]],":",G[[4]], sep="")[nZI] #as.character(G[[5]][nZI])
	colnames(d$counts) <- c("NH00", "LC00", "NHE2", "LCE2")
#	print(NROW(d$counts))
#	print(NROW(unique(rownames(d$counts))))
#       print(NROW(unique(as.character(rownames(d$counts)))))

	d$samples$files <- c("NH00", "LC00", "NHE2", "LCE2")
	d$samples$group <- as.factor(c("VEH", "VEH", "E2", "E2")) ## edgeR requires factor here...
	d$samples$description <- c("VEH", "VEH", "E2", "E2") 
	d$samples$lib.size <- c(nrNH00, nrLC00, nrNH10, nrLC10)

	## Create DGEList Object & analyze changes using edgeR
	d <- new("DGEList",d)

	if(debug) print("Estimating Common Dispersion.")
	d <- estimateCommonDisp(d)

	if(debug) print("exactTest.")
	de.com <- exactTest(d)

	if(debug) print("Adjusting p-values using fdr.")
	adj.p <- p.adjust(de.com$table$p.value, "fdr")
	print(summary(adj.p))
	pval <- rep(1, NROW(G))
	pval[nZI] <- adj.p
	CNG <- (pval <= pvalThresh)

	## Make the combined META-read counts for plotting and returning purposes.
	META00 <- log( (RefSeqNH00+RefSeqLC00)/(nrNH00+nrLC00), 10 )
	META40 <- log( (RefSeqNH10+RefSeqLC10)/(nrNH10+nrLC10), 10 )
# Make a plot (to check whether we believe the model-based approach).
	if(!is.null(FILENAME)) {
		png(FILENAME, height=600, width=600)

		plot(META00,META40, xlab="Log(Reads In Gene/ Total) -- 0 Min", ylab="Log(Reads In Gene/ Total) -- X Min")
		points(META00[CNG],META40[CNG], col="red", lwd=2)

		dev.off()
	}

	## Return the data in that nice table format.
	ChangedGenes <- data.frame(G[CNG,], FC=(10^META40[CNG])/(10^META00[CNG]), 
			C1=((RefSeqNH00[CNG]+RefSeqLC00[CNG])*expCounts), 
			C2=((RefSeqNH10[CNG]+RefSeqLC10[CNG])*expCounts), pVal=(pval[CNG]) )

	return(ChangedGenes)
}

######################################################################################################################
######################################################################################################################

