## Load RefSeq TSS.#########################################################################
load("Human.RefSeqGenes.Mar2006.BED.RData")
c_tss_indx <- rep(0,NROW(Gbed))
c_tss_indx[which(Gbed[[4]] == "+")] <- 2
c_tss_indx[which(Gbed[[4]] == "-")] <- 3
c_tss <- unlist(lapply(c(1:NROW(Gbed)), function(x) {
	Gbed[x,c_tss_indx[x]]
}))

# Each unique start site has a unique PosID.
#   Non-unique start sites will have the same PosID, and can be consolidated...
PosID= paste(Gbed[[1]], c_tss, Gbed[[4]], sep="")

# Calculate index of first unique tss.
INDX_FirstUnique <- unlist(lapply(unique(PosID), function(x) { which(PosID==x)[1] }))
FirstUnique <- rep(FALSE, NROW(PosID))
FirstUnique[INDX_FirstUnique] <- TRUE

# Put it all together and remove temporary variables.
# For all unique TSS, refseq_tss[refseq_tss$FirstUnique.
# Otherwise, the indices are the same as for Gbed.
refseq_tss <- data.frame(
	Chr= Gbed[[1]], TSS= c_tss, Str= Gbed[[4]], RefSeqID= Gbed[[5]], MGI= Gbed[[6]],
	PosID, FirstUnique
	)
remove(c_tss_indx)
remove(c_tss)
remove(INDX_FirstUnique)
remove(FirstUnique)
remove(PosID)

# Only keep unique start sites.
refseq_tss <- refseq_tss[refseq_tss$FirstUnique,] ## NOT YET!

# Remove chr*_random and chr*_hap*
KEEPINDX <- rep(TRUE, NROW(refseq_tss))
KEEPINDX[grep("random", refseq_tss$Chr)] <- FALSE
KEEPINDX[grep("hap", refseq_tss$Chr)] <- FALSE
KEEPINDX[grep("chrY", refseq_tss$Chr)] <- FALSE
refseq_tss <- refseq_tss[KEEPINDX,]

## Need to make a reverse strand array. #############################################
refseq_tss_rev <- refseq_tss
refseq_tss_rev$Str[refseq_tss[[3]] == "-"] <- "+"
refseq_tss_rev$Str[refseq_tss[[3]] == "+"] <- "-"

## Read cluster 40.  #############################################
require(GROseq)
load("SOAP.0-160.combined.noRNAreps.RData")
ExpReads <- 15000000

## Make meta gene profiles specifically for paused genes. ###########################
# Foreach cluster.
for(M in c(4,8)) {
CLUSTER40 <- read.table(paste("../G001/CLUSTER.",M,".Deviations", sep=""))
for(i in 1:M) {
	CLUSTER1 <- unlist(lapply(unique(as.character(CLUSTER40[[3]][CLUSTER40[[1]] == i])), function(x) { which(as.character(refseq_tss$RefSeqID)==x)[1] }))
	CLUSTER1 <- CLUSTER1[!is.na(CLUSTER1)]

	C10m_Plus <- MetaGene(f=refseq_tss[CLUSTER1,], p=E20m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C10m_Minus <- MetaGene(f=refseq_tss_rev[CLUSTER1,], p=E20m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C10m_Plus <- C10m_Plus/ NROW(CLUSTER1)
	C10m_Minus <- C10m_Minus/ NROW(CLUSTER1)
	## Normalize to expected read counts
	C10m_Plus <- C10m_Plus*ExpReads/NROW(E20m)
	C10m_Minus <- C10m_Minus*ExpReads/NROW(E20m)

	C110m_Plus <- MetaGene(f=refseq_tss[CLUSTER1,], p=E210m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C110m_Minus <- MetaGene(f=refseq_tss_rev[CLUSTER1,], p=E210m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C110m_Plus <- C110m_Plus/ NROW(CLUSTER1)
	C110m_Minus <- C110m_Minus/ NROW(CLUSTER1)
        ## Normalize to expected read counts
        C110m_Plus <- C110m_Plus*ExpReads/NROW(E210m)
        C110m_Minus <- C110m_Minus*ExpReads/NROW(E210m)

	C140m_Plus <- MetaGene(f=refseq_tss[CLUSTER1,], p=E240m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C140m_Minus <- MetaGene(f=refseq_tss_rev[CLUSTER1,], p=E240m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C140m_Plus <- C140m_Plus/ NROW(CLUSTER1)
	C140m_Minus <- C140m_Minus/ NROW(CLUSTER1)
        ## Normalize to expected read counts
        C140m_Plus <- C140m_Plus*ExpReads/NROW(E240m)
        C140m_Minus <- C140m_Minus*ExpReads/NROW(E240m)

	C1160m_Plus <- MetaGene(f=refseq_tss[CLUSTER1,], p=E2160m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C1160m_Minus <- MetaGene(f=refseq_tss_rev[CLUSTER1,], p=E2160m[,c(1:3,6)], size=100, up= 10000, down= 10000)
	C1160m_Plus <- C1160m_Plus/ NROW(CLUSTER1)
	C1160m_Minus <- C1160m_Minus/ NROW(CLUSTER1)
        ## Normalize to expected read counts
        C1160m_Plus <- C1160m_Plus*ExpReads/NROW(E2160m)
        C1160m_Minus <- C1160m_Minus*ExpReads/NROW(E2160m)

	# Determine max and min
	MAX <- max(C10m_Plus, C110m_Plus, C140m_Plus, C1160m_Plus)
	MIN <- -1*max(C10m_Minus, C110m_Minus, C140m_Minus, C1160m_Minus)

	POSITION <- c(-10000:+10000)
	png(paste("CLUSTER/Cluster",i,"-",M,".E2_000m.png", sep=""), width=750, height=750)
	 par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
	 plot(POSITION, C10m_Plus, col="red", type="h", cex.axis=3.5, cex.lab=4, lwd=3, xlim=c(-5000, 5000), xaxt='n', ylim=c(MIN,MAX), ylab="Number of reads per gene.", xlab="Position (relative to TSS).")
	 points(POSITION, (-1*rev(C10m_Minus)), col="blue", type="h")
	 abline(mean(C10m_Plus[5000:8000]), 0, lty="dotted")
	 axis(side=1, at=c(-4000, 0, 4000), cex.axis=3.5)
	dev.off()

	png(paste("CLUSTER/Cluster",i,"-",M,".E2_010m.png", sep=""), width=750, height=750)
	 par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
	 plot(POSITION, C110m_Plus, col="red", type="h", cex.axis=3.5, cex.lab=4, lwd=3, xlim=c(-5000, 5000), ylim=c(MIN,MAX), ylab="Number of reads per gene.", xlab="Position (relative to TSS).", xaxt='n')
	 points(POSITION, (-1*rev(C110m_Minus)), col="blue", type="h")
	 abline(mean(C110m_Plus[5000:8000]), 0, lty="dotted")
         axis(side=1, at=c(-4000, 0, 4000), cex.axis=3.5)
	dev.off()

	png(paste("CLUSTER/Cluster",i,"-",M,".E2_040m.png", sep=""), width=750, height=750)
	 par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
 	 plot(POSITION, C140m_Plus, col="red", type="h", cex.axis=3.5, cex.lab=4, lwd=3, xlim=c(-5000, 5000), ylim=c(MIN,MAX), ylab="Number of reads per gene.", xlab="Position (relative to TSS).", xaxt='n')
	 points(POSITION, (-1*rev(C140m_Minus)), col="blue", type="h")
	 abline(mean(C140m_Plus[5000:8000]), 0, lty="dotted")
         axis(side=1, at=c(-4000, 0, 4000), cex.axis=3.5)
	dev.off()

	png(paste("CLUSTER/Cluster",i,"-",M,".E2_160m.png", sep=""), width=750, height=750)
	 par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
	 plot(POSITION, C1160m_Plus, col="red", type="h", cex.axis=3.5, cex.lab=4, lwd=3, xlim=c(-5000, 5000), ylim=c(MIN,MAX), ylab="Number of reads per gene.", xlab="Position (relative to TSS).", xaxt='n')
	 points(POSITION, (-1*rev(C1160m_Minus)), col="blue", type="h")
	 abline(mean(C1160m_Plus[5000:8000]), 0, lty="dotted")
         axis(side=1, at=c(-4000, 0, 4000), cex.axis=3.5)
	dev.off()

        png(paste("CLUSTER/Cluster",i,"-",M,".AllTraces.png", sep=""), width=750, height=750)
         par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)

         plot(POSITION, C10m_Plus, col="red", type="l", cex.axis=3.5, cex.lab=4, lwd=3, 
		xlim=c(-5000, 5000), ylim=c(MIN,MAX), ylab="Number of reads per gene.", 
		xlab="Position (relative to TSS).", xaxt='n')
         points(POSITION, (-1*rev(C10m_Minus)), col="red", type="l", lwd=3)

         points(POSITION, C110m_Plus, col="orange", type="l", lwd=3)
         points(POSITION, (-1*rev(C110m_Minus)), col="orange", type="l", lwd=3)

         points(POSITION, C140m_Plus, col="green", type="l", lwd=3)
         points(POSITION, (-1*rev(C140m_Minus)), col="green", type="l", lwd=3)

         points(POSITION, C1160m_Plus, col="blue", type="l", lwd=3)
         points(POSITION, (-1*rev(C1160m_Minus)), col="blue", type="l", lwd=3)

         axis(side=1, at=c(-4000, 0, 4000), cex.axis=3.5)

        dev.off()


}

}
