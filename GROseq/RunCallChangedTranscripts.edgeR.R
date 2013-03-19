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
                                        pvalThresh=0.001, expCounts=1, G=Gbed, debug=TRUE) {
        ## Making a DGE list object will kill those with 0 counts -- pre-emptive elimination.
        nZI <- !(RefSeqNH00 == 0 & RefSeqLC00 == 0 & RefSeqNH10 == 0 & RefSeqLC10 == 0)
#       nZI <- rep(TRUE, NROW(RefSeqNH00))

        ## Use edgeR to determine p-values
        require(edgeR)
        if(debug) print("Building DGEList Object.")

        ## Create edgeR list object
        d <- list()
        d$counts <- as.matrix(data.frame(       NH00= RefSeqNH00[nZI],
                                                LC00= RefSeqLC00[nZI],
                                                NHE2= RefSeqNH10[nZI],
                                                LCE2= RefSeqLC10[nZI]))
        dim(d$counts) <- c(NROW(d$counts), NCOL(d$counts))
## Setting rownames to RefSeqIDs crashes b/c duplicate genes have same RefSeqID
#       rownames(d$counts) <- paste(G[[1]],":",G[[2]],":",G[[4]], sep="")[nZI] #as.character(G[[5]][nZI])
        colnames(d$counts) <- c("NH00", "LC00", "NHE2", "LCE2")
#       print(NROW(d$counts))
#       print(NROW(unique(rownames(d$counts))))
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
        lEXP   <- log( (((RefSeqNH00+RefSeqLC00)/(nrNH00+nrLC00)) + ((RefSeqNH10+RefSeqLC10)/(nrNH10+nrLC10)))/2*expCounts, 10)
        lFC    <- log( ((RefSeqNH10+RefSeqLC10)/(nrNH10+nrLC10))/ ((RefSeqNH00+RefSeqLC00)/(nrNH00+nrLC00)),10)

        ## Make a plot (to check whether we believe the model-based approach).
        if(!is.null(FILENAME)) {
                png(paste(FILENAME,".repBoth.png", sep=""), height=750, width=750)
                  par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
                  plot(log(((RefSeqNH00/nrNH00)+(RefSeqLC00/nrLC00))/2*expCounts,10),
                        log((RefSeqNH00/nrNH00)/(RefSeqLC00/nrLC00),10), ylim=c(-2,2),
                        cex.axis=3.5, cex.lab=4, lwd=3, pch=19, cex=1,
                        xlab="Expression Level (log 10)", ylab="Fold-Change (log 10)")
                  points(log(((RefSeqNH10/nrNH10)+(RefSeqLC10/nrLC10))/2*expCounts,10),
                        log((RefSeqNH10/nrNH10)/(RefSeqLC10/nrLC10),10), pch=19, cex=1)
                  abline(h=0,col="black", lwd=3,lty="dotted")
                dev.off()

                png(paste(FILENAME,".rep10.png", sep=""), height=750, width=750)
                  par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
                  plot(log(((RefSeqNH10/nrNH10)+(RefSeqLC10/nrLC10))/2*expCounts,10),
                        log((RefSeqNH10/nrNH10)/(RefSeqLC10/nrLC10),10), ylim=c(-2,2),
                        cex.axis=3.5, cex.lab=4, lwd=3, pch=19, cex=1,
                        xlab="Expression Level (log 10)", ylab="Fold-Change (log 10)")
                  abline(h=0,col="black", lwd=3,lty="dotted")
                dev.off()

                png(paste(FILENAME, sep=""), height=750, width=750)
                  par(font=3, font.lab=2, font.axis=2, mgp=c(5.5,2,0), mar=c(8.5, 8.5, 2, 2) + 0.3)
                  plot(lEXP,lFC, cex.axis=3.5, cex.lab=4, lwd=3, pch=19, cex=1, ylim=c(-2,2),
                        xlab="Expression Level (log 10)", ylab="Fold-Change (log 10)")
                  abline(h=0,col="black", lwd=3,lty="dotted")
                  points(lEXP[CNG],lFC[CNG], col="red", lwd=2, pch=19, cex=1)
                dev.off()
        }

        ## Return the data in that nice table format.
        ChangedGenes <- data.frame(G[CNG,], FC=(10^META40[CNG])/(10^META00[CNG]),
                        C1=((RefSeqNH00[CNG]+RefSeqLC00[CNG])*expCounts),
                        C2=((RefSeqNH10[CNG]+RefSeqLC10[CNG])*expCounts), pVal=(pval[CNG]) )

        return(ChangedGenes)
}



###########################################################
##################################### 0-10 minute analysis.
###########################################################

# 0-10 minute analysis.
require(GROseq)
source("FitModel.R")

print("Reducing to +1kb to +13kb")
Gbed <- read.table("../transcriptHMM/bedCG/MakeFinalVersion/FinalTranscripts.bed")
tID <- paste(Gbed[[1]],Gbed[[2]],Gbed[[6]], sep="")
Gbed <- data.frame(Chr=as.character(Gbed[[1]]), Start=as.integer(Gbed[[2]]), End=as.integer(Gbed[[3]]), 
			Str=as.character(Gbed[[6]]), RefSeqID=as.character(tID), MGI=as.character(tID))
print(head(Gbed))
G <- LimitToXkb(Gbed, SIZE=13000)

## Count the number of mapped reads.
print("Counting the number of mapped reads...")
load("../SOAP.0-160.NH.LC.RData")
nrNH00 <- NROW(NH00)
nrNH10 <- NROW(NH10)
nrNH40 <- NROW(NH40)
nrNH160<- NROW(NH160)
nrLC00 <- NROW(LC00)
nrLC10 <- NROW(LC10)
nrLC40 <- NROW(LC40)
nrLC160<- NROW(LC160)

# Count reads in +1 to +13kb; for the timecourse.
#load("../SOAP.0-160.NH.LC.noRNAreps.RData")

RefSeqNH00 <- CountReadsInInterval(f= G, p= NH00[,c(1:3,6)])
RefSeqLC00 <- CountReadsInInterval(f= G, p= LC00[,c(1:3,6)])
RefSeqNH10 <- CountReadsInInterval(f= G, p= NH10[,c(1:3,6)])
RefSeqLC10 <- CountReadsInInterval(f= G, p= LC10[,c(1:3,6)])

#Call changed genes.
ChangedGenes <- edgeRFindChangedGenes(RefSeqNH00,nrNH00,
					RefSeqLC00,nrLC00, 
					RefSeqNH10,nrNH10, 
					RefSeqLC10,nrLC10, FILENAME="T/T.ModelFit-10m.png", G=Gbed)

print("Number of Unique Genes 10 at Min")
print(paste("up:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC > 1)]))))
print(paste("dn:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC < 1)]))))

write.table(ChangedGenes, row.names=FALSE, quote=FALSE, file= "T/T.E2.10M.Changed_1-13kb.tsv", sep="\t")

###########################################################
##################################### 0-40 minute analysis.
###########################################################

#print("Reducing to +1kb to +40kb")
#load("../Human.RefSeqGenes.Mar2006.BED.RData")
#G <- LimitToXkb(Gbed, SIZE=40000)

# Count reads in +1 to +40kb.
RefSeqNH00 <- CountReadsInInterval(f= G, p= NH00[,c(1:3,6)])
RefSeqNH40 <- CountReadsInInterval(f= G, p= NH40[,c(1:3,6)])
RefSeqLC00 <- CountReadsInInterval(f= G, p= LC00[,c(1:3,6)])
RefSeqLC40 <- CountReadsInInterval(f= G, p= LC40[,c(1:3,6)])

# Call changed genes.
ChangedGenes <- edgeRFindChangedGenes(RefSeqNH00,nrNH00,
					RefSeqLC00,nrLC00, 
					RefSeqNH40,nrNH40, 
					RefSeqLC40,nrLC40, FILENAME="T/T.ModelFit-40m.png", G=Gbed)

#print("Number of Unique Genes 40 at Min")
print(paste("up:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC > 1)]))))
print(paste("dn:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC < 1)]))))

write.table(ChangedGenes, row.names=FALSE, quote=FALSE, file= "T/T.E2.40M.Changed_1-13kb.tsv", sep="\t")

###########################################################
##################################### 0-160 minute analysis.
###########################################################

#print("Reducing to +1kb to +160kb")
#load("../Human.RefSeqGenes.Mar2006.BED.RData")
#G <- LimitToXkb(Gbed, SIZE=160000)

# Count reads in +1 to +160kb.
RefSeqNH00 <- CountReadsInInterval(f= G, p= NH00[,c(1:3,6)])
RefSeqLC00 <- CountReadsInInterval(f= G, p= LC00[,c(1:3,6)])
RefSeqNH160 <- CountReadsInInterval(f= G, p= NH160[,c(1:3,6)])
RefSeqLC160 <- CountReadsInInterval(f= G, p= LC160[,c(1:3,6)])

#Call changed genes.
ChangedGenes <- edgeRFindChangedGenes(RefSeqNH00,nrNH00,
					RefSeqLC00,nrLC00, 
					RefSeqNH160,nrNH160, 
					RefSeqLC160,nrLC160, FILENAME="T/T.ModelFit-160m.png", G=Gbed)

print("Number of Unique Genes 160 at Min")
print(paste("up:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC > 1)]))))
print(paste("dn:",NROW(unique(ChangedGenes$MGI[which(ChangedGenes$FC < 1)]))))

write.table(ChangedGenes, row.names=FALSE, quote=FALSE, file= "T/T.E2.160M.Changed_1-13kb.tsv", sep="\t")

source("T.RunTimecourse_edgeR.R")
######################################################################################################################
######################################################################################################################

