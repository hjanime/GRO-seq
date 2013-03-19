####  The point of this code is to find the contour, or shell.

###########################################################
###
###	PointsAboveModel -- Identifies the number of points above the model.
###
###	Returns the indices for points above the line.
###	Also calculates the fraction of points above the line (for CI).
###
###	To calculate points below, just flip axes and call again :)
###
###	TODO: 
###	(1) 
###	(2) 
###	(3) 
###
###########################################################

PointsAboveModel <- function(MODEL, MAX00, MIN00) {
	O <- order(MAX00)
	gtModel <- rep(FALSE,NROW(MAX00))
	for(i in 1:NROW(MODEL)) {
		# Collect points that belong to this line segment (P).
		P <- which(MAX00 <= MODEL[i,5] & MAX00 >= MODEL[i,3])

		# Foreach point, is it above?
		gtModelTest <- unlist(lapply(P, function(x) {
			MIN00[x] > MODEL[i,1] + MAX00[x]*MODEL[i,2]
		}))

		# Make them big & red for testing purposes...
#		points(MAX00[P][gtModelTest], MIN00[P][gtModelTest], col="red", cex=5, lwd=2)

		# Assign ...
		gtModel[P] <- gtModelTest
	}

	# Mark them for testing purposes ...
#	points(MAX00[gtModel], MIN00[gtModel], col="green", cex=3, lwd=2)

	return(gtModel)
}

###########################################################
###
###	CalcIntersection -- Calculates intersection between two lines.
###
###########################################################

CalcIntersection <- function(b, m, b1, m1) {
	A <- NULL
	A$x= (b1-b)/(m-m1)
	A$y= (m1*A$x)+b1
	return(A)
}
#CalcIntersection(0,1,1,-1) # Should return A$x = 0.5; A$y = 0.5.


###########################################################
###
###	ModelFit -- Divides Into Equally Spaced Groups.
###
###	Returns a matrix of 2,n.  Columns represent Y-intercept, and slope.
###
###	TODO: 
###	--- Return the model as a matrix.
###	(2) Implement CI calling.
###	--- Smooth out the line, such that there are no jaggies between lines.
###	(4) Detect & Prevent backtrack!
###
###########################################################

ModelFit <- function(NH00, LC00, NGroups= 5, NUMOfGenesInBin= 10, MinTH= 1.5, CI= 0.05) {
	# Shift about the (b=0,m=1) line.
	MAX00 <- unlist(lapply( c(1:NROW(NH00)), function(x) { min(NH00[x], LC00[x]) }))
	MIN00 <- unlist(lapply( c(1:NROW(NH00)), function(x) { max(NH00[x], LC00[x]) }))

	I <- ((NH00 != -Inf) & (LC00 != -Inf))	# Identify the union of genes with > 0 reads in both samples.
	BinSize <- abs(min(MAX00[I])-
		max(MAX00[I]))/NGroups		# Find the number of genes in each bin.
	Sstart <- min(MAX00[I])			# Find the starting bin

	# Find index of top S genes in each bin.  Presently returns the index of points on the shell in 
	#   the vector that has been reduced to [I-space].
	SHELL <- lapply(c(1:NGroups), function(x) {
		i <- which( (((x-1)*BinSize+Sstart) <= MAX00[I]) &
				((x*BinSize+Sstart) >= MAX00[I]) )	# Find the index of genes in the bin.
		o <- order( (MAX00[I]/MIN00[I])[i], decreasing=TRUE )	# Order this bin by MAX/MIN.
		return(i[o[c(1:NUMOfGenesInBin)]])			# Get indicies of these.
	})

	## Now, fit line segments in each group.  Plot each line segment.
	# TODO: Make each line continue until it intersects the next line in the sequence.
	LineInfo <- as.matrix(data.frame(b= rep(0,NGroups),m= rep(0,NGroups)))
	for(i in 1:(NGroups-1)) {
	#	if(sum(!is.na(c(SHELL[[i]],SHELL[[i+1]]))) == NUMOfGenesInBin*2) { 
		line <- lm( formula=MIN00[I][ c(SHELL[[i]],SHELL[[i+1]]) ] ~ 
					MAX00[I][ c(SHELL[[i]],SHELL[[i+1]]) ] )
	#	}
		LineInfo[i,] <- c(line$coefficients[[1]], line$coefficients[[2]])
	}

	# Plots line info -- for function debugging.
	plot(MAX00[I], MIN00[I])
	points(MAX00[I][unlist(SHELL)], MIN00[I][unlist(SHELL)], col="RED", cex=1, lwd=2)
	abline(0, 1, col="blue", lwd=1)
	abline(log(MinTH,10),1,col="red",lwd=1)

	### Calculate the verticies & draw lines.
	### TODO: Detect & Prevent backtrack!

	# VertexInfo: x_left, y_left, x_right, y_right.
	VertexInfo <- as.matrix(data.frame(x_left= rep(0,NGroups), y_left=rep(0,NGroups),
						x_right=rep(0,NGroups), y_right=rep(0,NGroups)))
	P1 <- c( Sstart,  Sstart*LineInfo[1,2]+LineInfo[1,1] ) # Starting point...
	for(i in 2:(NGroups-1)) {
		A <- CalcIntersection(LineInfo[i,1], LineInfo[i,2], LineInfo[i-1,1], LineInfo[i-1,2])
		P2 <- c(A$x, A$y)
		VertexInfo[i-1,] <- c(P1[1],P1[2],P2[1],P2[2])
	#	print(paste("P1:",P1[1],P1[2],"P2:",P2[1],P2[2]))


		## Detect backtrack.
	#	if(P1[1] < P2[1]) {
			lines(c(P1[1],P2[1]), c(P1[2], P2[2]), col="blue", lwd=2)
			P1 <- P2
	#	}
	}

	### Assume convergence with MinTH?!
	# Calculate Last Points.
	A <- CalcIntersection(LineInfo[NGroups-1,1], LineInfo[NGroups-1,2], log(MinTH,10), 1)
	P2 <- c(A$x, A$y)

	# Draw Lines.
	lines(c(P1[1],P2[1]), c(P1[2], P2[2]), col="blue", lwd=2)
	lines(c(P2[1], 0), c(P2[2], log(MinTH,10)), col="blue", lwd=2)

	# Save information.
	VertexInfo[NGroups-1,] <- c(P1[1],P1[2],P2[1],P2[2])
	VertexInfo[NGroups,] <- c(P2[1],P2[2],0,log(MinTH,10))
	LineInfo[NGroups,] <- c(log(MinTH,10),1)

	MODEL <- cbind(LineInfo,VertexInfo)

	### Translate the model up & down to fit your desired confidence interval.
	# Measure the number of points that are "significant".
	SignificantPoints <- PointsAboveModel(MODEL,MAX00,MIN00)
	print(paste( sum(SignificantPoints),NROW(MAX00), sum(SignificantPoints)/NROW(MAX00) ))


	# Count the number of genes, g, that need to be added or removed to obtain CI.
	# Measure the residuals of each point to the model -- above or below.
	# Re-align the model such that it passes directly through point g.

	### Calculate the mirror-inversion of the model about the line (0,1)

	# Return the model as a data structure.
	return(MODEL)
}
#MODEL <- ModelFit(NH00, LC00)


###########################################################
###
###	PlotModel -- Plots the model... 
###
###########################################################

PlotModel <- function(M) {
	for(i in 1:NROW(M)) {
		lines(c(M[i,3], M[i,5]), c(M[i,4], M[i,6]), col="green", lwd=3)
		lines(c(M[i,4], M[i,6]), c(M[i,3], M[i,5]), col="green", lwd=3)
	}
}

###########################################################
###
###	CallChangedGenes -- Calls Changed Genes using two biological
###		replicates (R1, R2) across two biological conditions (C1, C2).
###
###	Since values are all calculated on a log10-scale, genes with 0 reads are set
###	to min.
###
###########################################################

CallChangedGenes <- function(C1R1rs, Nreads1, 
				C1R2rs, Nreads2, 
				C2R1rs, Nreads3, 
				C2R2rs, Nreads4, G, FILENAME=NULL, min=0.5) {
	print("Fitting Model")
## Change genes with 0 reads to min
	C1R1rs[which(C1R1rs == 0)] <- min
	C1R2rs[which(C1R2rs == 0)] <- min
	C2R1rs[which(C2R1rs == 0)] <- min
	C2R2rs[which(C2R2rs == 0)] <- min

	C1R1 <- log(C1R1rs/Nreads1,10)
	C1R2 <- log(C1R2rs/Nreads2,10)
	C2R1 <- log(C2R1rs/Nreads3,10)
	C2R2 <- log(C2R2rs/Nreads4,10)
	M <- ModelFit(c(C1R1,C2R1), c(C1R2, C2R2))

	# Put the two replicates together.
	META00 <- log( (C1R1rs+C1R2rs)/(Nreads1+Nreads2), 10 )
	META40 <- log( (C2R1rs+C2R2rs)/(Nreads3+Nreads4), 10 )

#	META00 <- unlist(lapply(c(1:NROW(C1R1)), function(x) {mean(C1R1[x], C1R2[x])}))
#	META40 <- unlist(lapply(c(1:NROW(C1R1)), function(x) {mean(C2R1[x], C2R2[x])}))

	# Plot the meta-model
	if(!is.null(FILENAME)) {
		png(FILENAME, height=1500, width=750)
	}

	## Plot the first 
	par(mfrow=c(2,1))
	plot(C1R1, C1R2, xlab="Log(Reads In Gene/ Total) -- Replicate 1", ylab="Log(Reads In Gene/ Total) -- Replicate 2")
	points(C2R1, C2R2)
	PlotModel(M)
	abline(0,1,col="red")
        UP <- PointsAboveModel(M,C1R1,C1R2)
        points(C1R1[UP],C1R2[UP], col="red", lwd=2)
        UP <- PointsAboveModel(M,C2R1,C2R2)
        points(C2R1[UP],C2R2[UP], col="red", lwd=2)
        DN <- PointsAboveModel(M,C1R2,C1R1)
        points(C1R1[DN],C1R2[DN], col="red", lwd=2)
        DN <- PointsAboveModel(M,C2R2,C2R1)
        points(C2R1[DN],C2R2[DN], col="red", lwd=2)


	plot(META00,META40, xlab="Log(Reads In Gene/ Total) -- 0 Min", ylab="Log(Reads In Gene/ Total) -- 40 Min")
	PlotModel(M)
	abline(0,1,col="red")

	# Find points above model...
	UP <- PointsAboveModel(M,META00,META40)
	points(META00[UP],META40[UP], col="red", lwd=2)
	DN <- PointsAboveModel(M,META40,META00)
	points(META00[DN],META40[DN], col="red", lwd=2)

	print(paste("UP:",sum(UP),"DN:",sum(DN)))
	print(paste("UP:",sum(UP)/sum(!is.na(META40)),"DN:",sum(DN)/sum(!is.na(META40))))

	# Upregulated > 1, 40/0
	# Find the symbol.
	CNG <- (UP | DN)
	ChangedGenes <- data.frame(G[CNG,], FC=(10^META40[CNG])/(10^META00[CNG]), 
		C1=(C1R1rs[CNG]+C1R2rs[CNG]), C2=(C2R1rs[CNG]+C2R2rs[CNG]))

	if(!is.null(FILENAME)) {
		dev.off()
	}

	return(ChangedGenes)
}

###########################################################
###
###	LimitTo40kb -- Calls Changed Genes using two biological
###		replicates (R1, R2) across two biological conditions (C1, C2).
###
###########################################################
LimitToXkb <- function(Gbed, SIZE=40000) {
#	SIZE <- 40000

	GSt <- unlist(lapply(c(1:NROW(Gbed)), function(x) {
		if(Gbed[[4]][x] == "+") {
			return (Gbed[[2]][x] + 1000)
		}
		if(Gbed[[4]][x] == "-") {
			if( (Gbed[[3]][x] - Gbed[[2]][x]) < SIZE ) {
				return( Gbed[[2]][x] )
			}
			else {
				return( Gbed[[3]][x] - SIZE )
			}
		}
	}))

	GEn <- unlist(lapply(c(1:NROW(Gbed)), function(x) {
		if(Gbed[[4]][x] == "+") {
			if( (Gbed[[3]][x] - Gbed[[2]][x]) < SIZE ) {
				return( Gbed[[3]][x] )
			}
			else {
				return( Gbed[[2]][x] + SIZE )
			}		
		}
		if(Gbed[[4]][x] == "-") {
			return(Gbed[[3]][x] - 1000)
		}
	}))

	G <- data.frame(Chr= Gbed[[1]], Start= GSt, End= GEn, Str= Gbed[[4]], ID= Gbed[[5]])
#	save.image("ImageG40KB.RData")
	return(G)
}

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

## 0-40 minute analysis.
#require(GROseq)
#print("Reducing to +1kb to +40kb")
#load("Human.RefSeqGenes.Mar2006.BED.RData")
#G <- LimitToXkb(Gbed)

#print("Loading data")
#load("LC.NH.0M.40M.Rdata")
## Count reads in +1 to +40kb.
#RefSeqNH00 <- CountReadsInInterval(f= G, p= NH00[,c(1:3,6)])
#RefSeqNH40 <- CountReadsInInterval(f= G, p= NH40[,c(1:3,6)])
#RefSeqLC00 <- CountReadsInInterval(f= G, p= LC00[,c(1:3,6)])
#RefSeqLC40 <- CountReadsInInterval(f= G, p= LC40[,c(1:3,6)])
#nrNH00 <- NROW(NH00)
#nrNH40 <- NROW(NH40)
#nrLC00 <- NROW(LC00)
#nrLC40 <- NROW(LC40)
#remove(LC00)
#remove(LC40)
#remove(NH00)
#remove(NH40)
#save.image(file="NReads-40kb.RData")

#load(file="NReads-40kb.RData")
#ChangedGenes <- CallChangedGenes(RefSeqNH00,nrNH00,
#					RefSeqLC00,nrLC00, 
#					RefSeqNH40,nrNH40, 
#					RefSeqLC40,nrLC40, G=Gbed)
#write.table(ChangedGenes, row.names=FALSE, quote=FALSE, file= "E2.Changed_1-40kb.RefSeq.tsv", sep="\t")

## 0-10 minute analysis.
#G <- LimitToXkb(Gbed, SIZE=10000)

## 0-160 minute analysis.

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

### Calculates the fraction of reads in a refseq genes.

#require(GROseq)
#print("Reducing to +1kb to +40kb")
#load("Human.RefSeqGenes.Mar2006.BED.RData")
#G <- LimitToXkb(Gbed)
#
#print("Loading data")
#load("LC.NH.0M.40M.Rdata")
#### Count reads in +1 to +40kb.
#RefSeqNH00 <- CountReadsInInterval(f= NH00[,c(1:3,6)], p= Gbed)
#RefSeqNH40 <- CountReadsInInterval(f= NH40[,c(1:3,6)], p= Gbed)
#RefSeqLC00 <- CountReadsInInterval(f= LC00[,c(1:3,6)], p= Gbed)
#RefSeqLC40 <- CountReadsInInterval(f= LC40[,c(1:3,6)], p= Gbed)

###########################################################
###########################################################
###########################################################
###########################################################

