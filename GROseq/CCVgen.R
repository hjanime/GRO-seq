######################################################################################################################
#
#    R script.  
#
#	Which/ Coorelation Matrix script ... makes heavy use of Pearson coorelation in microarray analysis and visulization.
#	Requires: WhichGenes.R, ExpressionColorScatterPlot.R
#	Suggests: GOfunctions.R 
#
#	Functions:
#	(1) OrderCluster(MATRIX, LMT_Sym=NULL) <- Returns (sorted) coorelation coeficients for each element in MATRIX.
#	    Returns the following rownames: SYMBOL, AFFX_ID, TITLE, COR.  
#	(2) CompleteDistance(CLU, INDX_A, INDX_B, INDX_L, INDX_R) <- Returns complete-link-style D(A,L), D(B,L), D(A,R), D(B,R)
#	(3) MeanColorCoorelationMatrix(MATRIX, PATTERN_VECTOR, LMT_SYM) <- Builds a color-matrix of expression data.
#	(4) ColorCoorelationMatrix(MATRIX, LMT_Sym) <- Builds a color-matrix of expression data.
#	(5) AddGOSquares(CUTREE, ) <- sdf; ASSUMES a graphics window that it has access to.
#
#
#	Todo:
#	(X) Write Which script.
#	(X) Add support for finding any gene via which.
#	(X) Make clustring based on the gene with MAX expression in GA.
#
#
######################################################################################################################
##	write.table(ANS, "CLUSTER.101.R2_before.Information")
##	write.table(data.frame(cutree(CLU, K)[CLU$order], AFFXID[CLU$order], SYMBOL[CLU$order]), "CLUSTER.101.R2_before.Deviations")

blah <- function(x, GA_EXP, A2S) {
	source("../bin/CorrelationClusterVisulization.R")
	#png("./figures/GAgt1k--K150.Nearest.png", height=25000, width=25000)
	png("./figures/gt5k-Rsquared-DIFF-NoNames-REDUCED.png", height=4000, width=4000)
	TMP <- MeanColorCorrelationMatrix(MATRIX= EB, PATTERN_VECTOR= PATTERN, LMT_SYM= x, LimitProbesToMax= GA_EXP, A2S= A2S)
	dev.off()
	return(TMP)
}


CreateCoorelationMatrix <- function(MATRIX, NAMES, k=NULL, METH="average", MODE="matrix") {
	require(cluster);
	choppedM <- MATRIX;

	# Create the coorelation!
	cc <- cor(t(choppedM))

##	Switch to R-squared before clustering
	cc <- cc*abs(cc)

#c("average", "single", "complete", "ward", "weighted", "flexible")
	CLU <- agnes(1-cc, diss=TRUE, method=METH);#complete");  
	CLU$data <- cc; # Does not automatically get assigned when passing a dissimalarity matrix.

	CLU <- OrderCluster(CLU)

	cc <- cc[CLU$order, CLU$order];
	choppedM<- choppedM[CLU$order,];
	NAMES <- NAMES[CLU$order];
	write.table(cbind(NAMES, rownames(choppedM)), "TMPORDER.file", row.names=FALSE, col.names=FALSE)


	# Draw the matrix.
#plot(-10, -10, axes=FALSE, xlab="", ylab="", xlim=c(0, NROW(cc)), ylim=c(0, NROW(cc)));
## Draw the proper size...
CANVASSIZE <- max(5000, NROW(cc))
plot(-10, -10, axes=FALSE, xlab="", ylab="", xlim=c(0, CANVASSIZE), ylim=c(0,CANVASSIZE))

if(MODE == "matrix") {
#	plot(-10, -10, axes=FALSE, xlab="", ylab="", xlim=c(0, NROW(cc)), ylim=c(0, NROW(cc)));
	for(i in 1:NROW(cc) ) {
		for(j in 1:NCOL(cc) ) {
			if(cc[i,j]>0) {
				points(i, j, col=rgb(cc[i,j], cc[i,j], 0), pch=15, cex=3.0);
			}
			else {
				points(i, j, col=rgb(0, abs(cc[i,j]), abs(cc[i,j])), pch=15, cex=3.0);
			}
		}
	
		# Labels on X and Y axis ...
		mtext(NAMES[i], side=1, at=i, line=-1, las=2);
		mtext(NAMES[i], side=2, at=i, line=-1, las=2);
	}
}

if(MODE == "heatmap") {
	print(MODE)
	## Center/scale
	SIZE = 250
	CS <- t(scale(t(choppedM), center=choppedM[,1]))
	CS <- CS/max(abs(CS))
	for(i in 1:NROW(CS) ) {
		## Change to [0/1].
#		ZeroOne <- ( CS[i,] - min(CS[i,]) )/max(CS[i,] - min(CS[i,]))
		for(j in 1:NCOL(CS)) {
			if(CS[i,j] > 0) {
				C <- rgb(CS[i,j],CS[i,j],0)
			}
			else {
				C <- rgb(0, 0, abs(CS[i,j]))
			}

                         points(j*SIZE+(0*SIZE/5),i, col=C, pch=15, cex=5.0)
                         points(j*SIZE+(1*SIZE/5),i, col=C, pch=15, cex=5.0)
                         points(j*SIZE+(2*SIZE/5),i, col=C, pch=15, cex=5.0)
                         points(j*SIZE+(3*SIZE/5),i, col=C, pch=15, cex=5.0)
                         points(j*SIZE+(4*SIZE/5),i, col=C, pch=15, cex=5.0)
                         points(j*SIZE+(5*SIZE/5),i, col=C, pch=15, cex=5.0)
		}
		points(NCOL(CS)*SIZE+SIZE,i,col=rgb(1,1,1), pch=15, cex=3.0)
                mtext(NAMES[i], side=2, at=i, line=-1, las=2);
	}
}

	if(!is.null(k)) {
	for(k in c(4:20)) {
		PAR <- AddGOSquares(CLU, K=k, include.labels=FALSE);
		write.table(PAR, paste("CLU/RES.",k, sep=""));
		write.table(data.frame(cutree(CLU, k=k)[CLU$order], NAMES, rownames(MATRIX)[CLU$order]), file=paste("CLUSTER.",k,".Deviations", sep=""), sep="\t", quote=FALSE)

		## Now, create plots.
		TIME <- c(0, 10, 40, 160)
		CLUSTER <- cutree(CLU, k=k)[CLU$order]
		for(num in c(1:k)) {
			print(paste("Cluster:",num,"of",k,"Size=",NROW(which(CLUSTER==num)) ))
			png(paste("IMG/k_",k,".",num,".CLUSTER.png",sep=""), height=750, width=750)
				par(font=3, font.lab=2, font.axis=2, mgp=c(5.5, 2, 0), mar=c(8.5, 8.5, 2, 2) + 0.3)
				MIN <- min(CS[which(CLUSTER==num),])
                                MAX <- max(CS[which(CLUSTER==num),])
				MEANS <- colMeans(CS[which(CLUSTER==num),])
				plot(TIME, MEANS, col="black", cex.axis=3.5, cex.lab=4, 
					type="l", xlab="E2 (min)", ylab="Relative Expression")#, ylim=c(MIN,0))
				for(row in which(CLUSTER==num))
					points(TIME, as.real(CS[row,]), col="light gray", type="l")
				points(TIME, MEANS, col="blue", type="l", lwd=6)
			dev.off()
		}
	}
	}

	return(CLU);
}
####################################################################################################
##	AddGOSquares(CLU, K, include.labels) <- Adds lines around selected positions; 
##						ASSUMES a graphics window is open & ready!
##
##	July 1, 2008: Added functions to find optimum level of chopping dendrogram.
##		Now reports the folowing paremeters:
##			Let Ng be the median number of genes per cluster.
##			Let Ybf be the ratio of genes that share the most significant biological
##			  function to those that are of a different (or unknown) function.
##			Let Zhy be the median p-value (of biological function enrichment) for all 
##			  clusters by the hypergeometric distribution.
##			Let C be the mean (median? product?) of mean (or min?) intra-cluster 
##			  correlation.
##
####################################################################################################
AddGOSquares <- function(CLU, K, include.labels=FALSE) {
	CT <- cutree(CLU, k=K);

	# Does not account for overlap between probes!
#	TMP <- hist(CT, K);

	N <- NULL;
	Y <- NULL;
	C <- NULL;
	P <- NULL;
	NAMES <- NULL;

	CT <- CT[CLU$order];
	if(include.labels) {
		GC <- GOCluster(CLU, CLU$order.lab, K=K, ontology="BP");
#		GC <- GOCluster(CLU, rownames(GA), K=K, ontology="BP");
	}
	for(i in 1:K) {
		# Find the index of the first & last.
		IF <- min(which(CT == i))-0.5;
		IL <- max(which(CT == i))+0.5;
		#print(paste(K, i, IF, IL));

		# Draw ... Square should go from: (INDXFirst-1, INDXFirst-1) -> (INDXFirst-1,INDXLast+1)
		## Removed this line -- only used in the correlation matrix visualization
#		segments(IF, IF, IF, IL, col=rgb(1, 1, 1), lwd=10)#, CEX=3.0);
#		segments(IF, IL, IL, IL, col=rgb(1, 1, 1), lwd=10)#, CEX=3.0);
#		segments(IL, IL, IL, IF, col=rgb(1, 1, 1), lwd=10)#, CEX=3.0);
#		segments(IL, IF, IF, IF, col=rgb(1, 1, 1), lwd=10)#, CEX=3.0);

		if(include.labels) {
			# Get the stats for the pathway ...
			if(!is.null(GC[[i*2]]) && (length(GC[[i*2-1]]) > 1)) {

				INDX <- which( max(GC[[i*2]][[4]]) ==  GC[[i*2]][[4]]);
				NInPathway <- GC[[i*2]][[4]][1];#[INDX];
				NTotal     <- NROW(GC[[i*2 - 1]]);
				PathName   <- GC[[i*2]][[2]][1];#[INDX];
				Pval       <- GC[[i*2]][[6]][1];#[INDX];
			}
			else {
				INDX <- 1;
				NInPathway <- 1;
				NTotal <- 1;
				Pval <- 1;
				PathName <- "";
			}

			N[i] <- NTotal;
			Y[i] <- NInPathway[1];
			C[i] <- as.real( IntraClusterDistanceMean(CLU,CLU$order[which(CT == i)]) );
			P[i] <- Pval;
			NAMES[i] <- PathName[1];

			if(((IF+IL)/2) < (NROW(CT)/2)) {
				ADJ=0;
				END=IL+1;
			}
			else {
				ADJ=1;
				END=IF-1;
			}

			# Write text to screen.
			text(END, (IF+IL)/2, 
				paste(PathName[1], " ", NInPathway[1], "/ ", NTotal, "\n", 
#				      PathName[2], " ", NInPathway[2], "/ ", NTotal, "\n",
#				      PathName[3], " ", NInPathway[3], "/ ", NTotal,  
					sep=""), 
				col=rgb(1, 1, 1), adj=ADJ, cex=5);
		}
		else {
			N[i] <- NROW(which(CT == i));
			Y[i] <- 0;
			C[i] <- as.real( IntraClusterDistanceMean(CLU,CLU$order[which(CT == i)]) );
			P[i] <- 1;
			NAMES[i] <- "NA";
	
		}
	}


	return(data.frame(1:K, N, Y, C, P, NAMES));
}

###############################################################################
##
##	Implementation of our own algorithm making 
##
##	MATRIX <- a matrix of coorelation coeficients.
##	Initial_Order <- index in order (think: CLU$order).
##
###############################################################################
OrderMatrix <- function(MATRIX, Initial_Order) {
	
}

## Calculates the optimal order across a cluster.  Returns modified CLU$order.
OrderCluster <- function(CLU, DEBUG = FALSE) {

	## Work from the top of the tree to the bottom.
	## use CLUSTER$height to to sort into groups.
	PREVIOUS_MAX <- max(CLU$height) + 1;
	INDX <- c(0,NROW(CLU$data));
	while(NROW(INDX) < NROW(CLU$height)) {
		# Gets the index of the next.
		M_CURRENT <- which(CLU$height == max(CLU$height[which(CLU$height < PREVIOUS_MAX)]));
		if(DEBUG) {
			print(paste("M_CURRENT:", M_CURRENT));
			print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");
		}

		for(LongComplexVariable in 1:NROW(M_CURRENT)) {
			CURRENT <- M_CURRENT[LongComplexVariable];

			INDX <- c(INDX, CURRENT);
			INDX <- sort(INDX);

			if(DEBUG) {
				print(paste("CURRENT:", CURRENT));
				if(NROW(CURRENT) > 1) {
					print(paste("WARNING: CURRENT has", NROW(CURRENT), "elements.  Using CURRENT =", CURRENT[1]));
					CURRENT <- CURRENT[1];
				}
			#	print(paste("INDX:", INDX));
			}

			## Get the new clusters below the split. 
			## If CLU$height referes to the space between clusters, and  then:
			## Note that the right side of couster is always +1

			## Also do error checking: Is there a left-more, or is A the first?
			if( which(INDX == CURRENT) > 1 ) { # Least it should be is 1 (-1+1) == 1, so no problem! However, problem during i==1.
				Group_A_INDX = (INDX[ which(INDX == CURRENT) -1 ]+1) : (CURRENT);
			}
			else {
				Group_A_INDX = ( 1 )                                 : (CURRENT);
			}

			if(DEBUG) {
				print(paste("A:", Group_A_INDX));
			}

			## Also do error checking: Is there a right-more?
			if( which(INDX == CURRENT) < NROW(INDX) ) {
				Group_B_INDX = (CURRENT + 1)                         : (INDX[which(INDX == CURRENT)+1]);
			}
			else {
				Group_B_INDX = (CURRENT + 1)                         : (NROW(CLU$data));
			}

			if(DEBUG) {
				print(paste("B:", Group_B_INDX));
			}		

			## Find the groups (L)eft and (R)ight of A and B, called, creatively, L and R.
			## If it is not on the bounds...
			if(Group_A_INDX[1] > 1) { # If the first element has INDX greater than 1
#				print( paste("Start:", (INDX[which(INDX == CURRENT) - 2] +1), "END:", Group_A_INDX[1] -1 ));
				START_VARIABLE_REALLY <- (INDX[which(INDX == CURRENT) - 2] +1);
				Group_L_INDX = (START_VARIABLE_REALLY +1) : Group_A_INDX[1] -1;
				
				## Odd bug!!  For some reason, the previous line gives n(expected)-1.  Bug in R??  Don't get it!
				if(Group_L_INDX[1] != START_VARIABLE_REALLY) {
					print("ERROR: Chosen wrong start variable");
					return(NULL);
				}
#				print(paste("L:", Group_L_INDX));
			}
			else {
				Group_L_INDX = NULL;
			}

			if(Group_B_INDX[NROW(Group_B_INDX)] < NROW(CLU$data)) { # If the last element has indx less than the last element.
				Group_R_INDX = (INDX[ which(INDX == CURRENT) +1] +1) : (INDX[ which(INDX == CURRENT) +2]);
			}
			else {
				Group_R_INDX = NULL;
			}
	
			if(DEBUG) {
				print(paste("L:", Group_L_INDX));
				print(paste("R:", Group_R_INDX));
			}

			## Calculate DISTANCE between both new branches & both placed branches on either side.
#			DIST_MATRIX <- CompleteDistance(CLU, Group_A_INDX, Group_B_INDX, Group_L_INDX, Group_R_INDX);
			DIST_MATRIX <- AverageDistance(CLU, Group_A_INDX, Group_B_INDX, Group_L_INDX, Group_R_INDX);
#			DIST_MATRIX <- NearestDistance(CLU, Group_A_INDX, Group_B_INDX, Group_L_INDX, Group_R_INDX);
			if(DEBUG) {
				#print(paste("DIST_MATRIX:",DIST_MATRIX));
				print(paste("Dal:",DIST_MATRIX$Dal,"Dbr:",DIST_MATRIX$Dbr,"Dar:",DIST_MATRIX$Dar,"Dbl:",DIST_MATRIX$Dbl));
			}

			## Take order that satisfies: MIN{ [D(A,L)+D(B,R)], [D(A,R)+D(B,L)] }
			## A starts on left, B on right
			## Therefore, if A shoudl be near L, then no switch is necessary.
			## Otherwise, if ... then make a switch ... otherwise, no action necessary.
			if( (DIST_MATRIX$Dal + DIST_MATRIX$Dbr) < (DIST_MATRIX$Dar + DIST_MATRIX$Dbl) ) {
	####### TODO: Write this function!!!!
				if(DEBUG) {
					print("Switching");
				}

				## Caluclate new order.
				if(Group_A_INDX[1] == 1) { ## If group A starts on the far right.
					delta_ORDER <- c(Group_B_INDX, 
								Group_A_INDX, 
								(Group_B_INDX[NROW(Group_B_INDX) ]+1):NROW(CLU$data));
					height_ORDER <- c(Group_B_INDX[0:(NROW(Group_B_INDX) -1)], 				#B heights, nothing left
								Group_A_INDX[NROW(Group_A_INDX)], 				#CENTER
								Group_A_INDX[0:(NROW(Group_A_INDX) -1)],			#A heights
								(Group_B_INDX[NROW(Group_B_INDX)]):(NROW(CLU$data)-1) );	#RIGHT of all
				}
				else if(Group_B_INDX[NROW(Group_B_INDX)] == NROW(CLU$data)) { ## If at the end!
					delta_ORDER <- c(1:(Group_A_INDX[1]-1), 
								Group_B_INDX, 
								Group_A_INDX);
					height_ORDER <- c(1:(Group_A_INDX[1]-1), 						#LEFT of all
								Group_B_INDX[0:(NROW(Group_B_INDX) -1)], 			#B heights
								Group_A_INDX[NROW(Group_A_INDX)], 				#CENTER
								Group_A_INDX[0:(NROW(Group_A_INDX) -1)] );			#A heights, nothing right.
				}
				else {
					delta_ORDER <- c(1:(Group_A_INDX[1]-1), 
								Group_B_INDX, 
								Group_A_INDX, 
								(Group_B_INDX[NROW(Group_B_INDX)]+1):NROW(CLU$data) );
						height_ORDER <- c(1:(Group_A_INDX[1]-1), 						#LEFT of all
									Group_B_INDX[0:(NROW(Group_B_INDX) -1)], 			#B heights
									Group_A_INDX[NROW(Group_A_INDX)], 				#CENTER
									Group_A_INDX[0:(NROW(Group_A_INDX) -1)],			#A heights
									(Group_B_INDX[NROW(Group_B_INDX)]):(NROW(CLU$data)-1) );	#RIGHT of all
#					print(paste( Group_A_INDX[1:(NROW(Group_A_INDX) -1)], Group_B_INDX[1:(NROW(Group_B_INDX) -1)], sep=" " ))
				}

				## Preform the switch inside the current CLU object, including: $data, $height, $order
				CLU$order  <- CLU$order[delta_ORDER];
				CLU$order.lab <- CLU$order.lab[delta_ORDER];
				CLU$SYM <- CLU$SYM[delta_ORDER];

# Heights naturally stay in the same order.  
#TODO: WRONG!  Only works out like that if the rows being switched are symemtric.
				CLU$height <- CLU$height[height_ORDER];
			
				## INDX represents the index of the current split ... 
				## This will change if the splitl is symmetric.
				## The new value should be: which(height_ORDER == Group_A_INDX[NROW(Group_A_INDX)]), units: INDEX
				INDX[which(INDX == CURRENT)] <- which(height_ORDER == Group_A_INDX[NROW(Group_A_INDX)]);

				#CLU$data   <- CLU$data[ORDER, ORDER]; # $data is NOT sorted to begin with!

				if(DEBUG) {
					print(paste("Aindx:",Group_A_INDX));
					print(paste("Bindx:",Group_B_INDX));
					print("delta_ORDER:");
					print(delta_ORDER);
					print("height_ORDER:");
					print(height_ORDER);
#					print("CLU$order:");
#					print(CLU$order);
				}
			}

			if(DEBUG) {
				print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");
			}
		}
		# Prepare for the next round!
		PREVIOUS_MAX <- CLU$height[ which( CLU$height == max(CLU$height[which(CLU$height < PREVIOUS_MAX)]) )[1] ];
		if(DEBUG) {
			print(paste("Previous Max:",PREVIOUS_MAX));
			print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");

		}


	}
	return(CLU);
}

############################################################
##
## Calculates distances between clusters in the MATRIX.
##
## Assumes that CLU$data is a matrix of Pearson correlations.
##
## CompleteDistance <- Max difference between single components taken to be distance.
## Since PC, min takes the complete link clustering ... if we assume a 
## difference matrix, we would instead take max.
##
## AverageDistance <- Average difference between all components taken to be distance.
##
############################################################
CompleteDistance <- function(CLU, INDX_A, INDX_B, INDX_L, INDX_R) {
	## Calculate: D(A, L) = max(D(member(A), member(L))), D(B,L), D(A,R), and D(B,R).
	if ( !is.null(INDX_L) ) {
		Dal <- min(CLU$data[CLU$order[INDX_A],CLU$order[INDX_L]]);
		Dbl <- min(CLU$data[CLU$order[INDX_B],CLU$order[INDX_L]]);
	}
	else {
		Dal <- 0;
		Dbl <- 0;
	}

	if ( !is.null(INDX_R) ) {
		Dar <- min(CLU$data[CLU$order[INDX_A],CLU$order[INDX_R]]);
		Dbr <- min(CLU$data[CLU$order[INDX_B],CLU$order[INDX_R]]);
	}
	else {
		Dar <- 0;
		Dbr <- 0;
	}

	return(data.frame(Dal = Dal, Dbl = Dbl, Dar = Dar, Dbr = Dbr));
}

AverageDistance <- function(CLU, INDX_A, INDX_B, INDX_L, INDX_R) {
	## Calculate: D(A, L) = max(D(member(A), member(L))), D(B,L), D(A,R), and D(B,R).
	if ( !is.null(INDX_L) ) {
		Dal <- mean(CLU$data[CLU$order[INDX_A],CLU$order[INDX_L]]);
		Dbl <- mean(CLU$data[CLU$order[INDX_B],CLU$order[INDX_L]]);
	}
	else {
		Dal <- 0;
		Dbl <- 0;
	}

	if ( !is.null(INDX_R) ) {
		Dar <- mean(CLU$data[CLU$order[INDX_A], CLU$order[INDX_R]]);
		Dbr <- mean(CLU$data[CLU$order[INDX_B], CLU$order[INDX_R]]);
	}
	else {
		Dar <- 0;
		Dbr <- 0;
	}

	return(data.frame(Dal = Dal, Dbl = Dbl, Dar = Dar, Dbr = Dbr));
}

NearestDistance <- function(CLU, INDX_A, INDX_B, INDX_L, INDX_R) {
	## Calculate: D(A, L) = max(D(member(A), member(L))), D(B,L), D(A,R), and D(B,R).
	if ( !is.null(INDX_L) ) {
		Dal <- max(CLU$data[CLU$order[INDX_A],CLU$order[INDX_L]]);
		Dbl <- max(CLU$data[CLU$order[INDX_B],CLU$order[INDX_L]]);
	}
	else {
		Dal <- 0;
		Dbl <- 0;
	}

	if ( !is.null(INDX_R) ) {
		Dar <- max(CLU$data[CLU$order[INDX_A], CLU$order[INDX_R]]);
		Dbr <- max(CLU$data[CLU$order[INDX_B], CLU$order[INDX_R]]);
	}
	else {
		Dar <- 0;
		Dbr <- 0;
	}

	return(data.frame(Dal = Dal, Dbl = Dbl, Dar = Dar, Dbr = Dbr));
}

############################################################
##
## Calculates intra-cluster variation.
##
## Assumes that CLU$data is a matrix of Pearson correlations.
##
## This is run on CUTREE, which means INDX references the order in
##   CLU$data!  This is DIFFERENT from above!
##
############################################################
IntraClusterDistanceMedian <- function(CLU, cutreeINDX) {
	TMP <- as.vector(CLU$data[cutreeINDX, cutreeINDX])
	if( !is.null(cutreeINDX) ) {
		D <- median(TMP[TMP <= 0.99]);
	}

	return(D);
}
IntraClusterDistanceMean <- function(CLU, cutreeINDX) {
	TMP <- as.vector(CLU$data[cutreeINDX, cutreeINDX])
	if( !is.null(cutreeINDX) ) {
		D <- mean(TMP[TMP <= 0.99])
	}

	return(D);
}

