### DataS7_GeneticMappingAnalyses.R
## All code for plotting figures in Schell, et al. 2021 are contained in this script
## Code is arranged in order of presented data
## Note, code compatible, tested under R version 3.5.1 (2018-07-02) -- 'Feather Spray" on x86_64-apple-darwin15.6.0 (64-bit) platform and R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" on x86_64-apple-darwin15.6.0 (64-bit)
# Necessary files for these analyses are DataS1_Genotype.txt', 'DataS2_Phenotype.txt', 'DataS3_ReciprocalHemizygosity.txt', 'DataS4_CloningExperiments.txt', 'DataS5_PetiteFrequency.txt', 'DataS6_ChromosomeIIDuplication.txt'
# Additionally, resultant linkage mapping files generated in DataS7.GeneticMappingAnalyses.R files are used for plotting, results, etc.


#### General overview/ layout of analyses in this file
## Figures are presented in the same order as discussed in main text, followed by the order in the supplement of Schell et al. 2021 
##Include all relevant files in the working directory

library(RColorBrewer) # v 1.1.2
library(scales) # v 1.0.0
library(zoo) # 1.8.4
library(RColorBrewer) # v 1.1.2
library(scales) # v 1.0.0
library(vioplot) # v 0.3.4

### Read in general files which are necessary for all subsequent sections
#read in general geno, pheno files from which subsets of data will be take for specific plots
genoAll <- read.table('DataS1_Genotype.txt', header=TRUE, stringsAsFactors=FALSE)
phenoAll <- read.table('DataS2_Phenotype.txt', header=TRUE, stringsAsFactors=FALSE)
#Add a universal position column
cls <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066, 85779) #lengths of chromosomes from S288C_reference_sequence_R64-2-1_20150113 genome release
clSums <-  cumsum(cls)
genoAll$gp <- unlist(lapply(1:17, function(x) {
	if (x != 1) {
		genoAll$p[genoAll$c==x]+clSums[x-1]
	} else {
		genoAll$p[genoAll$c==x]
	}
	
}))
colsBlOr <- c("#377EB8", "#FF7F00")




###### Below are distinct sections of analysis. ######




##### Figure 1


# Fig. 1B Compare MRP20 to mrp20-105E segregants

#Try violin plot and merge of data to be WT vs mrp20-105E

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, stringsAsFactors=FALSE)
result$gp <- genoAll$gp
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak4Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
mut <- as.factor(c('hos3', 'WT')[as.numeric(grepl('A', pheno$Sample)) +1]) #cross A is WTs, cross B is hos3∆s
listPops <- list( pheno$Ethanol[peak4Calls == 0 & mut=='WT'], pheno$Ethanol[peak4Calls == 1 & mut=='WT'], pheno$Ethanol[peak4Calls == 0 & mut=='hos3'], pheno$Ethanol[peak4Calls == 1 & mut=='hos3'])
newListPops <- list(c(listPops[[1]], listPops[[2]], listPops[[4]]), listPops[[3]]) 
listX <- seq(1,1.25, by=0.25)
dev.new(width=3.5, height=5, unit='cm') #1/3 width of one col fig
par(mar=c(2,2,1,1), oma=c(1,1,1,1), family='sans')
plot(0,0, xlim=c(0.85,1.4), ylim=c(0,max(pheno$Ethanol)*1.1), xlab='', ylab='', xaxt='n')
cols <- c('black', 'black', 'black', 'black')
xs <- c(rep(listX[1], length(newListPops[[1]])), rep(listX[2], length(newListPops[[2]])))
ys <- unlist(newListPops)
newData <- as.data.frame(cbind(xs, ys))
colnames(newData) <- c('xVal', 'growth')
vioplot(growth ~ xVal, data=newData, xlab='', xaxt='n')
mtext('Growth relative to BY', side=2, line=2)
#axis(1, at=listX, labels=c('MRP20', 'mrp201-05E'), cex.axis=0.7)





##### Figure 2

# Fig. 2A   Linkage mapping within hos3∆ random spores detects Chr XIV genetic interaction with mrp20-105E

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('mrp20*locus_withinhosMutOnly_natcommNew.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
result$gp <- genoAll$gp
genoMap <- result
cols <- rep(c('gray', 'black'), 9)[-18]
dev.new(width=7.5, height=5, unit='cm')
par(mar=c(2,2,0,0), oma=c(1,1,1,1), family='sans')
plot(-10,-10, xlim=c(0, max(genoMap$gp)), ylim=c(0, 16), xlab='', ylab='', xaxt='n', yaxt='n')
plotChrs <- sapply(1:17, function(x) {
	###add in dummy value at end of chromosome so that it connects to the first point of the next chromosome
	xChr <- genoMap$gp[genoMap$c == x]
	yChr <- genoMap$mrp20LocusInt[genoMap$c == x]
	if (x == 1) {
		prevX <- 0
	} else {
		prevX <-  max(genoMap$gp[genoMap$c == x-1])
	}
	if (x == 17) {
		nextX <- 0
	} else {
		nextX <- min(genoMap$gp[genoMap$c == x+1])
	}	
	prevY <- 1
	nextY <- 1
	xChr <- c(prevX, xChr, nextX)
	yChr<- c(prevY, yChr, nextY)
	points(x=xChr, y=-log10(yChr), pch=19, col=cols[x], type='l')
})
xLabVals <- c()
for (x in (1:17)) {
	xChr <- genoMap$gp[genoMap$c == x]
	xLabVals <- c(xLabVals, median(xChr))
}
axis(1, at=xLabVals, labels=as.roman(1:17), tick=TRUE, cex.axis=0.4)
axis(2, at=c(seq(0, 16, by=2)), labels=c(as.character(seq(0, 16, by=2))), tick=TRUE, cex.axis=0.5)
mtext('-log10 (pval) ', side=2, line=2, cex=0.7)
mtext('Genomic position', side=1, line=2, cex=0.7)


# Fig. 2B   Plot the growth effect of Chr XIV locus 

result <- read.table('mrp20*locus_withinhosMutOnly_natcommNew.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
#Subset out wt and hos3∆ segregants
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
#Now, subset only hos3∆s
phenoHos <- pheno[grepl('B', pheno$Sample),]
genoHos <- geno[, colnames(geno) %in% c('c', 'p') | colnames(geno) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
# 0 is BY allele which contains mrp20-105E mutation in hos3∆ segregants
peak14 <- result[which(result[,3] == min(result[,3], na.rm=TRUE)),]
peak14Pval <- min(result[,3], na.rm=TRUE) # 4.330511e-16
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak14Calls <- as.factor(as.numeric(genoHos[genoHos$c == peakMid$c & genoHos$p == peakMid$p, colnames(genoHos) %in% phenoHos$Sample]))
#plot 
listPops <- list( phenoHos$Ethanol[peak14Calls == 0 & mrpHos == 1], phenoHos$Ethanol[peak14Calls == 1 & mrpHos == 1], phenoHos$Ethanol[peak14Calls == 0 & mrpHos == 0], phenoHos$Ethanol[peak14Calls == 1 & mrpHos == 0])
listX <- seq(1,1.75, by=0.25)
dev.new(width=7, height=5, unit='cm') #1/3 width of one col fig
par(mar=c(2,2,1,1), oma=c(1,1,1,1), family='sans')
plot(0,0, xlim=c(0.85,1.9), ylim=c(0,max(pheno$Ethanol)*1.1), xlab='', ylab='', xaxt='n')
cols <- c('black', 'black', 'black', 'black')
plotPops <- lapply(1:4, function(x) {
	yvals <- listPops[[x]]
	xVals <- jitter(rep(listX[x], length(yvals)), amount=0.075)
	points(x=xVals, y=yvals, pch=19, cex=0.7, col=alpha(cols[x], alpha=0.2))
})
mtext('Growth relative to BY', side=2, line=2)
axis(1, at=listX, labels=c('BY', '3S', 'BY', '3S'), cex.axis=0.7)
mtext('Genotype at Chromsome XIV Locus', side=1, line=2)
points(x=c(listX[1], listX[2]), y=c(1.45, 1.45), type='l')
text(median(listX[1:2]), y=1.5, expression(italic('MRP20')))
points(x=c(listX[3], listX[4]), y=c(1.45, 1.45), type='l')
text(median(listX[3:4]), y=1.5, expression(italic('mrp20-105E')))


# Fig. 2D   F3_C cross resolves Chr XIV peak; linkage mapping significancne on chr 14
pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_intercrossLocusPvalNew.txt', header=TRUE, stringsAsFactors=FALSE)
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak14Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
result$gp <- genoAll$gp
genoPlot <- result[result$c == 14,]
genoPlot$locusPval[is.na(genoPlot$locusPval)] <- 1
cols <- rep(c('gray', 'black'), 9)[-18]
dev.new(width=4, height=5)
par(mar=c(2,2,0,0), oma=c(1,1,1,1), family='sans')
plot(-10,-10, xlim=c(min(genoPlot$p), max(genoPlot$p)), ylim=c(0, 45), xlab='', ylab='', xaxt='n', yaxt='n')
plotChrs <- sapply(14, function(x) {
	xChr <- genoPlot$p[genoPlot$c == x]
	yChr <- genoPlot$locusPval[genoPlot$c == x]
	points(x=xChr, y=-log10(yChr), pch=19, col=cols[x], type='l')
})
xLabVals <- c()
for (x in (1:17)) {
	xChr <- genoPlot$gp[genoPlot$c == x]
	xLabVals <- c(xLabVals, median(xChr))
}
axis(2, at=c(seq(0, 45, by=5)), labels=c(as.character(seq(0, 45, by=5))), tick=TRUE, cex.axis=0.5)
axis(1, at=seq(from=1,to=max(genoPlot$p), by=100000), labels=FALSE, tick=TRUE, cex.axis=0.4)
axis(1, at=c(467219, max(genoPlot$p)), labels=c('467,219', ''), tick=TRUE, cex.axis=0.5)
mtext('Position on Chr XIV', side=1, line=2)
mtext('-log10(pval)', side=2, line=2)


# Fig. 2E   Recombination breakpoints delimit Chr XIV peak in F3_C cross to a single SNP

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_intercrossLocusPvalNew.txt', header=TRUE, stringsAsFactors=FALSE)
peak14 <- result[which(result[,3] == min(result[,3], na.rm=TRUE)),]
peak14Pval <- min(result[,3], na.rm=TRUE) # 2.503996e-43
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
calls <- geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]
genoMap <- geno
allBy14 <- colnames((calls[,calls == 0]))
all3s14 <- colnames((calls[,calls == 1]))
bound <- 6000
region <- genoMap[genoMap$c == 14 & genoMap$p >= 467219-bound & genoMap$p <= 467219 + bound,]
sub <- region[, colnames(region) %in% pheno$Sample] #55 rows
subBy <- sub[,colnames(sub) %in% allBy14]
sub3s <- sub[,colnames(sub) %in% all3s14]
bounds <- c(462411, 471102)
region <- geno[geno$c == 14 & geno$p >= bounds[1] & geno$p <= bounds[2],] #make physical center of 4 minimum snps so that plot is centered on minimum BY region
sub <- region[, colnames(region) %in% pheno$Sample]
subBy <- sub[,colnames(sub) %in% allBy14]
sub3s <- sub[,colnames(sub) %in% all3s14]
phenoBy <- pheno[pheno$Sample %in% colnames(subBy),]
indsBy <- phenoBy$Sample[order(phenoBy$Ethanol)] #ordered worse to better growers
indsBy <- indsBy[c(1:25,45, 26:38, 52, 39:51)] #reorder inds if breaks arent visible
pheno3s <- pheno[pheno$Sample %in% colnames(sub3s),]
inds3s <- pheno3s$Sample[order(pheno3s$Ethanol)] #ordered worse to better growers
inds3s <- inds3s[c(11:20,1:10, 21:length(inds3s))]
genoBy <- subBy[match(indsBy, colnames(subBy))]
geno3s <- sub3s[match(inds3s, colnames(sub3s))]
#concatenate groups to put worse growing class on bottom, by then 3s on top
genoCalls <- list(genoBy, geno3s) #ok, so now if you move list1 (worst class iv-by) col1 to col n worst to best grower
snpLength <- nrow(genoBy)
indLength <- c(ncol(genoCalls[[1]]) + ncol(genoCalls[[2]]) )
bp <- max(bounds)-min(bounds)
snps <- cbind(region[,1:3], genoBy, geno3s)
#fill in rows of redundant calls not on snps and convert snp table ito physical base by base call table
y <- matrix(nrow=0, ncol=ncol(snps))
for (x in c(bounds[1]:c(bounds[2]))) {
	print(paste(bounds[2]- x, 'rows left'))
	if (x %in% snps$p) { #snp row exists
		y <- rbind(y, snps[snps$p == x,])
	} else { #non snp row fill with NAs
		blank <- matrix(data=rep(NA, ncol(snps)), ncol=ncol(snps))
		blank[,1] <- 14
		blank[,2] <- x
		blank[,3] <- NA 
		colnames(blank) <- colnames(snps)
		y <- rbind(y, blank)
	}
}
fillInCalls <- function(inputVector) {
	calls <- inputVector
	test2 <- na.fill(na.approx(inputVector, na.rm=FALSE), "extend")
	test2[test2 < 0.5 & test2 != 0] <- 0
	test2[test2 >= 0.5 & test2 != 1] <- 1
	#fill in leading or trailing NAs
	test2 <- na.locf(test2)
	test2
}
allPos <- do.call('cbind', lapply(4:ncol(y), function(x) {
	newCalls <- fillInCalls(y[,x])
}))
colnames(allPos) <- colnames(y)[4:ncol(snps)]
final14 <- cbind(y[,1:2], allPos)
matrixPos <- t(t(final14[,colnames(final14) %in% pheno$Sample]))
#plot
newCol <- brewer.pal(5, 'Set1')[c(2,5)]
dev.new(width=4.6, height=3.7, unit='cm') #if too small this warps the image and removes certain lines!!!! 4.4 x 4 works and so does 5x4
par(mar=c(3,2,0.5,0.5), oma=c(0,0,0,0))
image(matrixPos, xaxt='n', yaxt='n', col=newCol)
xLabPos <- ((final14$p[final14$p %in% snps$p]) - min(final14$p)) / bp
axis(1, at=xLabPos[-13], labels=FALSE, tick=TRUE)
axis(1, at=xLabPos[13], labels='467,219', tick=TRUE, cex=0.5)
abline(v=xLabPos[snps$p == 467219], col='black', lty=3)
mtext('Haplotype at peak position', side=2, line=1)
mtext('Chr XIV: 462,411 to 471,102', side=1, line=2)



# Fig. 2F   Cloning causal SNP at MKT1 in segregants

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
segsMkt <- cloning[cloning$Type == 'segregant' & cloning$Gene == 'MKT1',]
allYs <- list( segsMkt$Ethanol[segsMkt$Edit == 'from3StoBYcandidate'], segsMkt$Ethanol[segsMkt$Edit == 'from3StoBY-1SNP'], segsMkt$Ethanol[segsMkt$Edit == 'from3StoBY+1SNP'])
dev.new(width=4, height = 5)
par(mar=c(1,1,1,1), oma=c(2,2,2,2))
plot(x=-1, y=-1, xlim=c(0.5,2.5), ylim=c(0, 1.5), xaxt='n', xlab='', ylab='')
plotVals <- sapply(c(1:length(allYs)), function(x) {
	yVals <- allYs[[x]]
	xVals <- jitter(rep(c(1,2,2)[x], length(yVals)), amount=0.2)
	points(x=xVals, y=yVals, col=alpha('black', alpha=0.4), pch=19, cex=0.7)
})
mtext('Growth relative to BY', side=2, line=2)
mtext('Genotype', side=1, line=2)
axisLabels  <- c(expression(italic('MKT1'^'3S'*'to MKT1'^'BY'*'at causal')), expression(italic('MKT1'^'3S'*'to MKT1'^'BY'*'at +/- 1 SNP')))
axis(1, side=1, at=c(1:2), labels=axisLabels, cex.axis=0.5)





# Figure 3



# Fig. 3A   Cloning causal nucleotides at MRP20 and MKT1 into BY and 3S parental strains

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
parents <- cloning[cloning$Type == 'parent',]  
#get growth vlaues for MRP20 background then for mrp20-105E in this order:
# BY MKT1-BY, BY MKT1-3S, 3S MKT1-BY, 3S MKT1-3S 
wts <- list(parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'MRP20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'MRP20' & parents$MKT1 == 1],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'MRP20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'MRP20' & parents$MKT1 == 1])
muts <- list(parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$MKT1 == 1],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$MKT1 == 1])
#Get epected value of growth from observed segregant growth values in original BYx3S F2 wt and hos3∆ segregants
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
#Subset, only wt segs
phenoWt <- pheno[grepl('A', pheno$Sample),]
genoWt <- geno[, colnames(geno) %in% c('c', 'p') | colnames(geno) %in% phenoWt$Sample]
mktWt <- as.factor(as.numeric(genoWt[genoWt$c == 14 & genoWt$p == 467219, colnames(genoWt) %in% phenoWt$Sample]) )
#Now, subset only hos3∆s
phenoHos <- pheno[grepl('B', pheno$Sample),]
genoHos <- geno[, colnames(geno) %in% c('c', 'p') | colnames(geno) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
mktHos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample]) )
#put originally observed wt and mrp20-105E segregants into a list of WT MKT1 BY, WT MKT1 3S, mrp20-105E MKT1 BY, mrp20-105E MKT 3S
segGrowth <- list ((phenoWt$Ethanol[mktWt == 0]), (phenoWt$Ethanol[mktWt == 1]), c(phenoHos$Ethanol[mrpHos == 0 & mktHos == 0]), c(phenoHos$Ethanol[mrpHos == 0 & mktHos == 1]))
segStats <- do.call('rbind', lapply(1:length(segGrowth), function(x) {
	segSub <- segGrowth[[x]]
	subMean <- mean(segSub)
	subSd <- sd(segSub)
	result <- as.data.frame(matrix(data=c(subMean, subMean-2*subSd, subMean+2*subSd ), ncol=3))
}))
colnames(segStats) <- c('mean', 'boundLow', 'boundHigh')
#plot
dev.new(width=6, height=4)
par(mar=c(1,1,1,1), oma=c(1,1,1,1), family='sans')
plot(x=-100, y=-100, xlim=c(0.5, 4.5), ylim=c(0, max(unlist(wts, muts))*1.1), xlab='', ylab='', xaxt='n' )
#Note, bars in illustrator were used to replaced confidence intervals set by arrows, because thickening arrows via lwd makes the edge of the bars thicken, and CIs appear wider than they are
arrows(c(1,2,3,4), rep(c(segStats$boundLow[1], segStats$boundLow[2]), 2), c(1,2,3,4), rep(c(segStats$boundHigh[1], segStats$boundHigh[2]), 2), col=alpha('black', alpha=0.2), length=0.05, angle=90, code=3, lwd=2)
arrows(c(1,2,3,4), rep(c(segStats$boundLow[3], segStats$boundLow[4]), 2), c(1,2,3,4), rep(c(segStats$boundHigh[3], segStats$boundHigh[4]), 2), col=alpha('violetred3', alpha=0.2), length=0.05, angle=90, code=3, lwd=2)
plotVals <- lapply(1:4, function(x) {
	xValsWt <- jitter(rep(x, length(wts[[x]])), amount=0.2) #amount=0.075
	points(xValsWt, wts[[x]] , cex=0.7, col=alpha('black', alpha=0.4), pch=19)
	xValsMut <- jitter(rep(x, length(muts[[x]])), amount=0.2) #amount=0.075
	points(xValsMut, muts[[x]] , cex=0.7, col=alpha('violetred3', alpha=0.4), pch=19)
})
labelVals <- c(expression(italic('MKT1'^'BY')), expression(italic('MKT1'^'3S')), expression(italic('MKT1'^'BY')), expression(italic('MKT1'^'3S')) )
axis(1:4, at=c(1:4), labels=labelVals, cex.axis=0.7)
legendLabel <- c(expression(italic('MRP20')), expression(italic('mrp20-105E')))
legend('bottomright', bty='n', legend=legendLabel, pch=c(19), col=alpha(c('black', 'violetred3'), alpha=0.5), cex=0.8 )
points(x=c(1,2), y=c(1.4,1.4), type='l')
points(x=c(3,4), y=c(1.4,1.4), type='l')
text(expression(italic('BY')), x=1.5, y=1.45)
text(expression(italic('3S')), x=3.5, y=1.45)



# Fig. 3C   Phenotypic distribution of mrp20-105E segregants from engineered BY x 3S crosses fixed for mrp20-105E and  causal MKT1 SNP at 467,219
phenoHos <- phenoAll[grepl('B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),] #cross B random is original hos3∆s in which mrp20-105E was discovered on Chr IV-BY allele
genoHos <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% phenoHos$Sample])
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
mktHos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample]) )
phenoHosBy <- phenoHos[mrpHos == 0 & mktHos == 0,]
phenoHos3s <-  phenoHos[mrpHos == 0 & mktHos == 1,]
phenoNew <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E crosses
phenoNewBy <- pheno[grepl('D', pheno$Sample),]
phenoNew3s <- pheno[grepl('E', pheno$Sample),]
# plot
dev.new(width=6, height=4)
par(mar=c(2,1,1,0), oma=c(1,2,1,1), family='sans')
plot(x=-10, y=-10, xlim=c(0.5, 4.5), ylim=c(0, 1.2), main='', xlab='', ylab='', xaxt='n')
yVals <- list(phenoHosBy$Ethanol, phenoHos3s$Ethanol, phenoNewBy$Ethanol, phenoNew3s$Ethanol)
temp <- as.data.frame(unlist(yVals))
temp$group <- c(rep(1, length(yVals[[1]])), rep(2, length(yVals[[2]])), rep(3, length(yVals[[3]])), rep(4, length(yVals[[4]])))
boxplot(temp[,1] ~ temp[,2], data=temp, boxwex=0.5, add=TRUE, names=c('', '', '', ''), outline=FALSE)
plotPoints <- lapply(1:4, function(x) {
	xValsPoints <- jitter(rep(x, length(yVals[[x]])), amount=0.2) #amount=0.075
	points(xValsPoints, yVals[[x]] , cex=0.7, col=alpha('black', alpha=0.125), pch=19)
})
xLabels <- c(expression(italic('mrp20-105E MKT1'^'BY')), expression(italic('mrp20-105E MKT1'^'3S')))
mtext('Genotype', side=1, line=2)
axis(1, at=c(1:4), labels=rep(xLabels, 2), cex.axis=0.7)
points(x=c(1,2), y=c(1.1,1.1), type='l')
text(x=1.5, y=1.15, 'original')
points(x=c(3,4), y=c(1.1,1.1), type='l')
text(x=3.5, y=1.15, 'new')
mtext('Growth Relative to BY', side=2, line=2)



# Fig. 3D   Stacked linkage mapping plot for new crosses that show 16 QTL across for FR scans

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
# 4 forward regression scans detected peaks above multiple testing corrections
#pvals
pvals1 <- read.table('res~locus_scan1Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
pvals2 <- read.table('res1of7loci~locus_scan2Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
pvals3 <- read.table('res2of4loci~locus_scan3Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
pvals4 <- read.table('res3of2loci~locus_scan4Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
pvals <- cbind(pvals1[,3], pvals2[,3], pvals3[,3], pvals4[,3])
#peaks
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
loci$gp <- sapply(1:nrow(loci), function(x) {
	geno$gp[geno$c == loci$c[x] & geno$p == loci$p[x]]
})
###Plot stacked scans
cols <- rep(c('gray', 'black'), 9)[-18]
dev.new(width=7.5, height=6)
layout(1:4)
par(mar=c(0,1,0,1), oma=c(3,3,3,0), family='sans')
plotScans <- lapply(1:4, function(y) {
	yMax <- -log10(min(pvals[,y], na.rm=TRUE))*1.2
	plot(-10,-10, xlim=c(0, max(geno$gp)), ylim=c(0,yMax), xlab='', ylab='', xaxt='n', las=1)
	lociSub <- list(loci[1:7,], loci[8:11,], loci[12:13,], loci[14:16,])[[y]]
	abline(v= lociSub$gp, lty=3, col='darkgray')
	plotChrs <- lapply(1:17, function(x) {
		xVals <- geno$gp[geno$c == x]
		if (x ==2) {
			yVals <- rep(1, length(geno$gp[geno$c==2]))
		} else if (x == 17) {
			yVals <- rep(1, length(geno$gp[geno$c==17]))
		}else {
			yRaw <- pvals[,y][geno$c == x]
			yRaw[is.na(yRaw)] <- 1
			yVals <- yRaw
		}		
		prevX <- max(geno$gp[geno$c == x-1])
		prevY <- 1
		nextX <- min(geno$gp[geno$c == x+1])
		nextY <- 1
		xChr <- c(prevX, xVals, nextX)
		yChr<- c(prevY, yVals, nextY)
		points(x=xChr, y=-log10(yChr), pch=19, col=cols[x], type='l')
	})
	text(max(geno$gp)*0.975,yMax*0.9, labels=paste('iteration', y, sep=' '))
})
chrMids <- sapply(1:17, function(x) {
	vals <- geno$gp[geno$c == x]
	mid <- median(c(min(vals), max(vals)))
})
chrLabels <- as.roman(1:17)
axis(1, at=chrMids[1:17], labels= chrLabels, cex.axis=0.7)
mtext('Genomic position', side=1, line=2)
mtext(expression('-log'['10']*'(pval)'), side=2, line=1, outer=TRUE, las=0)




# Fig. 3E   Chr XIV MKT1 and SAL1 effects 

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
sal1 <- as.factor(geno[geno$c == 14 & geno$p == 473648, colnames(geno) %in% pheno$Sample])
#make new 4 factor level that accounts for both closely linked SNPs in MKT1 and SAL1
both <- rep(NA, length(mkt1))
both[mkt1 == 0 & sal1 == 0] <- 0
both[mkt1 == 0 & sal1 == 1] <- 1
both[mkt1 == 1 & sal1 == 0] <- 2
both[mkt1 == 1 & sal1 == 1] <- 3
both <- as.factor(both)
dev.new(width=6, height=4)
par(mar=c(2,1,0,0), oma=c(1,2,1,1), family='sans')
plot(x=-10, y=-10, xlim=c(0.5, 4.5), ylim=c(0, max(pheno$Ethanol)*1.1), main='', xlab='', ylab='', xaxt='n')
points(jitter(as.numeric(both), amount=0.25), pheno$Ethanol, pch=19, col=alpha('black', alpha=0.2), cex=0.7)
xLabels <- c(expression(italic('MKT1'^'BY'*'SAL1'^'BY')), expression(italic('MKT1'^'BY'*'SAL1'^'3S')), expression(italic('MKT1'^'3S'*'SAL1'^'BY')), expression(italic('MKT1'^'3S'*'SAL1'^'3S')))
axis(1, at=1:4, labels=xLabels, cex.axis=0.8)
mtext('Growth Relative to BY', side=2, line=2)
mtext('mrp20-105E Genotype', side=1, line=2)




# Fig. 3F   Aneuploidy

an <- read.table('DataS6_ChromosomeIIDuplication.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
an <- an$Chr2 # 0 is wt and 1 indicates presence of aneuploidy/ duplication event

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
sal1 <- as.factor(geno[geno$c == 14 & geno$p == 473648, colnames(geno) %in% pheno$Sample])
both <- rep(NA, length(mkt1))
both[mkt1 == 0 & sal1 == 0] <- 0
both[mkt1 == 0 & sal1 == 1] <- 1
both[mkt1 == 1 & sal1 == 0] <- 2
both[mkt1 == 1 & sal1 == 1] <- 3
both <- as.factor(both)
dev.new(width=6,height=3.5)
layout(matrix(data=c(1:4), ncol=4))
par(mar=c(1,0,1,0), oma=c(3,4,2,2))
plotPops <- lapply(1:4, function(x) {
	xVals <- an[which(both == levels(both)[x])]+1
	yVals <- pheno$Ethanol[which(both == levels(both)[x])]
	plot(x=-10, y=-10, xlab='', ylab='', xaxt='n', xlim=c(0.5,2.5), ylim=c(0, max(pheno$Ethanol)*1.1), yaxt='n')
	points(x=jitter(xVals, amount=0.2), y=yVals, pch=19, col=alpha('black', alpha=0.2), cex=0.7)
	axis(1, at=c(1,2), labels=c('WT', 'Aneuploid'))
	if (x == 1) {
		axis(2, at=c(seq(0,1.2, by=0.2)), labels=c(seq(0,1.2, by=0.2)))
	}
	plotLabels <- c(expression(italic('MKT1'^'BY'*'SAL1'^'BY')), expression(italic('MKT1'^'BY'*'SAL1'^'3S')), expression(italic('MKT1'^'3S'*'SAL1'^'BY')), expression(italic('MKT1'^'3S'*'SAL1'^'3S')))
	mtext(plotLabels[x], side=3, line=-1.5, cex=0.6)
})
mtext('Chr 2', side=1, outer=TRUE, line=2, cex=0.8)
mtext('Growth relative to BY', side=2, outer=TRUE, line=2, cex=0.8)






### Figure 4

# Fig. 4A   Petite frequency in parent strains

#First, get the petite frequency data
mito <- read.table('DataS5_PetiteFrequency.txt', header=TRUE, stringsAsFactors=FALSE)
delimit <- 0.001 # this was the maximum size for petite in BY WT and 3S WT control strains
mitoParent <- mito[mito$Type == 'parent',]
#get a list of all the replicates for mito 
mitoParentF <- sapply(1:nrow(mitoParent), function(x) {
	vals <- as.numeric(strsplit(mitoParent$ColonySizes[x], ',')[[1]])
	freq <- c(length(vals[vals <= delimit]) / length(vals))*100
})
#Now, vor each sample type, take replicates and sample across 1000 bootstraps to get CIs, in this order BY MRP20, 3S MRP20 BY mrp20, 3S mrp20
mitoParentBoot <- do.call('rbind', lapply(1:4, function(x) {
	parentStr <- c('BY','3S', 'BY', '3S')[x]
	mrpStr <- c('MRP20', 'MRP20', 'mrp20', 'mrp20')[x]
	indeces <- which(mitoParent$Sample == parentStr & mitoParent$MRP20 == mrpStr)
	vals <- mitoParentF[indeces]
	bootVals <- sapply(1:1000, function(x) {
		newVals <- sample(vals, length(vals), replace=TRUE)
		mean(newVals)
	})
	bootOrd <- bootVals[order(bootVals)]
	result <- as.data.frame(matrix(data=c(mean(vals), bootOrd[50], bootOrd[950]), ncol=3))
}))
colnames(mitoParentBoot) <- c('meanVal', 'boundLeft', 'boundR')
#Also, need the growth data from cloning seg file
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
parents <- cloning[cloning$Type == 'parent',]  
#get growth vlaues for MRP20 background then for mrp20-105E in this order:
wtB <- parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'MRP20' & parents$Gene == 'WT']
wtS <- parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'MRP20' & parents$Gene == 'WT']
mutB <- parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$Gene == 'MRP20']
mutS <- parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$Gene == 'MRP20']
growthVals <- list(wtB, wtS, mutB, mutS)
bootValsGrowth <- do.call('rbind', lapply(1:4, function(x) {
	vals <- growthVals[[x]]
	bootVals <- sapply(1:1000, function(x) {
		newVals <- sample(vals, length(vals), replace=TRUE)
		mean(newVals)
	})
	bootOrd <- bootVals[order(bootVals)]
	result <- as.data.frame(matrix(data=c(mean(vals), bootOrd[50], bootOrd[950]), ncol=3))
}))
colnames(bootValsGrowth) <- c('meanVal', 'boundLeft', 'boundR')
cutoff <- 0.3127594 #inclusive of worst original mrp20-105E mkt1-BY hos3∆ individual
dev.new(width=4, height=4.1)
par(mar=c(1,2.5,0,0), oma=c(3,1.5,2,2)) 
# dev.new(width=9, height=3.1)
# layout(matrix(1:3, ncol=3))
par(mar=c(1,2.5,0,0), oma=c(3,1.5,2,2)) #
#plot parents
parentPchs <- c(23, 19, 17, 15)
plot(-1,-1, xlim=c(0,1.3), ylim=c(0,100), xlab='', ylab='')
abline(v=cutoff, lty=3, col='darkgray')
##Add CIs on x axis: growth vals
addCIsGrowth <- sapply(1:4, function(x) {
	arrows(bootValsGrowth$boundLeft[x], mitoParentBoot$meanVal[x], bootValsGrowth$boundR[x],  mitoParentBoot$meanVal[x], col=alpha('black', alpha= 0.3), length=0.05, angle=90, code=3, lwd=2)
})
##Add CIs on y axis: mito instability F vals
addCIsPetiteF <- sapply(1:4, function(x) {
	arrows(bootValsGrowth$meanVal[x], mitoParentBoot$boundLeft[x], bootValsGrowth$meanVal[x], mitoParentBoot$boundR[x], col=alpha('black', alpha= 0.3), length=0.05, angle=90, code=3, lwd=2)
})
points(x= bootValsGrowth$meanVal, y=mitoParentBoot$meanVal, pch=parentPchs, cex=1.2, bg='black')
mtext('Growth relative to BY', side=1, line=2)
mtext('Petite Frequency (%)', side=2, line=2)
mtext('Parents and cloning strains', side=3, line=0)



# Fig. 4B   Petite frequency in BYx3S MRP20 segregants

mito <- read.table('DataS5_PetiteFrequency.txt', header=TRUE, stringsAsFactors=FALSE)
delimit <- 0.001 # this was the maximum size for petite in BY WT and 3S WT control strains
mitoSegWt <- mito[mito$Type == 'segregant' & mito$MRP20 == 'MRP20',]
pheno <- phenoAll[phenoAll$Sample %in% mitoSegWt$Sample,]
mitoSegWtF <- sapply(1:nrow(mitoSegWt), function(x) {
	vals <- as.numeric(unlist(strsplit(mitoSegWt$ColonySizes[x], ',')))
	freq <- c(length(vals[vals <= delimit]) / length(vals))*100
})
#make sure growth and petite freq vals are in same order
mitoSegWtFOrd <- mitoSegWtF[match(pheno$Sample, mitoSegWt$Sample)] 
# mitoSegWt$Sample[match(pheno$Sample, mitoSegWt$Sample)] == pheno$Sample #should all be true if in proper order
cutoff <- 0.3127594 #inclusive of worst original mrp20-105E mkt1-BY hos3∆ individual
dev.new(width=4, height=4.1)
par(mar=c(1,2.5,0,0), oma=c(3,1.5,2,2)) #
plot(-1,-1, xlim=c(0,1.3), ylim=c(0,100), xlab='', ylab='')
abline(v=cutoff, lty=3, col='darkgray')
points(pheno$Ethanol, mitoSegWtFOrd, pch=19, col=alpha('black', alpha = 0.3), cex=0.7)
mtext('Growth relative to BY', side=1, line=2)
mtext('Petite Frequency (%)', side=2, line=2)
mtext('BYx3S MRP20', side=3, line=0)



# Fig. 4C   Petite frequency in BYx3S mrp20 segregants

mito <- read.table('DataS5_PetiteFrequency.txt', header=TRUE, stringsAsFactors=FALSE)
delimit <- 0.001 # this was the maximum size for petite in BY WT and 3S WT control strains
mitoSegMut <- mito[mito$Type == 'segregant' & mito$MRP20 == 'mrp20',]
pheno <- phenoAll[phenoAll$Sample %in% mitoSegMut$Sample,]
mitoSegMutF <- sapply(1:nrow(mitoSegMut), function(x) {
	vals <- as.numeric(unlist(strsplit(mitoSegMut$ColonySizes[x], ',')))
	freq <- c(length(vals[vals <= delimit]) / length(vals))*100
})
#make sure growth and petite freq vals are in same order
mitoSegMutFOrd <- mitoSegMutF[match(pheno$Sample, mitoSegMut$Sample)]
# mitoSegMut$Sample[match(pheno$Sample, mitoSegMut$Sample)] == pheno$Sample #should all be true if in proper order
cutoff <- 0.3127594 #inclusive of worst original mrp20-105E mkt1-BY hos3∆ individual
dev.new(width=4, height=4.1)
par(mar=c(1,2.5,0,0), oma=c(3,1.5,2,2)) #
plot(-1,-1, xlim=c(0,1.3), ylim=c(0,100), xlab='', ylab='')
abline(v=cutoff, lty=3, col='darkgray')
points(pheno$Ethanol, mitoSegMutFOrd, pch=19, col=alpha('black', alpha = 0.3), cex=0.7)
mtext('Growth relative to BY', side=1, line=2)
mtext('Petite Frequency (%)', side=2, line=2)
mtext('BYx3S mrp20', side=3, line=0)



# Fig. 5A   Predicting growth and comparing to observed growth

#we model growth mrp20-105E segregants from  BYx3S crosses engineered for mrp20-105E and MKT1 SNP
pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#CDetected loci
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
##which allele is worse? (6) 3S are worse and (10) BY are worse
loci$worse <- sapply(1:nrow(loci), function(x) {
	calls <- geno[geno$c == loci$c[x] & geno$p == loci$p[x], colnames(geno) %in% pheno$Sample]
	avgB <- mean(pheno$Ethanol[calls == 0], na.rm=TRUE)
	avgS <- mean(pheno$Ethanol[calls == 1], na.rm=TRUE)
	c(0,1)[which(c(avgB, avgS) == min(c(avgB, avgS)))] #worse haplotype by avg
})
#put all detected loci in model
lociCalls <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- geno[geno$c == l$c & geno$p == l$p, colnames(geno) %in% pheno$Sample]
	as.numeric(genoRow)
})))
colnames(lociCalls) <- paste('l', 1:16, sep='_')
lociCalls$mkt1 <- as.numeric(as.character(mkt1))
an <- read.table('DataS6_ChromosomeIIDuplication.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
an <- an$Chr2 # 0 is wt and 1 indicates presence of aneuploidy/ duplication event
lociCalls$an <- an
#Put all necessary information for model in dataframe
removeInds <- unique(unlist(lapply(1:ncol(lociCalls), function(x) {
	which(is.na(lociCalls[,x]))
})))
lociCalls2 <- lociCalls[-removeInds,]
lociCalls2 <- lapply(lociCalls2, factor)
lociCalls2$eth <- pheno$Ethanol[-removeInds]
lociCalls2 <- as.data.frame(lociCalls2)
model2 <- lm(eth ~ mkt1 + l_1 + l_2 + l_3 + l_4 + l_5 + l_6 + l_7 + l_8 + l_9 + l_10 + l_11 + l_12 + l_13 + l_14 + l_15 + l_16 + an, data=lociCalls2) # full model with aneuploidy
predictVals <- predict(model2)
#original mrp20-105E segregants
phenoHos <- phenoAll[grepl('B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
genoHos <- genoAll[, colnames(genoAll) %in% c('c', 'p') | colnames(genoAll) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
mkt1Hos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample]) )
sal1Hos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 473648, colnames(genoHos) %in% phenoHos$Sample]) )
lociCallsHos <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- genoHos[genoHos$c == l$c & genoHos$p == l$p, colnames(genoHos) %in% phenoHos$Sample,]
	#result <- abs(l$worse - as.numeric(genoRow))
	as.numeric(genoRow)
})))
colnames(lociCallsHos) <- paste('l', 1:16, sep='_')
lociCallsHos$an <- rep(0, nrow(lociCallsHos))
lociCallsHos$mkt1 <- as.numeric(as.character(mkt1Hos))
newDataHos <- lociCallsHos
newDataHos <- lapply(newDataHos, factor)
levels(newDataHos$an) <- c(0,1)
predictHos <- predict(model2, newdata=newDataHos)
colsHos <- rep(NA, length(mrpHos))
colsHos[lociCallsHos$mkt1 == 0 & lociCallsHos$l_6 == 0] <- colsBlOr[1]
colsHos[lociCallsHos$mkt1 == 1 & lociCallsHos$l_6 == 1] <- colsBlOr[2]
colsHos[is.na(colsHos)] <- 'black'
#plot
library(scales)
colsBlOr <- c("#377EB8", "#FF7F00")
dev.new(width=7.5, height = 4)
layout(matrix(data=c(1:2), nrow=1))
par(mar=c(0,0,0,1), oma=c(3,3,2,1))
#color by XIV parental genotype blue for BY ahplotype at both XIV SNPs vs orange for 3S at both; recombinants are black
mkt1 <- lociCalls2$mkt1
sal1 <- lociCalls2$l_6
both <- rep(NA, length(mkt1))
both[mkt1 == 0 & sal1 == 0] <- 0
both[mkt1 == 0 & sal1 == 1] <- 1
both[mkt1 == 1 & sal1 == 0] <- 2
both[mkt1 == 1 & sal1 == 1] <- 3
both <- as.factor(both)
colorsVectordeNovo <- c(colsBlOr[1], 'black', 'black', colsBlOr[2])[as.numeric(both)]
colorsVectordeNovo <- colorsVectordeNovo[-removeInds] #only plot inds with predictions and observations
predictVals[predictVals < 0] <- 0
plot(-1,-1, xlim=c(0, max(pheno$Ethanol)*1.025), ylim=c(0, max(pheno$Ethanol)*1.025), xlab='', ylab='')
abline(a=0, b=1, lty=3, col=alpha('black', alpha=0.2))
points(x=pheno$Ethanol[-removeInds][colorsVectordeNovo == 'black'], y=predictVals[colorsVectordeNovo == 'black'], pch=19, col=alpha(colorsVectordeNovo[colorsVectordeNovo == 'black'], alpha=0.1), cex=0.7)
points(x=pheno$Ethanol[-removeInds][colorsVectordeNovo != 'black'], y=predictVals[colorsVectordeNovo != 'black'], pch=19, col=alpha(colorsVectordeNovo[colorsVectordeNovo != 'black'], alpha=0.3), cex=0.7)
mtext(expression('New'), side=3, line=0)
legendLabels <- c(expression(italic('mrp20-105E MKT1'^'BY'*'SAL1'^'BY')), expression(italic('mrp20-105E MKT1'^'3S'*'SAL1'^'3S')), expression(italic('mrp20-105E MKT1'^'BY'*'SAL1'^'3S')), expression(italic('mrp20-105E MKT1'^'3S'*'SAL1'^'BY')), expression(italic('BY mrp20-105E')), expression(italic('3S mrp20-105E')))
legend('topleft', bty='n', pch=c(19, 19, 19, 19, 17, 15), col=c(alpha(colsBlOr, alpha=0.3), rep(alpha('black', alpha=0.7), 4)), legend=legendLabels, cex=0.7, pt.bg='black')
#add hos3∆ mrp20-105E segregants
predictHos[predictHos < 0] <- 0
plot(-1,-1, xlim=c(0, max(pheno$Ethanol)*1.025), ylim=c(0, max(pheno$Ethanol)*1.025), xlab='', ylab='', yaxt='n')
abline(a=0, b=1, lty=3, col=alpha('black', alpha=0.2))
points(x= phenoHos$Ethanol[mrpHos == 0], y= predictHos[mrpHos == 0], pch=19, col=alpha(colsHos[mrpHos == 0], alpha=0.3), cex=0.7)
mtext(expression('Original and mutant parents'), side=3, line=0)
mtext('Growth relative to BY', side=1, line=2, outer=TRUE)
mtext('Predicted Growth', side=2, line=2, outer=TRUE)
###Add predicted parents
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
parents <- cloning[cloning$Type == 'parent',]  
muts <- list(parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$MKT1 == 1])
parentalMutGrowth <- c(mean(muts[[1]]), mean(muts[[2]]))
lociP <- as.data.frame(matrix(data=c(rep(0,16), rep(1,16)), nrow=2, byrow=TRUE)) #first row is BY its just the genotype thats worse
lociP[3,] <- as.numeric(lociP[1,])
lociP[4,] <- as.numeric(lociP[2,])
colnames(lociP) <- colnames(lociCalls2)[1:16]
lociP$mkt1 <- c(0,1, 1,0)
lociP$an <- c(0,0,0,0)
#lociP$an <- c(1,0,1,1)
lociP$l_15_Int <- c(0,1,0,1)
lociP$strain <- c('BY', '3S', 'BY_Mkt13S', '3S_Mkt1BY')
levels(lociP$an) <- c(0,1)
lociP <- lociP[1:2,]
lociP <- lapply(lociP, factor)
levels(lociP$an) <- c(0,1)
predictP <- predict(model2, newdata=lociP)
points(parentalMutGrowth, predictP, col='black', pch=c(17,15), cex=1.0)



# Fig. 5B   Sum of detrimental alleles vs growth values for mrp20-105E segregants from  BYx3S crosses engineered for mrp20-105E and MKT1 SNP

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#CDetected loci
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
##which allele is worse? (6) 3S are worse and (10) BY are worse
loci$worse <- sapply(1:nrow(loci), function(x) {
	calls <- geno[geno$c == loci$c[x] & geno$p == loci$p[x], colnames(geno) %in% pheno$Sample]
	avgB <- mean(pheno$Ethanol[calls == 0], na.rm=TRUE)
	avgS <- mean(pheno$Ethanol[calls == 1], na.rm=TRUE)
	c(0,1)[which(c(avgB, avgS) == min(c(avgB, avgS)))] #worse haplotype by avg
})
#put all detected loci in model
lociCalls <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- geno[geno$c == l$c & geno$p == l$p, colnames(geno) %in% pheno$Sample]
	as.numeric(genoRow)
})))
colnames(lociCalls) <- paste('l', 1:16, sep='_')
lociCalls$mkt1 <- as.numeric(as.character(mkt1))
an <- read.table('DataS6_ChromosomeIIDuplication.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
an <- an$Chr2 # 0 is wt and 1 indicates presence of aneuploidy/ duplication event
lociCalls$an <- an
#in this case keep inds that are missing geno calls just for color assignment
mkt1 <- lociCalls$mkt1
sal1 <- lociCalls$l_6
both <- rep(NA, length(mkt1))
both[mkt1 == 0 & sal1 == 0] <- 0
both[mkt1 == 0 & sal1 == 1] <- 1
both[mkt1 == 1 & sal1 == 0] <- 2
both[mkt1 == 1 & sal1 == 1] <- 3
both <- as.factor(both)
colorsVectordeNovo <- c(colsBlOr[1], 'black', 'black', colsBlOr[2])[as.numeric(both)]
#Put all necessary information for model in dataframe
####now count the number of detrimental alleles carried by each segregant
calls <- lociCalls[,1:16]
calls <- lociCalls[,1:16][,-6]
sumBad <- sapply(1:nrow(calls), function(x) {
	callsInd <- as.numeric(calls[x,])
	#result <- sum(as.numeric(callsInd == loci$worse))	
	#result<- sum(as.numeric(callsInd == loci$worse)) + as.numeric(mkt1[x]== '0')
	sum(as.numeric(callsInd == loci$worse[-6]))
})
sumBad[sumBad < 3] <- 3
sumBad[sumBad > 12] <- 12
cutoff <- 0.3127594
#plot
dev.new(width=8, height=4)
par(mar=c(0,1,0,0), oma=c(3,2,2,1))
plot(-10,-10, xlim=c(3-0.25,max(sumBad, na.rm=TRUE)+0.25), ylim=c(0, max(pheno$Ethanol)*1.025), xlab='', ylab='', xaxt='n',cex.axis=0.8)
plotXVals <- lapply(0:max(sumBad, na.rm=TRUE), function(x) {
	index <- which(sumBad == x)
	yVals <- as.numeric(pheno$Ethanol[index]) 
	colsInds <- colorsVectordeNovo[index]
	xVals <- jitter(rep(x, length(yVals)), amount=0.2)
	points(x=xVals[colsInds == 'black'], y=yVals[colsInds == 'black'], col=alpha(colsInds[colsInds == 'black'], alpha=0.1), pch=19, cex=0.7 )
	points(x=xVals[colsInds != 'black'], y=yVals[colsInds != 'black'], col=alpha(colsInds[colsInds != 'black'], alpha=0.3), pch=19, cex=0.7 )
} )
abline(h=cutoff, lty=3, col='darkgray')
#add parental prediction values
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
parents <- cloning[cloning$Type == 'parent',]  
muts <- list(parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$MKT1 == 1])
parentalMutGrowth <- c(mean(muts[[1]]), mean(muts[[2]]))
lociP <- as.data.frame(matrix(data=c(rep(0,16), rep(1,16)), nrow=2, byrow=TRUE)) 
colnames(lociP) <- colnames(lociCalls2)[1:16]
lociP$mkt1 <- c(0,1)
lociP$an <- c(0,0)
lociP$strain <- c('BYmrp20', '3Smrp20')
callsP <- lociP[,1:16][,-6]
parentalSumBad <- c(sum(as.numeric(callsP[1,] == loci$worse[-6])), sum(as.numeric(callsP[2,] == loci$worse[-6])))
points(parentalSumBad, parentalMutGrowth, col='black', pch=c(17,15), cex=1.0)
axisLabels <- c('2-3', as.character(4:11), '12-13')
axis(1, at=c(range(sumBad, na.rm=TRUE)[1] : range(sumBad, na.rm=TRUE)[2]), labels = axisLabels) 
mtext('Sum of detrimental alleles', side=1, line=2)
mtext('Growth relative to BY', side=2, line=2)
legendLabels <- c(expression(italic('mrp20-105E MKT1'^'BY'*'SAL1'^'BY')), expression(italic('mrp20-105E MKT1'^'3S'*'SAL1'^'3S')), expression(italic('mrp20-105E MKT1'^'BY'*'SAL1'^'3S')), expression(italic('mrp20-105E MKT1'^'3S'*'SAL1'^'BY')), expression(italic('BY mrp20-105E')), expression(italic('3S mrp20-105E')))
legend('topright', bty='n', pch=c(19, 19, 19, 19, 17, 15), col=c(alpha(colsBlOr, alpha=0.3), rep(alpha('black', alpha=0.3), 2), rep(alpha('black', alpha=0.8), 2)), legend=legendLabels, cex=0.7, pt.bg='black')




















#####Supplementary Figures

# fig. S1    Initial discovery of the mrp201-05E mutation

##### Figure 1

# fig. S1A   Variance in WT vs hos3∆
# subset BYx3S F2 wildtype and BYx3S F2 hos3∆ random segregants 
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
growthWt <- pheno$Ethanol[grepl('A', pheno$Sample)] # A cross is WT segs
growthHos <- pheno$Ethanol[grepl('B', pheno$Sample)] #B cross is hos3∆ segs
combined <- c(growthWt, growthHos)
grouping <- as.factor(c(rep('Wt', length(growthWt)), rep('hos3∆', length(growthHos)) ))
dev.new(width=6, height=5)
par(mar=c(1,3,1,1), oma=c(1,1,1,1), family='sans')
plot(0,0, xlim=c(0.75,1.5), ylim=c(0,max(pheno$Ethanol)*1.1), xlab='', ylab='', xaxt='n')
points(x=jitter(rep(1, length(growthWt)), amount=0.075), y=growthWt, pch=19, cex=0.7, col=c(alpha('black', alpha=0.2)))
points(x=jitter(rep(1.25, length(growthHos)), amount=0.075), y= growthHos, pch=19, cex=0.7, col=c(alpha('black', alpha=0.2)))
axis(1, at=c(1,1.25), labels=c('HOS3', 'hos3'), tick=TRUE, cex.axis=0.7)
mtext('Growth relative to BY', side=2, line=2)


# fig. S1B   Linkage Mapping in BYx3S F2 wildtype and BYx3S F2 hos3∆ segregants

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, stringsAsFactors=FALSE)
result$gp <- genoAll$gp
genoMap <- result
cols <- rep(c('gray', 'black'), 9)[-18]
dev.new(width=7.5, height=5, unit='cm')
par(mar=c(1,1,1,1), oma=c(1,1,1,1), family='sans')
plot(-10,-10, xlim=c(0, max(genoMap$gp)), ylim=c(0, 58), xlab='', ylab='', xaxt='n', yaxt='n')
plotChrs <- sapply(1:17, function(x) {
	###add in dummy value at end of chromosome so that it connects to the first point of the next chromosome
	xChr <- genoMap$gp[genoMap$c == x]
	yChr <- genoMap$locusMutPvals[genoMap$c == x]
	if (x == 1) {
		prevX <- 0
	} else {
		prevX <-  max(genoMap$gp[genoMap$c == x-1])
	}
	if (x == 17) {
		nextX <- 0
	} else {
		nextX <- min(genoMap$gp[genoMap$c == x+1])
	}	
	prevY <- 1
	nextY <- 1
	xChr <- c(prevX, xChr, nextX)
	yChr<- c(prevY, yChr, nextY)
	points(x=xChr, y=-log10(yChr), pch=19, col=cols[x], type='l')
})
xLabVals <- c()
for (x in (1:17)) {
	xChr <- genoMap$gp[genoMap$c == x]
	xLabVals <- c(xLabVals, median(xChr))
}
axis(1, at=xLabVals, labels=as.roman(1:17), tick=TRUE, cex.axis=0.4)
axis(2, at=c(seq(0, 50, by=10)), labels=c(as.character(seq(0, 50, by=10))), tick=TRUE, cex.axis=0.5)
mtext('-log10 (pval) ', side=2, line=2, cex=0.7)





# fig. S1C   Chromosome IV effect in HOS3 vs hos3∆

#middle of peak
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, stringsAsFactors=FALSE)
result$gp <- genoAll$gp
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak4Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
mut <- as.factor(c('hos3', 'WT')[as.numeric(grepl('A', pheno$Sample)) +1]) #cross A is WTs, cross B is hos3∆s
#plot 
listPops <- list( pheno$Ethanol[peak4Calls == 0 & mut=='WT'], pheno$Ethanol[peak4Calls == 1 & mut=='WT'], pheno$Ethanol[peak4Calls == 0 & mut=='hos3'], pheno$Ethanol[peak4Calls == 1 & mut=='hos3'])
listX <- seq(1,1.75, by=0.25)
dev.new(width=7, height=5, unit='cm') #1/3 width of one col fig
par(mar=c(2,2,1,1), oma=c(1,1,1,1), family='sans')
plot(0,0, xlim=c(0.85,1.9), ylim=c(0,max(pheno$Ethanol)*1.1), xlab='', ylab='', xaxt='n')
cols <- c('black', 'black', 'black', 'black')
plotPops <- lapply(1:4, function(x) {
	yvals <- listPops[[x]]
	xVals <- jitter(rep(listX[x], length(yvals)), amount=0.075)
	points(x=xVals, y=yvals, pch=19, cex=0.7, col=alpha(cols[x], alpha=0.2))
})
mtext('Growth relative to BY', side=2, line=2)
axis(1, at=listX, labels=c('BY', '3S', 'BY', '3S'), cex.axis=0.7)
mtext('Genotype at Chromsome IV Locus', side=1, line=2)
points(x=c(listX[1], listX[2]), y=c(1.45, 1.45), type='l')
text(median(listX[1:2]), y=1.5, 'WT')
points(x=c(listX[3], listX[4]), y=c(1.45, 1.45), type='l')
text(median(listX[3:4]), y=1.5, 'hos3∆')





# fig. S1D   Plot recombination at Chromosome IV peak

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, stringsAsFactors=FALSE)
result$gp <- genoAll$gp
genoMap <- geno
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, stringsAsFactors=FALSE)
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak4 <- geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]
mut <- as.factor(c('hos3', 'WT')[as.numeric(grepl('A', pheno$Sample)) +1]) #cross A is WTs, cross B is hos3∆s
wtBy <- colnames((peak4[, peak4 == 0 & grepl('A', colnames(peak4))])) #A cross is WTs
wt3s <- colnames((peak4[, peak4 == 1 & grepl('A', colnames(peak4))]))
hosBy <- colnames((peak4[, peak4 == 0 & grepl('B', colnames(peak4))])) #B cros is hos3∆
hos3s <- colnames((peak4[, peak4 == 1 & grepl('B', colnames(peak4))]))
###go into and plot the physical breaks of the hosBy vs hos3s
# 1272234 #start of Dit2 and left bound; 1283799 end of PDR15 and right bound; 11.5 kb region to examine
#peak is 4 snps from 1277231 to 1277378 (pr and coding Mrp20)
bounds <- c(1274602, 1283799)
bound <- 6000
region <- genoMap[genoMap$c == 4 & genoMap$p >= bounds[1] & genoMap$p <= bounds[2],]
sub <- region[, colnames(region) %in% pheno$Sample] #55 rows
subHosBy <- sub[,colnames(sub) %in% hosBy] #90 inds
subHos3s <- sub[,colnames(sub) %in% hos3s] #130 inds
####
sub <- region[, colnames(region) %in% pheno$Sample] #55 rows
subHosBy <- sub[,colnames(sub) %in% hosBy] #90 inds
subHos3s <- sub[,colnames(sub) %in% hos3s] #130 inds
phenoBy <- pheno[pheno$Sample %in% hosBy,]
indsBy <- phenoBy$Sample[order(phenoBy$Ethanol)] #ordered worse to better growers
indsBy <- indsBy[c(1:80,90,81:89)] #reorder since last BY has an informative recombination breakpoint
pheno3s <- pheno[pheno$Sample %in% hos3s,]
inds3s <- pheno3s$Sample[order(pheno3s$Ethanol)] 
genoBy <- subHosBy[match(indsBy, colnames(subHosBy))]
geno3s <- subHos3s[match(inds3s, colnames(subHos3s))]
#concatenate groups to put worse growing class on bottom, by then 3s on top
genoCalls <- list(genoBy, geno3s) #ok, so now if you move list1 (worst class iv-by) col1 to col n worst to best grower
snpLength <- nrow(genoBy)
indLength <- c(ncol(genoCalls[[1]]) + ncol(genoCalls[[2]]) )
#physical size of plot
bp <- max(bounds)-min(bounds)
#use image to get a nice graphic of the matrix of calls
snps <- cbind(region[,1:3], genoBy, geno3s)
#fill in rows of redundant calls not on snps and convert snp table ito physical base by base call table
y <- matrix(nrow=0, ncol=ncol(snps))
for (x in c(bounds[1]:c(bounds[2]))) {
	print(paste(bounds[2]- x, 'rows left'))
	if (x %in% snps$p) { #snp row exists
		y <- rbind(y, snps[snps$p == x,])
	} else { #non snp row fill with NAs
		blank <- matrix(data=rep(NA, ncol(snps)), ncol=ncol(snps))
		blank[,1] <- 4
		blank[,2] <- x
		blank[,3] <- NA 
		colnames(blank) <- colnames(snps)
		y <- rbind(y, blank)
	}
}
y <- y[1: 9197,]
##Add in the de Novo mutation at position 1277959 coding () (A in BY and C in ref and 3S) 
#all inds impute values here wtih excpetion of delimiter F2_B_MRP20_hos3_random_25 which breaks right after snp from BY haplotype to 3S haplotype
y[y$p == 1277959, grepl('F2_B_MRP20_hos3_random_25', colnames(y))] <- 1
#Now need to fill in NAs; force the breaks to be exactly in the middle of 
fillInCalls <- function(inputVector) {
	calls <- inputVector
	test2 <- na.fill(na.approx(inputVector, na.rm=FALSE), "extend")
	test2[test2 < 0.5 & test2 != 0] <- 0
	test2[test2 >= 0.5 & test2 != 1] <- 1
	#fill in leading or trailing NAs
	test2 <- na.locf(test2)
	test2
}
allPos <- do.call('cbind', lapply(4:ncol(y), function(x) {
	newCalls <- fillInCalls(y[,x])
}))
colnames(allPos) <- colnames(y)[4:ncol(snps)]
final4 <- cbind(y[,1:2], allPos)
matrixPos <- t(t(final4[,colnames(final4) %in% pheno$Sample]))
#plot physical breaks for 4
newCol <- brewer.pal(5, 'Set1')[c(2,5)]
dev.new(width=4.6, height=3.2)
par(mar=c(2,2,0.5,0.5), oma=c(0,0,0,0))
image(matrixPos, xaxt='n', yaxt='n', col=newCol)
#add two by sepcific snps back into this
newSnps <- c(snps$p, 1277324, 1277959)[order(c(snps$p, 1277324, 1277959))]
xLabPos <- ((final4$p[final4$p %in% newSnps]) - min(final4$p)) / bp
axis(1, at=xLabPos, labels=FALSE, tick=TRUE)
###add lines to clearly distinguish the delimit bounds in middle of snps just like haplotypes
pos1 <- (xLabPos[which(newSnps == 1277231)] - xLabPos[which(newSnps == 1277231)-1 ]) / 2 + xLabPos[which(newSnps == 1277231)-1 ]
pos2 <- (xLabPos[which(newSnps == 1277959)+1] - xLabPos[which(newSnps == 1277959)]) / 2 + xLabPos[which(newSnps == 1277959) ]
abline(v= pos1, col='black', lty=3)
abline(v= pos2, col='black', lty=3)
mtext('Haplotype at peak position', side=2, line=1)
mtext('Chr IV 1,274,602 to 1,283,799', side=1, line=1)


# fig. S1E    Reciprocal Hemizygosity experiments to delimit Chr IV to MRP20
rhVals <- read.table('DataS3_ReciprocalHemizygosity.txt', header=TRUE, stringsAsFactors=FALSE)
#get confidence intervals via 1000 boostraps per measurement (each gene each LoF)
rhBoot <- as.data.frame(do.call('rbind', lapply(1:length(unique(rhVals$Gene)), function(x) {
	gene <- unique(rhVals$Gene)[x]
	valsBy <- rhVals$Ethanol[rhVals$Gene == gene & rhVals$LossOfFunction == 0] # BY allele
	vals3s <- rhVals$Ethanol[rhVals$Gene == gene & rhVals$LossOfFunction == 1] # 3S allele
	iterations <- as.data.frame(do.call('rbind', lapply(1:1000, function(y) {
		sampleBy <- sample(valsBy, length(valsBy), replace=TRUE)
		sample3s <- sample(vals3s, length(vals3s), replace=TRUE)
		c(mean(sampleBy), mean(sample3s))
	})))
	estBy <- c(mean(iterations[,1]), iterations[,1][order(iterations[,1])][50], iterations[,1][order(iterations[,1])][950])
	est3s <- c(mean(iterations[,2]), iterations[,2][order(iterations[,2])][50], iterations[,2][order(iterations[,2])][950])
	result <- as.data.frame(rbind(estBy, est3s))
	result$gene <- gene
	result$lof <- c('b', 's')
	result
})))
rhBoot <- rhBoot[order(rhBoot$gene),] #order from physical left to right of peak, happens to be alphabetical
colnames(rhBoot) <- c('mean', 'lowerBound', 'upperBound')
dev.new(width=5, height=5) 
par(mar=c(2,2,0,0), oma=c(1,1,1,1), family='sans')
plot(-10,-10, xlim=c(0, c(nrow(rhBoot))+1), ylim=c(0,1.1), xlab='', ylab='', xaxt='n')
arrows(c(1:nrow(rhBoot)), as.numeric(rhBoot[,2]), c(1:nrow(rhBoot)), as.numeric(rhBoot[,3]), length=0.05, angle=90, code=3)
axisLabels <- rep(c('3S', 'BY'), 3)
axis(1, at=c(1:6), labels=axisLabels)
points(x=c(1:nrow(rhBoot)), y=rhBoot[,1], pch=19) #plot estimates/ means from bootstrapping
points(x=c(1,2), y=c(1.025,1.025), type='l')
text(x=1.5, y=1.05, expression(italic('DIT1')))
points(x=c(3,4), y=c(1.025,1.025), type='l')
text(x=3.5, y=1.05, expression(italic('MRP20')))
points(x=c(5,6), y=c(1.025,1.025), type='l')
text(x=5.5, y=1.05, expression(italic('PDR15')))
mtext('Growth relative to BY', side=2, line=2)
mtext('Functional Allele', side=1, line=2)


# fig. S1F    Cloning of causal mrp20-105E nucleotide in segregants

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
segsMrp <- cloning[cloning$Type == 'segregant' & cloning$Gene == 'MRP20',]
allYs <- list( segsMrp$Ethanol[segsMrp$Edit == 'fromWTtoMut'], segsMrp$Ethanol[segsMrp$Edit == 'fromMuttoWT'])
dev.new(width=4, height = 4, unit='cm')
par(mar=c(1,1,1,1), oma=c(2,2,2,2))
plot(x=-1, y=-1, xlim=c(0.5,2.5), ylim=c(0, 1.5), xaxt='n', xlab='', ylab='')
plotVals <- sapply(c(1:length(allYs)), function(x) {
	yVals <- allYs[[x]]
	xVals <- jitter(rep(x, length(yVals)), amount=0.2)
	points(x=xVals, y=yVals, col=alpha('black', alpha=0.4), pch=19, cex=0.7)
})
mtext('Growth relative to BY', side=2, line=2)
mtext('Genotype', side=1, line=2)
axis(1, side=1, at=c(1:2), labels=c('MRP20 to mrp20-105E', 'mrp20-105E to MRP20'), cex.axis=0.7)


# fig. S1G    Tetrad dissection of original BY/3S HOS3/hos3∆ MRP20/mrp20-105E diploid to examine the mrp20-105E vs MRP20 in WT vs hos3

pheno <- phenoAll[grepl('B', phenoAll$Sample) & grepl('dissected', phenoAll$Sample),] #diploid B dissected not random spore prep
geno <- genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp') | colnames(genoAll) %in% pheno$Sample ]
mrp <- as.factor(as.numeric(geno[geno$c == 4 & geno$p == 1277282, colnames(geno) %in% pheno$Sample]))
hos <- as.factor(c('hos3', 'HOS3')[as.numeric(grepl('HOS3', pheno$Sample))+1]) #match val is 1, add 1 pulls second index 
allYs <- list( pheno$Ethanol[mrp == 0 & hos == 'hos3'], pheno$Ethanol[mrp == 1 & hos == 'hos3'], pheno$Ethanol[mrp == 0 & hos == 'HOS3'], pheno$Ethanol[mrp == 1 & hos == 'HOS3'])
dev.new(width=7, height=5)
par(mar=c(1,1,1,1), oma=c(2,2,2,2))
plot(x=-1, y=-1, xlim=c(0.5,4.5), ylim=c(0, 1.55), xaxt='n', xlab='', ylab='')
plotVals <- sapply(c(1:length(allYs)), function(x) {
	yVals <- allYs[[x]]
	xVals <- jitter(rep(x, length(yVals)), amount=0.2)
	points(x=xVals, y=yVals, col=alpha('black', alpha=0.4), pch=19, cex=0.7)
})
mtext('Growth relative to BY', side=2, line=2)
mtext('Genotype', side=1, line=2)
axis(1, side=1, at=c(1:4), labels=c('mrp20-105E', 'MRP20', 'mrp20-105E', 'MRP20'), cex.axis=0.7)
points(x=c(1,2), y=c(1.525, 1.525), type='l')
text(x=1.55, y=1.57, 'HOS3')
points(x=c(3,4), y=c(1.525, 1.525), type='l')
text(x=3.5, y=1.57, 'hos3∆')





# fig. S3   Linkage mapping in F3 intercross to resolve Chr XIV locus 

#A) genome wide significance plot in F3 population 

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_intercrossLocusPvalNew.txt', header=TRUE, stringsAsFactors=FALSE)
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak14Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
result$gp <- genoAll$gp
genoPlot <- result
genoPlot$locusPval[is.na(genoPlot$locusPval)] <- 1
cols <- rep(c('gray', 'black'), 9)[-18]
dev.new(width=7.5, height=5, unit='cm')
plot(-10,-10, xlim=c(min(genoPlot$gp), max(genoPlot$gp)), ylim=c(0, 45), xlab='', ylab='', xaxt='n', yaxt='n')
plotChrs <- sapply(1:17, function(x) {
	xChr <- genoPlot$gp[genoPlot$c == x]
	yChr <- genoPlot$locusPval[genoPlot$c == x]
	points(x=xChr, y=-log10(yChr), pch=19, col=cols[x], type='l')
})
xLabVals <- c()
for (x in (1:17)) {
	xChr <- genoPlot$gp[genoPlot$c == x]
	xLabVals <- c(xLabVals, median(xChr))
}
axis(1, at=xLabVals, labels=as.roman(1:17), tick=TRUE, cex.axis=0.4)
axis(2, at=c(seq(0, 45, by=5)), labels=c(as.character(seq(0, 45, by=5))), tick=TRUE, cex.axis=0.5)
mtext('Genomic position', side=1, line=2)
mtext('-log10(pval)', side=2, line=2)


# B) Plot effect of Chr XIV locus on growth 

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_intercrossLocusPvalNew.txt', header=TRUE, stringsAsFactors=FALSE)
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak14Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
allYs <- list( pheno$Ethanol[peak14Calls == 0], pheno$Ethanol[peak14Calls == 1])
dev.new(width=4, height = 4.25)
par(mar=c(1,1,1,1), oma=c(2,2,2,2))
plot(x=-1, y=-1, xlim=c(0.5,2.5), ylim=c(0, 2), xaxt='n', xlab='', ylab='')
plotVals <- sapply(c(1:length(allYs)), function(x) {
	yVals <- allYs[[x]]
	xVals <- jitter(rep(x, length(yVals)), amount= 0.25)
	points(x=xVals, y=yVals, col=alpha('black', alpha=0.2), pch=19, cex=0.7)
})
axisLabels  <- c(expression(italic('MKT1'^'BY')), expression(italic('MKT1'^'3S')))
axis(1, side=1, at=c(1:2), labels=axisLabels, cex.axis=0.8)
mtext('Growth relative to BY', side=2, line=2)
mtext('Genotype at Chr XIV 467,219', side=1, line=2)





# fig. S4   16QTL in supplement by rank

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#CDetected loci
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
##which allele is worse? (6) 3S are worse and (10) BY are worse
loci$worse <- sapply(1:nrow(loci), function(x) {
	calls <- geno[geno$c == loci$c[x] & geno$p == loci$p[x], colnames(geno) %in% pheno$Sample]
	avgB <- mean(pheno$Ethanol[calls == 0], na.rm=TRUE)
	avgS <- mean(pheno$Ethanol[calls == 1], na.rm=TRUE)
	c(0,1)[which(c(avgB, avgS) == min(c(avgB, avgS)))] #worse haplotype by avg
})
loci$diff <- sapply(1:nrow(loci), function(x) {
	calls <- geno[geno$c == loci$c[x] & geno$p == loci$p[x], colnames(geno) %in% pheno$Sample]
	avgB <- mean(pheno$Ethanol[calls == 0], na.rm=TRUE)
	avgS <- mean(pheno$Ethanol[calls == 1], na.rm=TRUE)
	abs(avgB-avgS)
})
loci$rank <- rank(loci$diff) ##Note, rank puts smallest first but we want to plot from laregest to smallest effect size
#get actual genotype calls at each locus
lociCalls <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- geno[geno$c == l$c & geno$p == l$p, colnames(geno) %in% pheno$Sample]
	as.numeric(genoRow)
})))
colnames(lociCalls) <- paste('l', 1:16, sep='_')
lociCalls$mkt1 <- as.numeric(as.character(mkt1))
# plot
dev.new(width=7, height=6.5)
layout(matrix(data=c(1:16), nrow=4, byrow=TRUE))
par(mar=c(0,0,0,0), oma=c(3,4,1,1))
plotLoci <- sapply(16:1, function(x) { 	
	locus <- loci[loci$rank == x,]
	whichLocus <- which(loci$rank == x)
	xVals <- lociCalls[, whichLocus]
	calls <- xVals
	xVals[xVals == 0] <- 0.25
	xVals[xVals == 1] <- 0.35
	xVals <- jitter(xVals, amount=0.025) 
	yVals <- pheno$Ethanol
	cols <- 'black'
	alphaCols <- 0.07
	if (x %in% c(16,12,8,4)) {
		plot(x=-1, y=-1, xlim=c(0.18, 0.42), ylim=c(0, 1.1), xlab='', ylab='', xaxt='n')
	} else {		
		plot(x=-1, y=-1, xlim=c(0.18, 0.42), ylim=c(0, 1.1), xlab='', ylab='', xaxt='n', yaxt='n')
	} 
	points(x=xVals, y=yVals, pch=19, col=alpha(cols, alpha=alphaCols), cex=0.7)	
	popCoords <- list (which(calls == 0 ), which(calls == 1))
	
	temp <- as.data.frame(matrix(data=c(yVals[popCoords[[1]]], yVals[popCoords[[2]]]), ncol=1))
	temp[,2] <- as.factor(c(rep(0.25, length(popCoords[[1]])), rep(0.35, length(popCoords[[2]]))))
	boxplot(temp[,1] ~ temp[,2], data=temp, boxwex=0.05, add=TRUE, names=c('', ''), outline=FALSE, at=c(0.25, 0.35), yaxt='n', notch=FALSE)
	
	if (x < 5) {
		axis(1, at=c(0.25, 0.35), labels=c('BY', '3S'))
	}
})
mtext('Growth relative to BY', side=2, line=2, outer=TRUE)
mtext('Genotype', side=1, line=2, outer=TRUE)





# fig. S5   Sum of detrimental alleles at 16 loci for BYx3S WT segregants

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#CDetected loci
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
##which allele is worse? (6) 3S are worse and (10) BY are worse
loci$worse <- sapply(1:nrow(loci), function(x) {
	calls <- geno[geno$c == loci$c[x] & geno$p == loci$p[x], colnames(geno) %in% pheno$Sample]
	avgB <- mean(pheno$Ethanol[calls == 0], na.rm=TRUE)
	avgS <- mean(pheno$Ethanol[calls == 1], na.rm=TRUE)
	c(0,1)[which(c(avgB, avgS) == min(c(avgB, avgS)))] #worse haplotype by avg
})
#Get wildtype segregants
pheno <- phenoAll[grepl('A', phenoAll$Sample),] #cross A is original wildtype segrengats
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(geno[geno$c == 14 & geno$p == 467219, colnames(geno) %in% pheno$Sample]))
#put all detected loci in model
lociCalls <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- geno[geno$c == l$c & geno$p == l$p, colnames(geno) %in% pheno$Sample]
	as.numeric(genoRow)
})))
colnames(lociCalls) <- paste('l', 1:16, sep='_')
lociCalls$mkt1 <- as.numeric(as.character(mkt1))
#in this case keep inds that are missing geno calls just for color assignment
mkt1 <- lociCalls$mkt1
sal1 <- lociCalls$l_6
both <- rep(NA, length(mkt1))
both[mkt1 == 0 & sal1 == 0] <- 0
both[mkt1 == 0 & sal1 == 1] <- 1
both[mkt1 == 1 & sal1 == 0] <- 2
both[mkt1 == 1 & sal1 == 1] <- 3
colorsVectorWt<- c(colsBlOr[1], 'black', 'black', colsBlOr[2])[both + 1]
#Put all necessary information for model in dataframe
####now count the number of detrimental alleles carried by each segregant
calls <- lociCalls[,1:16]
calls <- lociCalls[,1:16][,-6]
sumBad <- sapply(1:nrow(calls), function(x) {
	callsInd <- as.numeric(calls[x,])
	#result <- sum(as.numeric(callsInd == loci$worse))	
	#result<- sum(as.numeric(callsInd == loci$worse)) + as.numeric(mkt1[x]== '0')
	sum(as.numeric(callsInd == loci$worse[-6]))
})
sumBad[sumBad < 3] <- 3
sumBad[sumBad > 12] <- 12
cutoff <- 0.3127594
#plot
dev.new(width=8, height=4)
par(mar=c(0,1,0,0), oma=c(3,2,2,1))
plot(-10,-10, xlim=c(3-0.25,max(sumBad, na.rm=TRUE)+0.25), ylim=c(0, max(pheno$Ethanol)*1.025), xlab='', ylab='', xaxt='n',cex.axis=0.8)
plotXVals <- lapply(0:max(sumBad, na.rm=TRUE), function(x) {
	index <- which(sumBad == x)
	yVals <- as.numeric(pheno$Ethanol[index]) 
	colsInds <- colorsVectorWt[index]
	xVals <- jitter(rep(x, length(yVals)), amount=0.2)
	points(x=xVals[colsInds == 'black'], y=yVals[colsInds == 'black'], col=alpha(colsInds[colsInds == 'black'], alpha=0.1), pch=19, cex=0.7 )
	points(x=xVals[colsInds != 'black'], y=yVals[colsInds != 'black'], col=alpha(colsInds[colsInds != 'black'], alpha=0.3), pch=19, cex=0.7 )
} )
abline(h=cutoff, lty=3, col='darkgray')
axisLabels <- c('2-3', as.character(4:11), '12-13')
axis(1, at=c(range(sumBad, na.rm=TRUE)[1] : range(sumBad, na.rm=TRUE)[2]), labels = axisLabels) 
mtext('Sum of detrimental alleles', side=1, line=2)
mtext('Growth relative to BY', side=2, line=2)
legendLabels <- c(expression(italic('MRP20 MKT1'^'BY'*'SAL1'^'BY')), expression(italic('MRP20 MKT1'^'3S'*'SAL1'^'3S')))
legend('topright', bty='n', pch=c(19, 19), col=c(alpha(colsBlOr, alpha=0.3)), legend=legendLabels, cex=0.7)
