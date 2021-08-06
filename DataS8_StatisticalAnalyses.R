### DataS8_StatisticalAnalyses.R
## All code relevant to statistical tests employed in Schell, et al. 2021 are contained in this script
## Code is arranged in order of presented data
## Note, code compatible, tested under R version 3.5.1 (2018-07-02) -- 'Feather Spray" on x86_64-apple-darwin15.6.0 (64-bit) platform and R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" on x86_64-apple-darwin15.6.0 (64-bit)
# Necessary files for these analyses are 'DataS1_Genotype.txt', 'DataS2_Phenotype.txt', 'DataS3_ReciprocalHemizygosity.txt', 'DataS4_CloningExperiments.txt', 'DataS5_PetiteFrequency.txt', 'DataS6_ChromosomeIIDuplication.txt'
# Additionally, resultant linkage mapping files generated in DataS7.GeneticMappingAnalyses.R files are used for plotting, results, etc.

#### General overview/ layout of analyses in this file
# Statistical tests are presented as distinct sections of code. Results are presented in the same order as discussed in main text and supplement of Schell et al. 2021 
# Note that peak calling, and CI calculation for the final linkage mapping popultions with BY x 3S crosses that were engineered for mrp20-105E and the causal SNP in MKT1 at position 467,219 (Diploids D and E) are instead conducted in DataS7.GeneticMappingAnalyses.R for clarity, as subsequent linkage mapping scans require those calculations
##Include all relevant files in the working directory

library(car) # v.3.0.6
library(mixtools) # v.1.2.0

#####General functions to be used in subsequent code
#function getConfidenceInterval
#inputs: pvalues (string) and inGeno (dataframe) with  'c' and 'p' columns; both sdhold have same length and number of rows
#note: while loops -- be careful!  
#outputs: data frame with peak bounds and 2 lod drop bounds
getConfidenceInterval <- function(pvals, inGeno) {
	pvals[is.na(pvals)] <- 1
	peakAll <- which(pvals == min(pvals))
	peakRange <- range(inGeno$p[peakAll]) 
	peak <- as.integer(median(which(pvals == min(pvals))))
	lodVals <- -log10(pvals)
	lodPeak <- -log10(pvals[peak])[1]
	
	#Now get 2 lod drop; radiate outward from peak toward L and R until pval is found that is >= peak -2
	checkL <- TRUE
	checkR <- TRUE
	boundL <- peak
	boundR <- peak
	while(checkL == TRUE) {
		boundL <- boundL -1
		# print(boundL)
		if (lodVals[boundL] <= (lodPeak -2)) {
			boundL <- boundL +1
			checkL <- FALSE
		} 
	}
	while(checkR == TRUE) {
		boundR <- boundR +1
		# print(boundR)
		if (lodVals[boundR] <= (lodPeak -2)) {
			boundR <- boundR -1
			checkR <- FALSE
		} 
	}
	result <- as.data.frame(matrix(data=c(peakRange[1], peakRange[2], inGeno$p[boundL], inGeno$p[boundR]), nrow=1), stringsAsFactors=FALSE)
	colnames(result) <- c('peakLBound', 'peakRBound', 'boundL', 'boundR')
	result
}

# Read in general files used for all subsequent analyses
genoAll <- read.table('DataS1_Genotype.txt', header=TRUE, stringsAsFactors=FALSE)
phenoAll <- read.table('DataS2_Phenotype.txt', header=TRUE, stringsAsFactors=FALSE)




###### Below are distinct sections of analysis. ######





# # # ### Linkage Mapping in BYx3S F2 wildtype and BYx3S F2 hos3∆ segregants

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
#Levene's test to compare variance between WT and hos3∆ segregants
growthWt <- pheno$Ethanol[grepl('A', pheno$Sample)] # A cross is WT segs
growthHos <- pheno$Ethanol[grepl('B', pheno$Sample)] #B cross is hos3∆ segs
combined <- c(growthWt, growthHos)
grouping <- as.factor(c(rep('Wt', length(growthWt)), rep('hos3∆', length(growthHos)) ))
test <- leveneTest(combined, grouping) 
#test[[3]][1] # p = 1.420127e-22

#Linkage mapping 
#1) First, genome-wide scan to find loci with different effects in hos3∆s compared to wt segregants

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_locusMutPvalsNew.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
peak4 <- result[which(result[,3] == min(result[,3], na.rm=TRUE)),]
peak4Pval <- min(result[,3], na.rm=TRUE) # 2.537756e-54
#middle of peak
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak4Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
mut <- as.factor(c('hos3', 'WT')[as.numeric(grepl('A', pheno$Sample)) +1]) #cross A is WTs, cross B is hos3∆s
#examine how much variance is explained within hos3∆s by the peak
rsq <- summary(lm(pheno$Ethanol[mut == 'hos3'] ~ peak4Calls[mut == 'hos3']))$r.squared #0.74
p <- summary.aov(lm(pheno$Ethanol[mut == 'hos3'] ~ peak4Calls[mut == 'hos3']))[[1]][,5][1] # p = 4.317845e-66
#does this locus explain variance in wts?
pWt <- summary.aov(lm(pheno$Ethanol[mut == 'WT'] ~ peak4Calls[mut == 'WT']))[[1]][,5][1] # p = 0.91
###Get confidence interval as 2 lod drop around the peak
peak4Details <- getConfidenceInterval(result[,3], geno) #Note, peak breaks before de novo mutation since its not a SNP it isn't included in peak marker



##)Tetrad dissection of original BY/3S HOS3/hos3∆ MRP20/mrp20-105E diploid to see if hos3∆ is required for increased variance in originally observed hos3∆ mrp20-105E segregants

pheno <- phenoAll[grepl('B', phenoAll$Sample) & grepl('dissected', phenoAll$Sample),] #diploid B dissected not random spore prep
geno <- genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp') | colnames(genoAll) %in% pheno$Sample ]
mrp <- as.factor(as.numeric(geno[geno$c == 4 & geno$p == 1277282, colnames(geno) %in% pheno$Sample]))
hos <- as.factor(c('hos3', 'HOS3')[as.numeric(grepl('HOS3', pheno$Sample))+1]) #match val is 1, add 1 pulls second index 
#does effect depend on hos3?
pvalInt <- summary.aov(lm(pheno$Ethanol ~ mrp*hos))[[1]][,5][3] # intrxn term pval p = 0.76

# 1) Examining MRP20 vs mrp20-105E segregants

#phenotypic variation
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mrp <- as.factor(as.numeric(geno[geno$c == 4 & geno$p == 1277378, colnames(geno) %in% pheno$Sample]) )
growthWt <- pheno$Ethanol[grepl('A', pheno$Sample) | (grepl('B', pheno$Sample) & mrp == 1) ]
growthMut <- pheno$Ethanol[ (grepl('B', pheno$Sample) & mrp == 0) ]
combined <- c(growthWt, growthMut)
grouping <- as.factor(c(rep('Wt', length(growthWt)), rep('mrp', length(growthMut)) ))
#use Levene's test to see if variance differs between the two groupings
test <- leveneTest(combined, grouping)
pval <- test[[3]][1] #  5.891813e-22
fit <- normalmixEM(growthMut)
fit$loglik # 29.7
fit$mu #0.09786326 0.57461625
# 2) Linkage mapping in BYx3S F2 hos3∆ segregants detects Chr XIV effect

pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
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
#examine how much variance is explained within hos3∆s by the peak
rsq <- summary(lm(phenoHos$Ethanol[mrpHos == 0] ~ peak14Calls[mrpHos == 0]))$r.squared #0.79
p <- summary.aov(lm(phenoHos$Ethanol[mrpHos == 0] ~ peak14Calls[mrpHos == 0]))[[1]][,5][1] # p = 3.208216e-31
###Get confidence interval as 2 lod drop around the peak
peak14Details <- getConfidenceInterval(result[,3], geno)  



### Linkage Mapping in F3 Cross C to further resolve Chr XIV locus

pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
result <- read.table('geno_intercrossLocusPvalNew.txt', header=TRUE, stringsAsFactors=FALSE)
peak14 <- result[which(result[,3] == min(result[,3], na.rm=TRUE)),]
peak14Pval <- min(result[,3], na.rm=TRUE) # 2.503996e-43
peakMid <- result[as.integer(median(which(result[,3] == min(result[,3], na.rm=TRUE)))),]
peak14Calls <- as.factor(as.numeric(geno[geno$c == peakMid$c & geno$p == peakMid$p, colnames(geno) %in% pheno$Sample]))
###Get confidence interval as 2 lod drop around the peak
peak14Details <- getConfidenceInterval(result[,3], geno)



### Cloning causal nucleotides at MRP20 and MKT1 into BY and 3S parental strains

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
parents$MKT1 <- as.factor(parents$MKT1)
#What is the effect of mrp20-105E in each parent strain? Yes
pvalByMut <- summary.aov(lm(parents$Ethanol[parents$Sample =='BY' & parents$MKT1 == 0] ~ parents$MRP20[parents$Sample =='BY'& parents$MKT1 == 0])) [[1]][,5][1] # 4.313196e-25
pval3sMut <- summary.aov(lm(parents$Ethanol[parents$Sample =='3S' & parents$MKT1 == 1] ~ parents$MRP20[parents$Sample =='3S'& parents$MKT1 == 1])) [[1]][,5][1] # 4.0 x 10 -4
#Does Mkt1 influence the effect of mrp20-105E in each parent strain? Yes in 3S but not in BY
pvalByInt <- summary.aov(lm(parents$Ethanol[parents$Sample =='BY'] ~ parents$MRP20[parents$Sample =='BY']*parents$MKT1[parents$Sample =='BY']))[[1]][,5][3] # 0.99
pval3sInt <- summary.aov(lm(parents$Ethanol[parents$Sample =='3S'] ~ parents$MRP20[parents$Sample =='3S']*parents$MKT1[parents$Sample =='3S']))[[1]][,5][3] # 0.01



### Examine BYx3S mrp20-105E MKT1 engineered crosses

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
phenoBy <- pheno[grepl('D', pheno$Sample),]
pheno3s <- pheno[grepl('E', pheno$Sample),]
# Are means of these distributions different?
meanPval <- t.test(phenoBy$Ethanol, pheno3s$Ethanol)$p.value # 4.753501e-34

### Linkage mapping in  detected in new crosses and how much new loci explain repsonse to mrp20-105E

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
# how much phenotypic variance does the engineered mkt1 SNP account for in the new crosses?
rsqMkt1 <- summary(lm(pheno$Ethanol ~ mkt1))$r.squared # 0.1768972
#Cumulative and individual effects of the 16 loci detected
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
model <- lm(eth ~ mkt1 + l_1 + l_2 + l_3 + l_4 + l_5 + l_6 + l_7 + l_8 + l_9 + l_10 + l_11 + l_12 + l_13 + l_14 + l_15 + l_16, data=lociCalls2) # no aneuploidy in model
model2 <- lm(eth ~ mkt1 + l_1 + l_2 + l_3 + l_4 + l_5 + l_6 + l_7 + l_8 + l_9 + l_10 + l_11 + l_12 + l_13 + l_14 + l_15 + l_16 + an, data=lociCalls2) # full model with aneuploidy
#rsqModel <- summary(model2)$r.squared # 0.7269463/0.9296 = 78% (0.7819492); 0.9296 is the heritability of these new populations calculated by growth ~ replicates + error
#how much do loci, individually and collectively explain penothypic variation among mrp20-105E segregants?
modelInfo <- summary.aov(model)[[1]]
sumsqLoci <- modelInfo[,2][grepl('l_', rownames(modelInfo))]
all <- (sumsqLoci/sum(modelInfo[,2])) *100 #divide by total sum sq
#range(all[-(length(all))]) # 0.79 % to 14%  #excludes mkt1 and an; only examining 16 new loci
# most influential new locus
whichMax <- rownames(modelInfo)[grepl('l_', rownames(modelInfo))][which(all == max(all))] 
locusMax <- loci[as.numeric(strsplit(whichMax, '_')[[1]][2]),] # c 14 pos 473648; 'SAL1' locus
#Does the combined affect of mkt1 and sal1 deviate from additivity? Are they genetically interacting?
pvalInt <- summary.aov(lm(eth ~ mkt1*l_6, data=lociCalls2))[[1]][,5][3] # interaction term p = 0.77
# collective loci effects
rsqMkt <- modelInfo[,2][grepl('mkt1', rownames(modelInfo))] / sum(modelInfo[,2])*100
rsqSal <- modelInfo[,2][grepl('l_6', rownames(modelInfo))] / sum(modelInfo[,2])*100
restLoci <- sum(all[-which(all == max(all))]) # 35.78003 %
largeEffectXIV <- rsqMkt + rsqSal # 31.80887 %
# So, 15 other loci collectively have more phenotypic infleunce vs the two large effect loci on XIV at Mkt1 and Sal1
# Aneuploidy
an <- read.table('DataS6_ChromosomeIIDuplication.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
an <- an$Chr2 # 0 is wt and 1 indicates presence of aneuploidy/ duplication event
#difference in abundance between mkt1 by vs mkt1 3s
#sum(an[mkt1==0]) #21 w/ aneuploid vs 332 wt
#sum(an[mkt1==1]) #196 w/ aneuploidy vs 200 wt
counts <- as.data.frame(matrix(data=c(196, 200, 21, 332), ncol=2, byrow=TRUE))
colnames(counts) <- c('aneuploid', 'wt')
rownames(counts) <- c('MKT13S','MKT1BY')
fisher.test(counts)$p.value # 1.47555e-43
#aneuploidy significance 
#model accounting for mkt1 and an
rsqAn <- summary.aov(lm(pheno$Ethanol ~ as.factor(mkt1) +  as.factor(an)))[[1]][,5][2] # 4.759283e-19 for an term and model srsq is 0.26
modelPval <- pf(87.81, 3, 745, lower.tail=FALSE) # 1.171566e-48



### Test for interactions among the 16 loci detected in the BYx3S mrp20-105E MKT1 engineered crosses

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p', 'gp')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt1 <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
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

#Put all necessary information for model in dataframe
removeInds <- unique(unlist(lapply(1:ncol(lociCalls), function(x) {
	which(is.na(lociCalls[,x]))
	
})))
lociCalls2 <- lociCalls[-removeInds,]
lociCalls2 <- lapply(lociCalls2, factor)
lociCalls2$eth <- pheno$Ethanol[-removeInds]
lociCalls2 <- as.data.frame(lociCalls2)
pairs <- as.data.frame(t(combn(c(1:16), 2))) #120 combos
trios <- as.data.frame(t(combn(c(1:16), 3))) #560
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt1)))
testPairs <- sapply(1:nrow(pairs), function(x) {
	pair <- pairs[x,]
	l1 <- lociCalls2[,pair[1,1]]
	l2 <- lociCalls2[,pair[1,2]]
	pvalInt <-  summary.aov(aov(res[-removeInds] ~ l1*l2))[[1]][,5][3]
})
# use bonferroni correction as threshold
sigPairs <- which(testPairs <= (0.05/length(testPairs)) ) #none
testTrios <- sapply(1:nrow(pairs), function(x) {
	trio <- trios[x,]
	l1 <- lociCalls2[,trio[1,1]]
	l2 <- lociCalls2[,trio[1,2]]
	l3 <- lociCalls2[,trio[1,3]]
	pvalInt <-  summary.aov(aov(res[-removeInds] ~ l1*l2*l3))[[1]][,5][7]
})
#min is only 0.007067837 so won't pass bonferroni correction
pairsMkt1 <- sapply(1:16, function(x) {
	l1 <- lociCalls2[,x]
	l2 <- mkt1[-removeInds]
	pvalInt <-  summary.aov(aov(res[-removeInds] ~ l1*l2))[[1]][,5][3]
}) #min is 0.007593622
triosMkt1 <- sapply(1:16, function(x) {
	pair <- pairs[x,]
	l1 <- lociCalls2[,pair[1,1]]
	l2 <- lociCalls2[,pair[1,2]]
	l3 <- mkt1[-removeInds]
	pvalInt <-  summary.aov(aov(res[-removeInds] ~ l1*l2*l3))[[1]][,5][7]
}) #min is 0.03512427





### Examine petite frequency in parents, BYx3S MRP20 and BYx3S mrp20-105E segregants

mito <- read.table('DataS5_PetiteFrequency.txt', header=TRUE, stringsAsFactors=FALSE)
delimit <- 0.001 # this was the maximum size for petite in BY WT and 3S WT control strains
mitoParent <- mito[mito$Type == 'parent',]
#get a list of all the replicates for mito 
mitoParentF <- sapply(1:nrow(mitoParent), function(x) {
	vals <- as.numeric(strsplit(mitoParent$ColonySizes[x], ',')[[1]])
	freq <- c(length(vals[vals <= delimit]) / length(vals))*100
})
#make a list of the values in this order: BY MRP20, BY mrp20, 3S MRP20, 3S mrp20
mitoParentValsList <- list(mitoParentF[which(mitoParent$Sample == 'BY' & mitoParent$MRP20 =='MRP20')], mitoParentF[which(mitoParent$Sample == 'BY' & mitoParent$MRP20 =='mrp20')], mitoParentF[which(mitoParent$Sample == '3S' & mitoParent$MRP20 =='MRP20')], mitoParentF[which(mitoParent$Sample == '3S' & mitoParent$MRP20 =='mrp20')]
)
pvalBy <- t.test(mitoParentValsList[[1]], mitoParentValsList[[2]])$p.value # p = 0.01335452
pval3s <- t.test(mitoParentValsList[[3]], mitoParentValsList[[4]])$p.value # p = 0.3910022
#segregants values are clearly non normally distributed, so use nonparametric t. test
mito <- read.table('DataS5_PetiteFrequency.txt', header=TRUE, stringsAsFactors=FALSE)
delimit <- 0.001 # this was the maximum size for petite in BY WT and 3S WT control strains
mitoSegWt <- mito[mito$Type == 'segregant' & mito$MRP20 == 'MRP20',]
pheno <- phenoAll[phenoAll$Sample %in% mitoSegWt$Sample,]
mitoSegWtF <- sapply(1:nrow(mitoSegWt), function(x) {
	vals <- as.numeric(unlist(strsplit(mitoSegWt$ColonySizes[x], ',')))
	freq <- c(length(vals[vals <= delimit]) / length(vals))*100
})
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
segVals <- as.data.frame(matrix(c(mitoSegWtF, mitoSegMutF)))
segVals[,2] <- as.factor(c(rep('MRP20', length(mitoSegWtF)), rep('mrp20', length(mitoSegMutF))))
# Are mrp20-105E segregants differnt than MRP20 segregants? Yes
pvalWtvsMut <- wilcox.test(segVals[,1] ~ segVals[,2], exact=FALSE)$p.value # p = 0.02313469





### Model growth in mrp20-105E individuals with detected loci

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
# how well does th emodel do?
modelrsq <- summary(model2)$r.squared # 0.7269463 / /0.9296 = 78% (0.7819492); 0.9296 is the heritability of these new populations calculated by growth ~ replicates + error
modelPval <- pf(106, 18, 717, lower.tail=FALSE) # 5.183933e-188
# correlations between observed and predicted
predictVals <- predict(model2)
test <- cor.test(c(pheno$Ethanol[-removeInds]), c(as.numeric(predictVals)), method='pearson')
r <- test$estimate # 0.8526114
p <- test$p.value # 4.3755e-209
###get predictions for original mrp20-105E segregants found in hos3∆ segregants
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci <- rbind(scan1Peaks, scan2Peaks, scan3Peaks, scan4Peaks)
phenoHos <- phenoAll[grepl('B', phenoAll$Sample),]
genoHos <- genoAll[, colnames(genoAll) %in% c('c', 'p') | colnames(genoAll) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
mkt1Hos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample]) )
lociCallsHos <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- genoHos[genoHos$c == l$c & genoHos$p == l$p, colnames(genoHos) %in% phenoHos$Sample,]
	as.numeric(genoRow)
})))
colnames(lociCallsHos) <- paste('l', 1:16, sep='_')
lociCallsHos$an <- rep(0, nrow(lociCallsHos))
lociCallsHos$mkt1 <-as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample,])
newDataHos <- lociCallsHos
newDataHos <- lapply(newDataHos, factor)
levels(newDataHos$an) <- c(0,1)
predictHos <- predict(model2, newdata=newDataHos)
#We only want to examine the mrp201-05E individuals, remember in this cross the BY haplotype on IV has the spontaneous mutation
test <-cor.test(predictHos[mrpHos ==0], growthHos[mrpHos == 0])
test$estimate # 0.8980902
test$p.value # 1.605759e-39
##predictions for parent mrp20-105E strains 
cloning <- read.table('DataS4_CloningExperiments.txt', header=TRUE, stringsAsFactors=FALSE)
parents <- cloning[cloning$Type == 'parent',]  
muts <- list(parents$Ethanol[parents$Sample == 'BY' & parents$MRP20 == 'mrp20' & parents$MKT1 == 0],
parents$Ethanol[parents$Sample == '3S' & parents$MRP20 == 'mrp20' & parents$MKT1 == 1])
lociP <- as.data.frame(matrix(data=c(rep(0,16), rep(1,16)), nrow=2, byrow=TRUE)) 
colnames(lociP) <- colnames(lociCalls2)[1:16]
lociP$mkt1 <- c(0,1)
lociP$an <- c(0,0)
lociP$strain <- c('BYmrp20', '3Smrp20')
lociP <- lapply(lociP, factor)
levels(lociP$an) <- c(0,1)
predictP <- predict(model2, newdata=lociP)



### comparing model to wt segregants; use original random spore pops from diploids A and B

# read in mutants to generate model
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
# how well does th emodel do?
modelrsq <- summary(model2)$r.squared # 0.7269463 / /0.9296 = 78% (0.7819492); 0.9296 is the heritability of these new populations calculated by growth ~ replicates + error
modelPval <- pf(106, 18, 717, lower.tail=FALSE) # 5.183933e-188
# correlations between observed and predicted
predictVals <- predict(model2)
#WT segregants (diploid A)
phenoWt <- phenoAll[grepl('A', phenoAll$Sample),]
genoWt <- genoAll[, colnames(genoAll) %in% c('c', 'p') | colnames(genoAll) %in% phenoWt$Sample]
growthWt <- phenoWt$Ethanol
mkt1Wt <- as.factor(as.numeric(genoWt[genoWt$c == 14 & genoWt$p == 467219, colnames(genoWt) %in% phenoWt$Sample]) )
lociCallsWt <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- genoWt[genoWt$c == l$c & genoWt$p == l$p, colnames(genoWt) %in% phenoWt$Sample,]
	as.numeric(genoRow)
})))
colnames(lociCallsWt) <- paste('l', 1:16, sep='_')
lociCallsWt$an <- rep(0, nrow(lociCallsWt))
lociCallsWt$mkt1 <- mkt1Wt
newDataWt <- lociCallsWt
newDataWt <- lapply(newDataWt, factor)
levels(newDataWt$an) <- c(0,1)
predictWt <- predict(model2, newdata=newDataWt)
#How well do we predict in WT?
est <- cor.test(predictWt, growthWt)$estimate # 0.6924901
pval <- cor.test(predictWt, growthWt)$p.value # 9.560303e-25


#hos3∆ mrp20-105E segregants (diploid B)
phenoHos <- phenoAll[grepl('B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
genoHos <- genoAll[, colnames(genoAll) %in% c('c', 'p') | colnames(genoAll) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
mkt1Hos <- as.factor(as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample]) )
lociCallsHos <- as.data.frame(do.call('cbind', lapply(1:nrow(loci), function(x) {
	l <- loci[x,]
	genoRow <- genoHos[genoHos$c == l$c & genoHos$p == l$p, colnames(genoHos) %in% phenoHos$Sample,]
	as.numeric(genoRow)
})))
colnames(lociCallsHos) <- paste('l', 1:16, sep='_')
lociCallsHos$an <- rep(0, nrow(lociCallsHos))
lociCallsHos$mkt1 <-as.numeric(genoHos[genoHos$c == 14 & genoHos$p == 467219, colnames(genoHos) %in% phenoHos$Sample,])
newDataHos <- lociCallsHos
newDataHos <- lapply(newDataHos, factor)
levels(newDataHos$an) <- c(0,1)
predictHos <- predict(model2, newdata=newDataHos)
#run a model to see if departure between observed growth and predicted growth intearcts with MRP20 status
#Now hos3∆s are comprised of mrp20-105E and MRP20 individuals since mrp20-105E spontaneously arrose on IV BY allele in that releavnt diploid
temp <- as.data.frame(matrix(data=c(growthWt, growthHos), ncol=1))
colnames(temp) <- 'observed'
temp$predicted <- c(predictWt, predictHos)
temp$MRP20 <- as.factor(c(rep('MRP20', length(growthWt)), c('mrp20', 'MRP20')[as.numeric(mrpHos)]))
#Does departure from predicted and observed depend on MRP20 status?
pval <- summary.aov(lm(observed ~ predicted*MRP20, data=temp))[[1]][,5][3] #interaction is significant  2.834955e-23 
