### DataS7_GeneticMappingAnalyses.R
## All code relevant to genetic mapping and associated multiple testing correction Schell, et al. 2021 are contained in this script
## Code is arranged in order of presented data
## Note, code compatible, tested under R version 3.5.1 (2018-07-02) -- 'Feather Spray" on x86_64-apple-darwin15.6.0 (64-bit) platform and R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night" on x86_64-apple-darwin15.6.0 (64-bit)
# Necessary files for these analyses are the 'DataS1_Genotype.txt' and 'DataS2_Phenotype.txt' files
# Note, resulting files generated in this file are used for satistical analyses and figure plotting in DataS8.StatisticalAnalyses.R and DataS9.FigurePlotting.R files

#### General overview/ layout of analyses in this file
##Include all relevant files in the working directory
## 1) Linkage Mapping in BYx3S F2 wildtype and BYx3S F2 hos3∆ segregants
## 2) Linkage Mapping in BYx3S F2 hos3∆ segregants
## 3) Linkage Mapping in F3 mrp20-105E segregants
## 4) Linkage Mapping in BYx3S F2 mrp20-105E MKt1BY and MKT1 3S engineered crosses, permutations and peak calling
#Note, that genome wide scans generally take ~5-10 minutes
#However, permutation scans take much longer as they consist of 1000 genome wide scans. We broke these into subsets and ran them on usc's high performance computing cluseter, but they are included here for those interested. BE aware if you attempt to run them serially, this will take days to run.
# Note that peak calling, and CI calculation for the final linkage mapping popultions with BY x 3S crosses that were engineered for mrp20-105E and the causal SNP in MKT1 at position 467,219 (Diploids D and E) are conducted in DataS7.GeneticMappingAnalyses.R for clarity, as subsequent linkage mapping scans require those calculations. All other peak analysis is conducted n DataS8_StatisticalAnalyses.R


#read in general files, used for all code blocks below
genoAll <- read.table('DataS1_Genotype.txt', header=TRUE, stringsAsFactors=FALSE)
phenoAll <- read.table('DataS2_Phenotype.txt', header=TRUE, stringsAsFactors=FALSE)

###General functions

#callPeaks function for linakge mapping in BYx3S mrp20-105E fixed crosses D and E
#for calling peaks gneome wide, calls multiple peaks per chromsome, remove if peaks are too nearby after peak calling
# calculates 2 lod drop 99% CI around each peak by radiated outward until hitting first test below 2 lod from peak test
#note, excludes chr 2 due to aneuploidy of this chromosome in these crosses
callPeaks <- function(inGeno, inPvals, thresh, lodVal) {
	chrPeaks <- lapply(c(1, 3:16), function(x) {
		print(paste('scanning chr', x))
		peaksChr <- data.frame()
		checkChr <- TRUE
		while (checkChr == TRUE) {
			mostSig <- min(inPvals$genoRow[inGeno$c == x], na.rm=TRUE)
			if (-log10(mostSig) < thresh) {
				checkChr <- FALSE
			} else {
				pos <- which(inPvals$genoRow[inGeno$c == x] == mostSig)
				delimitVals <- range(inGeno$p[inGeno$c == x][pos])
				center <- pos[which(pos == as.integer(median(pos)))]
				diff <- -log10(mostSig) - c(-log10(inPvals$genoRow[inGeno$c == x]))
				diff[is.na(diff)] <- 0
				findL <- TRUE
				boundL <- min(pos)
				while (findL == TRUE) {
					if (boundL == 1) { 
						findL <- FALSE
					} else {
						boundL <- boundL - 1
					}
					if (diff[boundL] >= lodVal) {
						boundL <- boundL + 1
						findL <- FALSE
					}
				} 
				findR <- TRUE
				boundR <- max(pos)
				while (findR == TRUE) {
					if (boundR == length(diff)) {
						findR <- FALSE
					} else {
						boundR <- boundR + 1
					}
					if (diff[boundR] >= lodVal) {
						boundR <- boundR - 1
						findR <- FALSE
					}
				}
				# newPeak <- data.frame(matrix(data=c(x, inGeno$p[inGeno$c == x][center], inGeno$p[inGeno$c == x][boundL], inGeno$p[inGeno$c == x][boundR], mostSig, -log10(mostSig)), nrow=1), stringsAsFactors=FALSE)
				# colnames(newPeak) <- c('c', 'p', 'boundL', 'boundR', 'pval', 'lod')
				newPeak <- data.frame(matrix(data=c(x, inGeno$p[inGeno$c == x][center], delimitVals[1], delimitVals[2],	 inGeno$p[inGeno$c == x][boundL], inGeno$p[inGeno$c == x][boundR], mostSig, -log10(mostSig)), nrow=1), stringsAsFactors=FALSE)
				colnames(newPeak) <- c('c', 'p', 'peakL', 'peakR','boundL', 'boundR', 'pval', 'lod')
				
				if (nrow(peaksChr) == 0) {
					peaksChr <- newPeak
				} else {
					peaksChr <- rbind(newPeak, peaksChr)
				}
				checkChr <- FALSE
			}						
		}
		peaksChr
		})
	result <-do.call('rbind', chrPeaks)			
	}
threshScans <- c(4.305027, 4.325596, 4.274288, 4.332096) #alpha 0.05 or 50th smallest pval of 1000 permuted genome wide scans 




# 1)  Linkage Mapping in BYx3S F2 wildtype and BYx3S F2 hos3∆ segregants

#Subset out wt and hos3∆ segregants
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
#First, genome-wide scan to find loci with different effects in hos3∆s compared to wt segregants
growthEth <- pheno$Ethanol
mut <- as.factor(c('hos3', 'WT')[as.numeric(grepl('A', pheno$Sample)) +1]) #cross A is WTs, cross B is hos3∆s
genoMap <- geno

locusVals <- sapply(1:nrow(genoMap), function(x) {
	print(paste(x, 'of', nrow(genoMap)))
	row <- genoMap[x,][,colnames(genoMap) %in% pheno$Sample]
	genoVals <- as.factor(as.numeric(genoMap[x,][,colnames(genoMap) %in% pheno$Sample]))
	if (length(unique(genoVals)) > 1) {
		locusPval <- summary.aov(lm(growthEth ~ mut* genoVals))[[1]][,5][3]
	} else {
		locusPval <- NA
	}
	locusPval
})
out <- cbind(genoMap[,1:2], locusVals)
colnames(out)[ncol(out)] <- 'locusMutPvals'
#write.table(out, file='geno_locusMutPvalsNew.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')


#2) Linkage Mapping in BYx3S F2 hos3∆ segregants: genome wide scan within hos3∆s for loci that have affects that depend upon MRP20 state (wt or mrp20-105E)
#Subset out wt and hos3∆ segregants
pheno <- phenoAll[grepl('F2_A', phenoAll$Sample) | grepl('F2_B', phenoAll$Sample) & !grepl('dissected', phenoAll$Sample),]
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
#Now, only use hos3∆s, since we know they have mrp20-105E on IV-BY allele but WT do not
phenoHos <- pheno[grepl('B', pheno$Sample),]
genoHos <- geno[, colnames(geno) %in% c('c', 'p') | colnames(geno) %in% phenoHos$Sample]
growthHos <- phenoHos$Ethanol
mrpHos <- as.factor(as.numeric(genoHos[genoHos$c == 4 & genoHos$p == 1277378, colnames(genoHos) %in% phenoHos$Sample]) )
hosMrp20Int <- sapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.factor(as.numeric(genoHos[x, colnames(genoHos) %in% phenoHos$Sample]))
	if (length(unique(genoRow)) > 1) {
		result <- summary.aov(lm(growthHos ~ mrpHos*genoRow))[[1]][,5][3]
	} else {
		result <- NA
	}
	result
})
out <- cbind(genoMap[,1:2], hosMrp20Int)
colnames(out)[ncol(out)] <- 'mrp20LocusInt'
#write.table(out, file='mrp20*locus_withinhosMutOnly_natcommNew.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')



#3) Linkage mapping in F3 mrp20-105E segregants


pheno <- phenoAll[grepl('F3_C', phenoAll$Sample),] #Cross C is only cross to produce F3s
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
genoMap <- geno
growthE <- pheno$Ethanol
locusPvals <- sapply(1:nrow(genoMap), function(x) {
	print(paste(x, 'of', nrow(genoMap), sep=' '))
	row <- genoMap[x, colnames(genoMap) %in% pheno$Sample]
	locus <- as.factor(as.numeric(row))
	if (length(unique(locus)) > 1) {
		pval <-  summary.aov(lm(growthE ~ locus))[[1]][,5][1]
	} else {
		NA
	}
})
genoResult <- cbind(genoMap[,1:2], locusPvals)
colnames(genoResult)[ncol(genoResult)] <- 'locusPval'
# write.table(genoResult, file='geno_intercrossLocusPvalNew.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)



#4) Linkage mapping in BYx3S crosses fixed for mrp20-105-E and causal MKT1 variant at 467,219

#Scan1
pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#first thing, correct for effect of causal MKT1 variant
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt)))
#Scan 1 
#Now, run genome wide scan looking for loci affecting growth
scan1 <- do.call('rbind', lapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))
	if (length(unique(genoRow)) == 1) {
		pvals <- c(NA)
	} else {
		pvals <- summary.aov(lm(res ~ as.factor(genoRow)))[[1]][,5][1]
	} 
	if (length(pvals) != 1) {
		pvals <- c(NA)
	} else {
		pvals
	}
}))
scan1 <- as.data.frame(scan1)
out <- cbind(geno[,1:2], scan1[,1])
write.table(out, file='res~locus_scan1Pvals_New.txt', col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
#Scan 1 Permutations
# ####scaramble phenos 1000 for 1000 genome wide scans
# permPhenos <- do.call('cbind', lapply(1:1000, function(x) {
	# sample(res, length(res), replace=TRUE)
# }))
# permP <- as.data.frame(permPhenos)
# # write.table(permP, file='permutedPhenos_1000.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
# permRes <- read.table('permutedResPhenos_1000.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)
# scan1PermutedPvals <- do.call('cbind', lapply(1:1000, function(y) {
	# print(paste('permutation genome scan', y))
	# permVals <- permRes[,y]
	# scan <- do.call('rbind', lapply(1:nrow(geno[geno$c != 17,]), function(x) {
		# genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))
		# if (length(unique(genoRow)) == 1) {
			# c(NA)
		# } else {
			# pvals <- summary.aov(lm(permVals ~ genoRow))[[1]][,5][1]
		# } 
		# if (length(pvals) != 1) {
			# c(NA)
		# } else {
			# pvals
		# }
		
	# }))
# }))
# scan1PermutedPvals <- as.data.frame(scan1PermutedPvals)
# # write.table(scan1PermutedPvals, file='res~locus_scan1PermutedPvals.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
scan1_PermutedPvals <- read.table('res~locus_scan1PermutedPvals.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
#Remove Chr 2 since we are not doing linkage mapping on the chromsome due to the aneuploidy
valsSet <- scan1_PermutedPvals[geno$c != 2,]
# #Get the minimum p value from each genome wide (each column) scan off permuted growth values 
# getMins <- sapply(1:ncol(valsSet), function(y) {
	# min(valsSet[,y], na.rm=TRUE)
# })
# thresh1 <- -log10(minVals[order(minVals)][50]) # 4.305027
#Scan 1 Peaks
# pvals1 <- read.table('res~locus_scan1Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
# colnames(pvals1)[3] <- 'genoRow'
# scan1Peaks <- callPeaks(genoAll[genoAll$c != 17,], pvals1, threshScans[1], 2 )
# # write.table(scan1Peaks, file='scan1Peaks.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)


#Scan 2   Incorporate scan1 peaks, take residuals and run more genome wide linkage mapping

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
#first thing, correct for effect of causal MKT1 variant
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt)))
#scan1 peaks
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci1 <- scan1Peaks[,1:2]
lociCalls <- lapply(1:nrow(loci1), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci1$c[x] & geno$p == loci1$p[x], colnames(geno) %in% pheno$Sample]))
})
temp <- unique(unlist(lapply(1:length(lociCalls), function(x) {
	which(is.na(lociCalls[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop <- ncol(lociCalls) + 1
} else {
	toDrop <- temp
}
#remove individuals lacking a call at a given locus
lociCalls <-  lapply(1:nrow(loci1), function(x) {
	lociCalls[[x]][-toDrop]
})
res1 <- as.numeric(summary(lm(res[-toDrop] ~ lociCalls[[1]] + lociCalls[[2]] + lociCalls[[3]] + lociCalls[[4]] + lociCalls[[5]] + lociCalls[[6]] + lociCalls[[7]] ))$residuals) #try without chr 2 peak
#scan2 
scan2 <- do.call('rbind', lapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop]
	if (length(unique(genoRow)) == 1) {
		pvals <- c(NA)
	} else {
		pvals <- summary.aov(lm(res1 ~ as.factor(genoRow)))[[1]][,5][1]
	} 
	if (length(pvals) != 1) {
		pvals <- c(NA)
	} else {
		pvals
	}
}))
scan2 <- as.data.frame(scan2)
out <- cbind(geno[,1:2], scan2[,1])
write.table(out, file='res1of7loci~locus_scan2Pvals_New.txt', col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
#Scan 2 Permutations
# ####scaramble phenos 1000 for 1000 genome wide scans
# res1 <- as.numeric(summary(lm(res[-toDrop] ~ lociCalls[[1]] + lociCalls[[2]] + lociCalls[[3]] + lociCalls[[4]] + lociCalls[[5]] + lociCalls[[6]] + lociCalls[[7]] ))$residuals) #try without chr 2 peak
# permPhenos <- do.call('cbind', lapply(1:1000, function(x) {
	# sample(res1, length(res1), replace=TRUE)
# }))
# permP <- as.data.frame(permPhenos)
# write.table(permP, file='permutedRes1Phenos_1000.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
# permRes <- read.table('permutedRes1Phenos_1000.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)
# scan2PermutedPvals <- do.call('cbind', lapply(1:1000, function(y) {
	# print(paste('permutation genome scan', y))
	# permVals <- permRes[,y]
	# scan <- do.call('rbind', lapply(1:nrow(geno[geno$c != 17,]), function(x) {
		# genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop]
		# if (length(unique(genoRow)) == 1) {
			# c(NA)
		# } else {
			# pvals <- summary.aov(lm(permVals ~ genoRow))[[1]][,5][1]
		# } 
		# if (length(pvals) != 1) {
			# c(NA)
		# } else {
			# pvals
		# }
		
	# }))
# }))
# scan2PermutedPvals <- as.data.frame(scan2PermutedPvals)
# write.table(scan2PermutedPvals, file='res1of7loci~locus_scan1PermutedPvals.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
# scan2_PermutedPvals <- read.table('res1of7loci~locus_scan1PermutedPvals.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
#Remove Chr 2 since we are not doing linkage mapping on the chromsome due to the aneuploidy
# valsSet <- scan2_PermutedPvals[geno$c != 2,]
# #Get the minimum p value from each genome wide (each column) scan off permuted growth values 
# getMins <- sapply(1:ncol(valsSet), function(y) {
	# min(valsSet[,y], na.rm=TRUE)
# })
# thresh2 <- -log10(minVals[order(minVals)][50]) # 4.325596
#Scan 2 Peaks
# pvals2 <- read.table('res1of7loci~locus_scan2Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
# colnames(pvals2)[3] <- 'genoRow'
# scan2Peaks <- callPeaks(genoAll[genoAll$c != 17,], pvals2, threshScans[2], 2 )
# write.table(scan2Peaks, file='scan2Peaks.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)




#Scan 3   Incorporate scan2 peaks, take residuals and run more genome wide linkage mapping

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt)))
#scan1 peaks removed
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci1 <- scan1Peaks[,1:2]
lociCalls <- lapply(1:nrow(loci1), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci1$c[x] & geno$p == loci1$p[x], colnames(geno) %in% pheno$Sample]))
})
temp <- unique(unlist(lapply(1:length(lociCalls), function(x) {
	which(is.na(lociCalls[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop <- ncol(lociCalls) + 1
} else {
	toDrop <- temp
}
#remove individuals lacking a call at a given locus
lociCalls <-  lapply(1:nrow(loci1), function(x) {
	lociCalls[[x]][-toDrop]
})
res1 <- as.numeric(summary(lm(res[-toDrop] ~ lociCalls[[1]] + lociCalls[[2]] + lociCalls[[3]] + lociCalls[[4]] + lociCalls[[5]] + lociCalls[[6]] + lociCalls[[7]] ))$residuals) #try without chr 2 peak
#scan 2 peaks remvoved
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci2 <- scan2Peaks[,1:2]
inds <- pheno$Sample[-toDrop]
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci2$c[x] & geno$p == loci2$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls2), function(x) {
	which(is.na(lociCalls2[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop2 <- ncol(lociCalls2) + 1
} else {
	toDrop2 <- temp
}
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	lociCalls2[[x]][-toDrop2]
})
res2 <- as.numeric(summary(lm(res1[-toDrop2] ~ lociCalls2[[1]] + lociCalls2[[2]] + lociCalls2[[3]] + lociCalls2[[4]]) )$residuals)
#Scan3 
scan3 <- do.call('rbind', lapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop][-toDrop2]
	if (length(unique(genoRow)) == 1) {
		pvals <- c(NA)
	} else {
		pvals <- summary.aov(lm(res2 ~ as.factor(genoRow)))[[1]][,5][1]
	} 
	if (length(pvals) != 1) {
		pvals <- c(NA)
	} else {
		pvals
	}
}))
scan3 <- as.data.frame(scan3)
out <- cbind(geno[,1:2], scan3[,1])
write.table(out, file='res2of4loci~locus_scan3Pvals_New.txt', col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
#Scan 3 Permutations
# ####scaramble phenos 1000 for 1000 genome wide scans
# permPhenos <- do.call('cbind', lapply(1:1000, function(x) {
	# sample(res2, length(res2), replace=TRUE)
# }))
# permP <- as.data.frame(permPhenos)
# write.table(permP, file='permutedRes2Phenos_1000.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
# permRes <- read.table('permutedRes2Phenos_1000.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)
# scan3PermutedPvals <- do.call('cbind', lapply(1:1000, function(y) {
	# print(paste('permutation genome scan', y))
	# permVals <- permRes[,y]
	# scan <- do.call('rbind', lapply(1:nrow(geno[geno$c != 17,]), function(x) {
		# genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop][-toDrop2]
		# if (length(unique(genoRow)) == 1) {
			# c(NA)
		# } else {
			# pvals <- summary.aov(lm(permVals ~ genoRow))[[1]][,5][1]
		# } 
		# if (length(pvals) != 1) {
			# c(NA)
		# } else {
			# pvals
		# }
		
	# }))
# }))
# scan3PermutedPvals <- as.data.frame(scan3PermutedPvals)
# write.table(scan3PermutedPvals, file='res2of4loci~locus_scan3PermutedPvals.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
# scan3_PermutedPvals <- read.table('res2of4loci~locus_scan3PermutedPvals.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
#Remove Chr 2 since we are not doing linkage mapping on the chromsome due to the aneuploidy
# valsSet <- scan3_PermutedPvals[geno$c != 2,]
# #Get the minimum p value from each genome wide (each column) scan off permuted growth values 
# getMins <- sapply(1:ncol(valsSet), function(y) {
	# min(valsSet[,y], na.rm=TRUE)
# })
# thresh3 <- -log10(minVals[order(minVals)][50]) # 4.274288
#Scan 3 Peaks
# pvals3 <- read.table('res2of4loci~locus_scan3Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
# colnames(pvals3)[3] <- 'genoRow'
# scan3Peaks <- callPeaks(genoAll[genoAll$c != 17,], pvals3, threshScans[3], 2 )
# write.table(scan3Peaks, file='scan3Peaks.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)



#Scan4   Incorporate scan3 peaks, take residuals and run more genome wide linkage mapping

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt)))
#scan1 peaks removed
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci1 <- scan1Peaks[,1:2]
lociCalls <- lapply(1:nrow(loci1), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci1$c[x] & geno$p == loci1$p[x], colnames(geno) %in% pheno$Sample]))
})
temp <- unique(unlist(lapply(1:length(lociCalls), function(x) {
	which(is.na(lociCalls[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop <- ncol(lociCalls) + 1
} else {
	toDrop <- temp
}
#remove individuals lacking a call at a given locus
lociCalls <-  lapply(1:nrow(loci1), function(x) {
	lociCalls[[x]][-toDrop]
})
res1 <- as.numeric(summary(lm(res[-toDrop] ~ lociCalls[[1]] + lociCalls[[2]] + lociCalls[[3]] + lociCalls[[4]] + lociCalls[[5]] + lociCalls[[6]] + lociCalls[[7]] ))$residuals) #try without chr 2 peak
#scan 2 peaks remvoved
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci2 <- scan2Peaks[,1:2]
inds <- pheno$Sample[-toDrop]
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci2$c[x] & geno$p == loci2$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls2), function(x) {
	which(is.na(lociCalls2[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop2 <- ncol(lociCalls2) + 1
} else {
	toDrop2 <- temp
}
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	lociCalls2[[x]][-toDrop2]
})
res2 <- as.numeric(summary(lm(res1[-toDrop2] ~ lociCalls2[[1]] + lociCalls2[[2]] + lociCalls2[[3]] + lociCalls2[[4]]) )$residuals)
#scan3 peaks removed
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci3 <- scan3Peaks[,1:2]
inds <- pheno$Sample[-toDrop][-toDrop2]
lociCalls3 <- lapply(1:nrow(loci3), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci3$c[x] & geno$p == loci3$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls3), function(x) {
	which(is.na(lociCalls3[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop3 <- ncol(lociCalls3) + 1
} else {
	toDrop3 <- temp
}
lociCalls3 <- lapply(1:nrow(loci3), function(x) {
	lociCalls3[[x]][-toDrop3]
})
res3 <- as.numeric(summary(lm(res2[-toDrop3] ~ lociCalls3[[1]] + lociCalls3[[2]]))$residuals )
#Scan 4
scan4 <- do.call('rbind', lapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop][-toDrop2][-toDrop3]
	if (length(unique(genoRow)) == 1) {
		pvals <- c(NA)
	} else {
		pvals <- summary.aov(lm(res3 ~ as.factor(genoRow)))[[1]][,5][1]
	} 
	if (length(pvals) != 1) {
		pvals <- c(NA)
	} else {
		pvals
	}	
}))
scan4 <- as.data.frame(scan4)
out <- cbind(geno[,1:2], scan4[,1])
write.table(out, file='res3of2loci~locus_scan4Pvals_New.txt', col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
#Permutations scan 4
# ####scaramble phenos 1000 for 1000 genome wide scans
# permPhenos <- do.call('cbind', lapply(1:1000, function(x) {
	# sample(res3, length(res3), replace=TRUE)
# }))
# permP <- as.data.frame(permPhenos)
# write.table(permP, file='permutedRes3Phenos_1000.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
# permRes <- read.table('permutedRes3Phenos_1000.txt', sep='\t', header=FALSE, stringsAsFactors=FALSE)
# scan4PermutedPvals <- do.call('cbind', lapply(1:1000, function(y) {
	# print(paste('permutation genome scan', y))
	# permVals <- permRes[,y]
	# scan <- do.call('rbind', lapply(1:nrow(geno[geno$c != 17,]), function(x) {
		# genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop][-toDrop2][-toDrop3]
		# if (length(unique(genoRow)) == 1) {
			# c(NA)
		# } else {
			# pvals <- summary.aov(lm(permVals ~ genoRow))[[1]][,5][1]
		# } 
		# if (length(pvals) != 1) {
			# c(NA)
		# } else {
			# pvals
		# }
		
	# }))
# }))
# scan4PermutedPvals <- as.data.frame(scan4PermutedPvals)
# write.table(scan4PermutedPvals, file='res3of2loci~locus_scan4PermutedPvals.txt', quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
# scan4_PermutedPvals <- read.table('res3of2loci~locus_scan4PermutedPvals.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
#Remove Chr 2 since we are not doing linkage mapping on the chromsome due to the aneuploidy
# valsSet <- scan4_PermutedPvals[geno$c != 2,]
# #Get the minimum p value from each genome wide (each column) scan off permuted growth values 
# getMins <- sapply(1:ncol(valsSet), function(y) {
	# min(valsSet[,y], na.rm=TRUE)
# })
# thresh4 <- -log10(minVals[order(minVals)][50]) # 4.332096
#Scan 4 Peaks
pvals4 <- read.table('res3of2loci~locus_scan4Pvals_New.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
colnames(pvals4)[3] <- 'genoRow'
scan4Peaks <- callPeaks(genoAll[genoAll$c != 17,], pvals4, threshScans[4], 2 )
write.table(scan4Peaks, file='scan4Peaks.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)



# Scan 5 

pheno <- phenoAll[grepl('D', phenoAll$Sample) | grepl('E', phenoAll$Sample),] #crosses D and E are mrp20-105E 
geno <- cbind(genoAll[, colnames(genoAll) %in% c('c', 'p')], genoAll[,colnames(genoAll) %in% pheno$Sample])
mkt <- as.factor(as.numeric(grepl('E', pheno$Sample))) #match to 3S allele yields a 1, match haplotype encoding
res <- as.numeric(residuals(lm(pheno$Ethanol ~ mkt)))
#scan1 peaks removed
scan1Peaks <- read.table('scan1Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci1 <- scan1Peaks[,1:2]
lociCalls <- lapply(1:nrow(loci1), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci1$c[x] & geno$p == loci1$p[x], colnames(geno) %in% pheno$Sample]))
})
temp <- unique(unlist(lapply(1:length(lociCalls), function(x) {
	which(is.na(lociCalls[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop <- ncol(lociCalls) + 1
} else {
	toDrop <- temp
}
#remove individuals lacking a call at a given locus
lociCalls <-  lapply(1:nrow(loci1), function(x) {
	lociCalls[[x]][-toDrop]
})
res1 <- as.numeric(summary(lm(res[-toDrop] ~ lociCalls[[1]] + lociCalls[[2]] + lociCalls[[3]] + lociCalls[[4]] + lociCalls[[5]] + lociCalls[[6]] + lociCalls[[7]] ))$residuals) #try without chr 2 peak
#scan 2 peaks remvoved
scan2Peaks <- read.table('scan2Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci2 <- scan2Peaks[,1:2]
inds <- pheno$Sample[-toDrop]
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci2$c[x] & geno$p == loci2$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls2), function(x) {
	which(is.na(lociCalls2[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop2 <- ncol(lociCalls2) + 1
} else {
	toDrop2 <- temp
}
lociCalls2 <- lapply(1:nrow(loci2), function(x) {
	lociCalls2[[x]][-toDrop2]
})
res2 <- as.numeric(summary(lm(res1[-toDrop2] ~ lociCalls2[[1]] + lociCalls2[[2]] + lociCalls2[[3]] + lociCalls2[[4]]) )$residuals)
#scan3 peaks removed
scan3Peaks <- read.table('scan3Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci3 <- scan3Peaks[,1:2]
inds <- pheno$Sample[-toDrop][-toDrop2]
lociCalls3 <- lapply(1:nrow(loci3), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci3$c[x] & geno$p == loci3$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls3), function(x) {
	which(is.na(lociCalls3[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop3 <- ncol(lociCalls3) + 1
} else {
	toDrop3 <- temp
}
lociCalls3 <- lapply(1:nrow(loci3), function(x) {
	lociCalls3[[x]][-toDrop3]
})
res3 <- as.numeric(summary(lm(res2[-toDrop3] ~ lociCalls3[[1]] + lociCalls3[[2]]))$residuals )
# scan 4 peaks removed

scan4Peaks <- read.table('scan4Peaks.txt', header=TRUE, stringsAsFactors=FALSE, sep='\t')
loci4 <- scan4Peaks[,1:2]
inds <- pheno$Sample[-toDrop][-toDrop2][-toDrop3]
lociCalls4 <- lapply(1:nrow(loci4), function(x) {
	call <- as.character(as.numeric(geno[geno$c == loci4$c[x] & geno$p == loci4$p[x], colnames(geno) %in% inds]))
})
temp <- unique(unlist(lapply(1:length(lociCalls4), function(x) {
	which(is.na(lociCalls4[[x]]))
})))
temp <- temp[order(temp)]
if (length(temp) == 0) {
	toDrop4 <- ncol(lociCalls4) + 1
} else {
	toDrop4 <- temp
}
lociCalls4 <- lapply(1:nrow(loci4), function(x) {
	lociCalls4[[x]][-toDrop4]
})
res4 <- as.numeric(summary(lm(res3[-toDrop4] ~ lociCalls4[[1]] + lociCalls4[[2]] + lociCalls4[[3]] ))$residuals)
# Scan 5
scan5 <- do.call('rbind', lapply(1:nrow(geno), function(x) {
	print(x)
	genoRow <- as.character(as.numeric(geno[x, colnames(geno) %in% pheno$Sample]))[-toDrop][-toDrop2][-toDrop3][-toDrop4]
	if (length(unique(genoRow)) == 1) {
		pvals <- c(NA)
	} else {
		pvals <- summary.aov(lm(res4 ~ as.factor(genoRow)))[[1]][,5][1]
	} 
	if (length(pvals) != 1) {
		pvals <- c(NA)
	} else {
		pvals
	}	
}))
scan5 <- as.data.frame(scan5)
out <- cbind(geno[,1:2], scan5[,1])
write.table(out, file='res4of3loci~locus_scan5Pvals_New.txt', col.names=TRUE, quote=FALSE, row.names=FALSE, sep='\t')
#min val is not near previous permuted threhsolds of ~ lod 4.2-4.3, so don't run permutations
# -log10(min(out[,3][geno$c != 2], na.rm=TRUE))  # 3.681641