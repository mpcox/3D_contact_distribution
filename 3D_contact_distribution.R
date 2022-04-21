#!/usr/bin/env Rscript


# 3D contact distribution


# size of sparse matrix strongly affects runtime
# consider reducing resolution (increasing bin sizes)


# libraries
require(data.table)


# global variables
title              <- "RXLRs"
all_gene_file      <- "example_all_genes.gff3"
subset_gene_file   <- "example_gene_subset.gff3"
sparse_matrix_file <- "example_matrix.dat"
bin.size           <- 10000
metric             <- "mean"   # default is "mean"; alternative is "median"
iterations         <- 100      # number of iterations in Monte Carlo simulations
max.number         <- 10000    # maximum number of pairwise comparisons to consider (otherwise exponential)


# functions
# functions
# find number of contacts between two genes
find_contacts <- function(contact.matrix, gene1.chr, gene1.start, gene1.stop, gene2.chr, gene2.start, gene2.stop, metric="mean"){
	
	# examples
	# find_contacts(m, 1, 27526, 32376, 1, 59914, 60699)
	# find_contacts(m, 10, 71000, 75000, 10, 21000, 26000)

	# calculate matrix bins that genes fall within
	gene1.bin1 <- gene1.start - (gene1.start %% bin.size)
	gene1.bin2 <- gene1.stop - (gene1.stop %% bin.size)
	
	gene2.bin1 <- gene2.start - (gene2.start %% bin.size)
	gene2.bin2 <- gene2.stop - (gene2.stop %% bin.size)
	
	# generate empty contacts vector
	contacts <- vector()
	
	# calculate contacts for all pairs of bins (including absent values given the matrix is sparse)
	val1 <- contact.matrix[.(gene1.chr, gene1.bin1, gene2.chr, gene2.bin1)]$contacts
	if( is.na(val1) ){ contacts <- 0 } else { contacts <- val1 }
	
	if( gene2.bin2 != gene2.bin1 ){
		val2 <- contact.matrix[.(gene1.chr, gene1.bin1, gene2.chr, gene2.bin2)]$contacts
		if( is.na(val2) ){ contacts <- c(contacts, 0) } else { contacts <- c(contacts, val2) }
	}
	
	if( gene1.bin2 != gene1.bin1 ){
		val3 <- contact.matrix[.(gene1.chr, gene1.bin2, gene2.chr, gene2.bin1)]$contacts
		if( is.na(val3) ){ contacts <- c(contacts, 0) } else { contacts <- c(contacts, val3) }
	}
	
	if( gene1.bin2 != gene1.bin1 & gene2.bin2 != gene2.bin1 ){
		val4 <- contact.matrix[.(gene1.chr, gene1.bin2, gene2.chr, gene2.bin2)]$contacts
		if( is.na(val4) ){ contacts <- c(contacts, 0) } else { contacts <- c(contacts, val4) }
	}
	
	# return mean or median value of contacts
	if( metric == "mean" ){
		return( mean(contacts) )
	} else {
		return( median(contacts) )
	}
}

# find mean or median number of contacts between all (or a reduced set of) non-self pairs
# note: function is *slow*
#       profiling shows ~99% of runtime is identifying the required row in the sparse matrix
#       in data.table, this function is already heavily optimized, so it is hard to improve this
#       the only feasible option seems to be to reduce the size of the matrix, by using larger window bins
summary_contacts <- function(gene.dataset, contact.matrix, max.number, metric="mean"){

	# calculate all gene set pairs
	combinatorics <- t(combn(length(gene.dataset[,1]), 2))
	
	# downsample if number larger than user-set threshold
	if( length(combinatorics[,1]) > max.number ){
		rows <- sort(sample(1:length(combinatorics[,1]), max.number, replace=F))
		combinatorics <- combinatorics[rows,]
	}
	combinatorics <- as.data.frame(combinatorics)
	colnames(combinatorics) <- c("first", "second")
	
	# calculate pairwise contacts
	# contacts <- mapply(find_contacts,
	#        contact.matrix = contact.matrix,
	#        gene1.chr = gene.dataset$chr[combinatorics$first], 
	#        gene1.start = gene.dataset$start[combinatorics$first],
	#        gene1.stop = gene.dataset$stop[combinatorics$first],
	#        gene2.chr = gene.dataset$chr[combinatorics$second],
	#        gene2.start = gene.dataset$start[combinatorics$second],
	#        gene2.stop = gene.dataset$stop[combinatorics$second])
	
	# calculate pairwise contacts
	contacts <- vector(length=length(combinatorics$first))
	
	for( i in 1:length(combinatorics$first) ){
		contacts[i] <- find_contacts(contact.matrix = contact.matrix,
									 gene1.chr = gene.dataset$chr[combinatorics$first][i], 
									 gene1.start = gene.dataset$start[combinatorics$first][i],
									 gene1.stop = gene.dataset$stop[combinatorics$first][i],
									 gene2.chr = gene.dataset$chr[combinatorics$second][i],
									 gene2.start = gene.dataset$start[combinatorics$second][i],
									 gene2.stop = gene.dataset$stop[combinatorics$second][i],
									 metric)
	}
	
	# calculate mean or median number of contacts
	if( metric == "mean" ){
		return( mean(contacts) )
	} else {
		return( median(contacts) )
	}
}

# monte carlo simulator
monte_carlo <- function(allgene.dataframe, sparse.matrix, sample.size, max.number, iterations, observed.summary, metric="mean"){
	
	# exclude mitochondrial genes
	mt.genes <- which(allgene.dataframe$chr == 'mt')
	allgene.dataframe <- allgene.dataframe[-mt.genes,]
	
	# make distribution vector
	summaries <- vector(length=iterations)
	
	# run iterations
	for (i in 1:iterations) {
		
		sim <- allgene.dataframe[sample(nrow(allgene.dataframe), sample.size), ]
		sim$chr <- as.integer(as.character(sim$chr))
		
		summaries[i] <- summary_contacts(sim, sparse.matrix, max.number, metric)
	}
	
	# calculate how many simulated datasets have median values larger/smaller than the observed number
	calc_prop(summaries, observed.summary)
	
	return( summaries )
}

# calculate probabilities
calc_prop <- function(vector, observed_value){
	
	n_greater <- length(which(vector >= observed_value))
	n_less    <- length(which(vector <= observed_value))
	
	message(paste('Greater than', observed_value, ': p =', n_greater / length(vector)))
	message(paste('Less than', observed_value, ': p =', n_less / length(vector)))
}


# data acquisition
# extract genes only from the GFF3 file (all genes)
command <- paste0("grep \"gene\" ", all_gene_file, " > all_genes.gff3")
cat(command,"\n")
try(system(command))

# read all genes GFF3
d <- read.table("all_genes.gff3", header=F)
colnames(d) <- c("chr", "software", "class", "start", "stop", "null1", "strand", "null2", "comments")
head(d)


# extract genes only from the GFF3 file (subset genes)
command <- paste0("grep \"gene\" ", subset_gene_file, " > subset_genes.gff3")
cat(command,"\n")
try(system(command))

# read subset genes GFF3
s <- read.table("subset_genes.gff3", header=F)
colnames(s) <- c("chr", "software", "class", "start", "stop", "null1", "strand", "null2", "comments")
head(s)


# read sparse matrix and change to indexed data.table
m <- read.table(sparse_matrix_file, header=F)
colnames(m) <- c("first.chr", "second.chr", "first.start", "second.start", "contacts")
m <- data.table(m)
setkey(m, first.chr, first.start, second.chr, second.start, contacts)
head(m)


# calculate mean or median number of contacts between subset genes
start_time <- Sys.time()
summary.subset <- summary_contacts(s, m, max.number, metric=metric)
end_time <- Sys.time()
end_time - start_time
summary.subset


# run monte carlo
start_time <- Sys.time()
summary.values <- monte_carlo(d, m, length(s$chr), max.number, iterations, summary.subset, metric=metric)
end_time <- Sys.time()
end_time - start_time

