library(ape)
library(Biostrings)

#############################
# SETTING UP THE INPUT DATA #
#############################

# Load the gff. Put the path to your .gff3 file here. 
gff.path <- 'Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
gff <- read.gff(gff.path, GFF3 = T)

# Load the genome. Put the path to your genome here. 
genome.path <- 'Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa'
genome <- readDNAStringSet(genome.path)
names(genome) <- sub(' .*', '', names(genome))


############################
# SETTING UP THE FUNCTIONS # 
############################

# The user doesn't need to do anything in this section. We're just creating the functions that we will need. 

# Let's write a function to extract useful information from a gene of interest 
extract.gene.info <- function(agap.number, gff.table, isoform = 'A'){
	if (!(isoform %in% LETTERS))
		stop('"isoform" should be a single capitalised alphabetical letter')
	isoform.num <- which(LETTERS == isoform)
	gff.focus <- subset(gff.table, grepl(agap.number, attributes) & (!grepl(paste(agap.number, '-[PR][', paste(LETTERS[-isoform.num], collapse = ''), ']', sep = ''), attributes)))
	# Get the vector of indices that will reconstitute the gene sequence
	cds.ranges <- subset(gff.focus, type == 'CDS')[, c('start', 'end')]
	# Get the strand direction for the gene and do a sanity check that only one strand is obtained
	strand <- unique(as.character(gff.focus$strand))
	if (length(strand) > 1) stop('More than one strand direction detected')
	# Same for chromosome
	chrom <- unique(as.character(gff.focus$seqid))
	if (length(chrom) != 1) stop('The gene should be found on one and only one chromosome.')
	# Get the sequence for each of the CDS segments. 
	sequence.list <- apply(cds.ranges, 1, function(x) substr(genome[[chrom]], x[1], x[2]))
	# c() can be used to paste the sequences together, but this doesn't seem to work with do.call (which just
	# concatenates everything like the normal c() behaviour), so I don't know how to combine everything together 
	# except with this clusmy for loop. 
	full.cds.sequence <- sequence.list[[1]]
	for (s in sequence.list[2:length(sequence.list)]){
		full.cds.sequence <- c(full.cds.sequence, s)
	}
	# If necessary, reverse commplement the sequence
	if (strand == '-') 
		full.cds.sequence <- reverseComplement(full.cds.sequence)
	# Translate the sequence
	full.aa.sequence <- translate(full.cds.sequence)
	# Split the cds into codons. 
	full.codon.list <- sapply(seq(1, length(full.cds.sequence), by = 3), function(i) substr(full.cds.sequence, i, i+2))
	# Translate each of those codons
	full.aa.list <- sapply(full.codon.list, translate, no.init.codon = T)
	# And use that to translate the sequence. This is redundant as we've already done it, but it's just a sanity
	# check.
	full.aa.sequence.2 <- do.call(c, full.aa.list)
	if (full.aa.sequence != full.aa.sequence.2) stop('Codon-by-codon translation didn\'t match full sequence translation.')
	
	# For a given SNP genomic position, we need to work out in which part of the sequence it is and, if it's a CDS,
	# where on the CDS it sits
	gene.range <- c(min(gff.focus$start), max(gff.focus$end))
	sequence.indices <- rep('intron', gene.range[2] - gene.range[1] + 1)
	for (i in 1:nrow(gff.focus)){
		feature.type <- as.character(gff.focus[i, 'type'])
		if (feature.type %in% c('five_prime_UTR', 'CDS', 'three_prime_UTR'))
			sequence.indices[(gff.focus[i, 'start']:gff.focus[i, 'end']) - gene.range[1] + 1] <- feature.type
	}
	names(sequence.indices) <- as.character(gene.range[1]:gene.range[2])
	
	list(agap = agap.number,
	     sequence.indices = sequence.indices, 
	     cds.ranges = cds.ranges,
	     strand = strand, 
		 # To save on memory, store the CDS as character, not DNAStringSet
	     full.codon.list = sapply(full.codon.list, as.character), 
	     full.aa.list = sapply(full.aa.list, as.character))
}

# Write a function that will assess the relevance of a SNP, given its genomic position
assess.SNP <- function(snp, gene.data){
	snp.id <- strsplit(snp, ':')[[1]]
	genome.position <- snp.id[2]
	element.type <- gene.data$sequence.indices[genome.position]
	if (is.na(element.type))
		output <- 'outside'
	else if (element.type != 'CDS')
		output <- as.character(element.type)
	else{
		# Which nucleotide of the sequence is it from
		cds.indices <- unlist(apply(gene.data$cds.ranges, 1, function(x) x[1]:x[2]))
		if (gene.data$strand == '-')
			cds.indices <- rev(cds.indices)
		nuc.index <- which(cds.indices == genome.position)
		codon.index <- ((as.numeric(nuc.index) - 1) %/% 3) + 1
		within.codon.position <- ((as.numeric(nuc.index) - 1) %% 3) + 1
		ref.codon <- gene.data$full.codon.list[[codon.index]]
		ref.aa <- gene.data$full.aa.list[[codon.index]]
		ac <- gene.data$allele.counts[(genome.position), paste('counts', 0:3, sep = '')]
		if (length(snp.id) > 2) {
			alleles = strsplit(snp.id[3], '/')[[1]]
			# If the gene is on the negative strand, need to adjust the coding sequence alleles
			if (gene.data$strand == '-'){
				revcomp.table <- setNames(c('A','C','T','G'), c('T','G','A','C'))
				ref.allele <- revcomp.table[alleles[1]]
				mut.allele <- revcomp.table[alleles[2]]
			}
			else {
				ref.allele <- alleles[1]
				mut.allele <- alleles[2]
			}
			# Sanity check with ref allel
			ref.allele.2 <- substr(ref.codon, within.codon.position, within.codon.position)
			if (ref.allele != ref.allele.2) stop('Reference allele sanity check failed.')
			# Check whether the mutation results in a change in the amino acid
			mut.codon <- replaceAt(DNAString(ref.codon), IRanges(within.codon.position, width = 1), mut.allele)
			mut.aa <- translate(mut.codon, no.init.codon = T)
			if (as.character(mut.aa) == ref.aa)
				output <- paste('SNP(syn) in codon ', codon.index, ' position ', within.codon.position, sep = '')
			else 
				output <- paste('SNP(', ref.aa, '->', mut.aa, ') in codon ', codon.index, ' position ', within.codon.position, sep = '')
		}
		else
			output <- paste('Locus is in codon ', codon.index, ' position ', within.codon.position, sep = '')
	}
	print(paste(snp, output, sep = '   '))
	invisible(output)
}

# If you want to assess more than one SNP at once in a given gene, you can use this function instead.
assess.SNPs <- function(snp.ids, gene.data){
	invisible(sapply(snp.ids, assess.SNP, gene.data))
}

##################################
# EXAMPLE USAGE OF THE FUNCTIONS #
##################################

# For any gene that we want to look for SNPs in, we first extract the information for that gene from the gff, 
# using the gene's AGAP number. For example, let's say we want to look into SNPs from kdr (AGAP004707) and 
# Ace1 (AGAP001356). 
kdr.gene.info <- extract.gene.info('AGAP004707', gff)
ace1.gene.info <- extract.gene.info('AGAP001356', gff)

# Next we get the loci that were are interested in. Let's look at the following three SNPs. 
#kdr-995F	2L	2422652	[A/T]	
#VGC-1874S	2L	2430880	[C/T]	
#Ace1-280S	2R	3492074	[G/A]	
# For each SNP, we want a string in the format, CHROMOSOME:POSITION:REF/ALT.
# So, for example, if there is a SNP at position 2422652 on chromosome 2L, where the reference allele is A and
# the mutant allele in T, then we code this as 2L:2422652:A/T
kdr.snps.of.interest <- c('2L:2422652:A/T', '2L:2430880:C/T')
ace1.snps.of.interest <- c('2R:3492074:G/A')

# Then we run the assess.SNPs function (if there is only one SNP of interest, such as in the case of Ace1 here,
# you can also use the assess.SNP function instead). 
assess.SNPs(kdr.snps.of.interest, kdr.gene.info)
assess.SNPs(ace1.snps.of.interest, ace1.gene.info)

# You can see that the output tells you that there is a SNP, followed by the resulting amino-acid change, which 
# codon the SNP is in, and in which codon position. 

# By default, the extract.gene.info function uses isoform A for whatever gene you give it. But you can also tell
# it which isoform you want instead. For example, in kdr, the isoform used to determine SNP nomenclature is 
# isoform D (hence why in the output above, the 995F allele was reported in position 980, because it was using
# isoform A). We can rerun this with isoform D and see that it now reports it as codon 995:
kdr.D.info <- extract.gene.info('AGAP004707', gff, isoform = 'D')
assess.SNPs(kdr.snps.of.interest, kdr.D.info)

# If we don't know the reference / mutant alleles for the SNP in question, we can just code the SNP as 
# CHROMOSOME:POSITION, and the function will just report codon number and position:
kdr.snps.of.interest.noalleles <- c('2L:2422652', '2L:2430880')
assess.SNPs(kdr.snps.of.interest.noalleles, kdr.D.info)

# If we ask about a synonymous mutation, the function will report that it as synonymous
kdr.snps.non.synonymous <- c('2L:2422652:A/G', '2L:2430882:T/A')
assess.SNPs(kdr.snps.non.synonymous, kdr.D.info)

# If we give it mutations in non-coding regions, the function will report the type of mutation it is
ace1.snps.non.coding <- c('2R:3484116:G/A', '2R:3491479:C/A')
assess.SNPs(ace1.snps.non.coding, ace1.gene.info)

# And finally, if we give it a SNP that is outside of the region of the gene entirely, the function just reports
# "outside", like here if we run the function with the kdr SNPs but using the Ace1 gene information.
assess.SNPs(kdr.snps.of.interest, ace1.gene.info)




