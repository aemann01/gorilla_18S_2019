####################################
#Dada2 pipeline for Gorilla Samples
####################################



####Libraries####
library(dada2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(seqinr)
library(ShortRead)
library(Biostrings)

####Environment Setup####
theme_set(theme_bw())
setwd("~/Desktop/lab_member_files/morien/blastocystis_colitis_1_2019_18s/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/Desktop/lab_member_files/morien/blastocystis_colitis_1_2019_18s/raw_data/demultiplexed/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE)) #change the pattern to match all your R1 files
sample.names <- sapply(strsplit(basename(fnFs), "fastq.gz_"), `[`, 2) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name
sample.names <- gsub(".fastq.gz", "", sample.names)

####fastq Quality Plots####
plotQualityProfile(fnFs[1:20]) #quality variable, but most samples look good

####Primer Removal####
####identify primers####
FWD <- "CCAGCASCYGCGGTAATTCC"  ## CHANGE ME to your forward primer sequence
REV <- "ACTTTCGTTCTTGATYRA"  ## CHANGE ME to your reverse primer sequence
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]),  #add the index of the sample you'd like to use for this test (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]))

####OPTIONAL!!!!####
#REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details. it is linked at the top of this guide.

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[5]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[5]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
#150 should be well below the lower bound for V4 data
#if you are working with V9 data, I have found that a minLen of 80bp is appropriate. Giardia sequences are ~95bp in V9
out <- filterAndTrim(cutFs, filtFs,truncLen=c(0), minLen = c(150),
                     maxN=c(0), maxEE=c(8), truncQ=c(2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #error rates a bit weird, but probably okay

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs[[1]]

####construct sequence table####
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 1)) #how many are there with N reads or fewer? (N=1 in example)

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.18s_R1.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.18s.R1.txt", row.names=FALSE, quote=F, sep="\t")

#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.18s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

####assign taxonomy####
#note, this takes ages if you have a large dataset. strongly recommend doing on a multi-core machine (zoology cluster, or entamoeba in the lab). another option: saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function (their implementation of RDP).
#for 18s mammalian gut experiments
taxa_18s <- assignTaxonomy(seqtab.nosingletons.nochim, "~/Desktop/lab_member_files/taxonomy_databases/silva_for_dada2/v128_for_parfreylab/18s/silva_128.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa_18s[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa_18s[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa_18s[NAs,1] <- "Unassigned" #apply new label to identified indices
#set column ranks for the taxa table. should be Rank1 through RankN, depending on how many ranks you have. last column is the accession for our parfreylab custom-formatted databases, and species for other publicly available databases.
colnames(taxa_18s) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Rank8", "Accession") #for parfreylab version of SIVLA 128

####saving taxonomy data####
write.table(data.frame("row_names"=rownames(taxa_18s),taxa_18s),"taxonomy_table.18s_R1.txt", row.names=FALSE, quote=F, sep="\t")

####remove high-confidence prokaryotic taxa####
taxa_dada2 <- assignTaxonomy(seqtab.nosingletons.nochim, "~/Desktop/lab_member_files/taxonomy_databases/silva_for_dada2/v132/silva_nr_v132_train_set.fa.gz", multithread=TRUE) #for filtering out bacterial and archeal sequences

a <- which(taxa_dada2[,1] == "Archaea") #which have a bacterial or archaeal assignment
b <- which(taxa_dada2[,1] == "Bacteria")
c <- union(a,b)
d <- which((is.na(taxa_dada2[,2]) == F)) #which taxa have an assignment at phylum level (we can be more confident that these are really in the kingdoms they say they are in)
e <- intersect(c,d)

f <- which((is.na(taxa_18s[,2]) == T)) #which taxa do not have a phylum assignment from SILVA
bacterial_taxa <- Reduce(intersect, list(e,f)) #which taxa share these two conditions, i think we can safely remove them #i blasted to check these for eukaryotic hits, there were none

taxa_18s <- taxa_18s[-bacterial_taxa,]

#save taxonomy table again
write.table(data.frame("row_names"=rownames(taxa_18s),taxa_18s),"taxonomy_table.18s_R1.noprok.txt", row.names=FALSE, quote=F, sep="\t")

#### hand off to PhyloSeq ####
#load sample data
rawmetadata <- read_delim(file = file.path("~/Desktop/lab_member_files/morien/blastocystis_colitis_1_2019_18s/", "metadata_rat-blastocystis-colitis1_info-plate-TareukREV3-TareukFWD454.barcodes_added.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

#set row names as the sample IDs for the metadata
rownames(rawmetadata) <- rawmetadata$`#SampleID`

#record samples absent in either metadata or OTU table
notinmeta <- setdiff(row.names(seqtab.nosingletons.nochim), row.names(rawmetadata))
notinraw <- setdiff(row.names(rawmetadata), row.names(seqtab.nosingletons.nochim))

#create phyloseq object from "seqtab.nosingletons.nochim", "rawmetadata", and "taxa"
ps.dada2_join <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(taxa_18s))

# now replace the long ASV names (the actual sequences) with human-readable names, and save the new names and sequences as a .fasta file in your project working directory
my_otu_table <- t(as.data.frame(unclass(otu_table(ps.dada2_join)))) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ntaxa(ps.dada2_join)), sep='') #create new names
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "blastocystis_colitis_1_2019.18s_ASV_sequences.fasta") #save sequences with new names in fasta format
taxa_names(ps.dada2_join) <- ASV.num #rename your sequences in the phyloseq object

# at this point you i recommend saving your complete experiment as a phyloseq object, for ease of picking up here if you wish to make changes to your filtering criteria later on
saveRDS(ps.dada2_join, "blastocystis_colitis_1_2019.18s.full_dataset.phyloseq_format.RDS")

# Remove samples with less than N reads (N=100 in example. adjust per experiment)
ps.dada2_join <- prune_samples(sample_sums(ps.dada2_join) >= 100, ps.dada2_join)

#OPTIONAL: Remove OTUs with less than N total reads. (N = 50 for example. adjust per experiment) 
ps.dada2_join <- prune_taxa(taxa_sums(ps.dada2_join) >= 50, ps.dada2_join)

# Set aside unassigned taxa #with dada2, there may not be any unassigned taxa as dada2's RDP classifier usually ends up classifying everything. You may want to adjust this command to set aside ASVs that have no assignment beyond Rank1 instead of no assignment at Rank1.
ps.dada2_join.unassigned <- ps.dada2_join %>%
  subset_taxa(Rank1 == "Unassigned")
# Remove unassigned taxa
ps.dada2_join <- ps.dada2_join %>%
  subset_taxa(Rank1 != "Unassigned")

# Remove counts of 1 from OTU table
otu <- as.data.frame(otu_table(ps.dada2_join))
otu_table(ps.dada2_join)[otu <= 1] <- 0

#after denoising you can also remove ASVs that are in groups you wish to exclude (i.e. mammalia, embryophyta, etc.)
#to do this, just determine which rank the clade is captured by, and filter like so:
# Remove mammalian and plant OTUs #this is just an example for a human gut study where we're not interested in anything but what's living in the gut
ps.dada2_join <- ps.dada2_join %>%
  subset_taxa(Rank5 != "Mammalia") %>% #you can chain as many of these subset_taxa calls as you like into this single command using the R pipe (%>%)
  subset_taxa(Rank5 != "Embryophyta")

#with your filtered phyloseq object, you are now ready to move on to whatever statistical/diversity analyses you're planning to do. please see other R script guides in the lab github.
