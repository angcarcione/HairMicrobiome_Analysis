#urls from SRA for project PRJDB8176
urls <- c(
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175005/DRR175005_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175005/DRR175005_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175019/DRR175019_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175019/DRR175019_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175006/DRR175006_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175006/DRR175006_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175018/DRR175018_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175018/DRR175018_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175007/DRR175007_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175007/DRR175007_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175011/DRR175011_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175011/DRR175011_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175013/DRR175013_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175013/DRR175013_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175017/DRR175017_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175017/DRR175017_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175012/DRR175012_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175012/DRR175012_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175008/DRR175008_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175008/DRR175008_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175009/DRR175009_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175009/DRR175009_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175016/DRR175016_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175016/DRR175016_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175015/DRR175015_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175015/DRR175015_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175010/DRR175010_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175010/DRR175010_2.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175014/DRR175014_1.fastq.gz",
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR175/DRR175014/DRR175014_2.fastq.gz"
)

#function that "w-gets" the urls so i dont have to do it one by one
download_files <- function(urls) {
  for (url in urls) {
    destfile <- basename(url)
    system(paste("wget -O", destfile, url))
  }
}

download_files(urls)

#############################Primer Trimming 
## trimmed one by one lol pain THIS IS WRONG IGRNORE

# 515f, 5′-GTG CCA GCM GCC GCG GTA A-3′
# 806r, 5′-GGA CTA CHV GGG TWT CTA AT-3′
# 
# 1-515F, 5′-TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG GTG CCA GCM GCC GCG GTA A-3′ 
# 1-806R, 5′-GTC TCG TGG GCT CGG AGA TGT GTA TAA GAG ACA GGG ACT ACH VGG GTW TCT AAT-3
# 
# Forward primer, 5′-AAT GAT ACG GCG ACC ACC GAG ATC TAC AC-Index sequence-TCG TCG GCA GCG TC-3′
# Reverse primer 5′-CAA GCA GAA GAC GGC ATA CGA GAT-Index sequence-GTC TCG TGG GCT CGG-3′
# 
# # Loop through each file ending in _1.fastq
# foreach ($read1 in Get-ChildItem -Filter "*_1.fastq") {
#   # Define the corresponding _2 file
#   $read2 = $read1.Name -replace "_1.fastq", "_2.fastq"
#   
#   # Define the output file names
#   $out1 = $read1.Name -replace "_1", "_1.trimmed"
#   $out2 = $read2 -replace "_2", "_2.trimmed"
#   
#   # Run cutadapt command
#   cutadapt -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -o $out1 -p $out2 $read1.Name $read2 | Tee-Object -FilePath "cutadapt.log.txt" -Append
# }




############################
#check my directory so im not working in the wrong folder
system('pwd')
system('ls')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

library(dada2); packageVersion("dada2")

path <- "/cloud/project/HairReads"
######
#sorting my trimmed fastq files that i trimmed using cut adapt
fnFs <- sort(list.files(path, pattern="_1.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.trimmed.fastq", full.names = TRUE))

#extracting the sample names 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#checking for quality of all of my forward (fnFs) and reverse (fnRs) reads
plotQualityProfile(fnFs[1:15])
plotQualityProfile(fnRs[1:15])

#putting all of the filtered files into their own directory 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering parameters 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
#On Windows set multithread=FALSE (i have a windows but also this is an environment so it doesn't matter)
head(out)

#learning the error for forward and reverse reads and then plotting the estimated error rates for transitions from nucleotide to nucleotide 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#applying "core sample inference algorithm" to our filtered and trimmed reads to find the number of unique reads in sequences
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#shows us the denoising results for the forward and reverse reads 
dadaFs[[1]]
dadaRs[[1]]

####################Merging paired reads

#merge the forward and the reverse reads now that they are denoised 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#creating my ASV table and looking at how many samples vs original ASVs 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#looking at how many sequences each sample has 
table(nchar(getSequences(seqtab)))

##chimera check and removing them
#identifies chimeras and removes them. i had 180 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#tells me what percent of the chimeras were removed 
sum(seqtab.nochim)/sum(seqtab)

#tracks made through pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#allows us to visualize how many reads made it through each step of the pipeline
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#just shows a few
head(track)
#shows all 
track
##############################################Naive Bayesian Classifier method for Assigning taxonomy 
list.files("/cloud/project/HairReads")

#assigning taxonomy using the most up to date training det silva nr99 v138.1
taxa <- assignTaxonomy(seqtab.nochim, "/cloud/project/HairReads/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#assinging species with the most up to date species classification of v138.1
taxa <- addSpecies(taxa, "/cloud/project/HairReads/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#### no mock! so i cant evaluate accuracy 

########################################################PHYLOSEQ
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq")


BiocManager::install("Biostrings")
library(Biostrings)
packageVersion("Biostrings")

BiocManager::install("ggplot2")
library(ggplot2)
packageVersion("ggplot2")

#Extracting the sample names and subjects
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

#Reading in the metadata 
sample_data <- read.csv("/cloud/project/HairReads/SraRunTable_hair.csv", row.names = 1)
sample_data

#Create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(sample_data), 
               tax_table(taxa))
#didnt need to remove mock sample
print(ps)

#rename samples with ASV so its easier to interact with 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Check column names in sample data
print(colnames(sample_data(ps)))

#interested in Sample_name!!

##################################################SAMPLE STATISTICS FOR FINAL 
#reads per sample
reads_per_sample <- sample_sums(ps)
reads_per_sample


#visualizing reads per sample with a bar plot
barplot(reads_per_sample, las = 2, main = "Total Reads per Sample", ylab = "Read Count", xlab = "Samples")

##################################################
#adding a new ps variable "first word" so that HairRoot and other names arent repeated 3 times like it had it in the metadata 
sample_data(ps)$First_Word <- sub("(\\D+).*", "\\1", sample_data(ps)$Sample_name)
unique(sample_data(ps)$First_Word)

#make sure that the new categorizing worked 
sample_data(ps)

#Check if the merged_physeq object 'ps' has OTU counts and merged correctly 
otu_table(ps)

#Ensuring First_Word is a factor
sample_data(ps)$First_Word <- as.factor(sample_data(ps)$First_Word)
####################################################
#alpha diversity 
plot_richness(ps, measures = c("Shannon", "Simpson"), color="First_Word")

#the same but harder to visualize sample richness and diversity as their groupings
plot_richness(ps, measures = c("Shannon", "Simpson"), color="Sample_name")

#beta diversity to check for sample clustering for species difference 
#Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="First_Word", title="Bray NMDS")

#############top 50 instead of top 20 for abundance plotting 
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
top50
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50
ps.top50 <- prune_taxa(top50, ps.top50)
ps.top50

#visualize abundance 

#clumping them together
plot_bar(ps.top50, x="First_Word", fill="Phylum")
#how the paper had it 
plot_bar(ps.top50, x="Sample_name", fill="Phylum")

##########################binning majors and minors like the paper did for the phylogenetic tree 
#Selecting the top 50 ASVs
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]

#Transform the sample counts to relative abundances and prune to keep only the top 50 ASVs
ps.top50 <- prune_taxa(top50, transform_sample_counts(ps, function(OTU) OTU/sum(OTU)))

#Classifying the ASVs into major and minor groups
#Using MARGIN=2 to apply across columns (ASVs)
major <- top50[apply(otu_table(ps.top50), 2, max) > 0.01]
minor <- top50[apply(otu_table(ps.top50), 2, max) < 0.01]

#checking the size of both groups to see if i agree with paper
length(major)
length(minor)
#got different numbers but may be due to the fact that they are looking at OTUs and me ASVs 
#################################phylogenetic tree
#https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html#:~:text=In%20the%20ape%20package%20(Paradis%20and%20Schliep,of%20data%20(amino%20acid%2C%20nucleotides%2C%20morphological%20data).&text=parsimony%20performs%20tree%20rearrangements%20to%20find%20trees%20with%20a%20lower%20parsimony%20score.
if (!requireNamespace("phyloseq", quietly = TRUE)) install.packages("phyloseq")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")
if (!requireNamespace("msa", quietly = TRUE)) BiocManager::install("msa")
library(phyloseq)
library(ape)
library(phangorn)
library(msa)

#rebinning for fun i guess, cant hurt##################################
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]

top50# Transform the sample counts to relative abundances and prune to keep only the top 50 ASVs
ps.top50 <- prune_taxa(top50, transform_sample_counts(ps, function(OTU) OTU/sum(OTU)))

major <- top50[apply(otu_table(ps.top50), 2, max) > 0.01]
minor <- top50[apply(otu_table(ps.top50), 2, max) <= 0.01]

length(major)
length(minor)
#####################################################################

#only plotting the major ASVs
ps.major <- prune_taxa(major, ps.top50)
ps.major

#using MSA to do MSA on my sequences to find diversion between samples 
seqs <- refseq(ps.major)
alignment <- msa(seqs)

#Building a phylogenetic tree using phangorn 
phang.align <- as.phyDat(alignment, type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Neighbor-Joining tree
fit <- pml(treeNJ, data = phang.align)

#using GTR model to infer the tree 
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic")

#remove tips with NA labels (removing the specific asvs anyway if this fails)
tree_clean <- drop.tip(fitGTR$tree, which(is.na(fitGTR$tree$tip.label)))

#Removing ASVs that were not assigned taxonomy so i dont get ugly NAs ruining the phylogenetic tree 
to_remove <- c("ASV25", "ASV11", "ASV39", "ASV44")  #these asvs were all coming up as NA, ruining my tree
tree_clean <- drop.tip(tree_clean, to_remove)

#keeps my cleaned tree 
ps.major <- prune_taxa(taxa_names(ps.major) %in% tree_clean$tip.label, ps.major)

#merges the clean tree and the normal tree 
ps.major <- merge_phyloseq(ps.major, phy_tree(tree_clean))

#plotting the tree! using phylum as their seperator and genus as their labels (since i coudlnt do species like in the paper)
plot_tree(ps.major, color = "Phylum", label.tips = "Genus", ladderize = TRUE)

#Check the taxonomy table, indeed lacking many species 
tax_table(ps.major)

############################## DeSeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

#this again in case it forgot for some reason
################################################################################
sample_data(ps)$First_Word <- sub("(\\D+).*", "\\1", sample_data(ps)$Sample_name)
unique(sample_data(ps)$First_Word)
sample_data(ps)
length(unique(sample_data(ps)$First_Word))
sample_variables(ps)
otu_table(ps)
sample_data(ps)$First_Word <- as.factor(sample_data(ps)$First_Word)
####################################################################################

#Converting phyloseq object to DESeq2 object using First_word to differentiate them
dds <- phyloseq_to_deseq2(ps, ~ First_Word)  

#Perform DESeq2 analysis on my DeSeq object 
dds <- DESeq(dds)

#Viewing the results of the differential abundance testing
res <- results(dds, alpha = 0.05)  #Setting alpha to 0.05 for adjusted p-value cutoff (i tried 0.1 instead and it didnt change anything)
summary(res)

#Extracting significantly differentially abundant taxa (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05), ]
not_sig_res <- res[which(res$padj > 0.05), ]

sig_res <- sig_res[order(sig_res$log2FoldChange, decreasing = TRUE), ]
not_sig_res <- not_sig_res[order(not_sig_res$log2FoldChange, decreasing = TRUE), ]
sig_res
not_sig_res
#not finding many significant taxa!

#Convert significant results to a dataframe
output <- as.data.frame(sig_res)
View(output) #only 1 loool 

output2 <- as.data.frame(not_sig_res)
View(output2)

#Visualize differentially abundant taxa
# Plotting the log2 fold changes of significantly differentially abundant taxa
# Create an MA plot
plotMA(res, main = "DESeq2 - Differentially Abundant Taxa", ylim = c(-5, 5))


#Accessing and plotting specific taxa
top_taxa <- rownames(sig_res)[1:50]  # Top 50 differentially abundant taxa
print(top_taxa) #only 1 still 

#Check how many significant taxa are available
num_sig_taxa <- nrow(sig_res)
print(paste("Number of significant taxa:", num_sig_taxa))

#Extract the top significant taxa based on availability
top_taxa <- rownames(sig_res)[1:min(50, num_sig_taxa)]

if (length(top_taxa) > 0) {
  plotCounts(dds, gene = top_taxa[1], intgroup = "First_Word")
} else {
  print("No significant taxa to visualize.")
}

dds$Location <- as.factor(dds$First_Word)

plotCounts(dds, gene = top_taxa[1], intgroup = "First_Word")

# Create a volcano plot using ggplot2 
library(ggplot2)

#Convert DESeq2 results to a dataframe for visualization
res_df <- as.data.frame(res)
res_df

#Add a column to indicate significance based on adjusted p-value because its multiple testing 
res_df$significance <- ifelse(res_df$padj <= 0.05, "Significant", "Not Significant")

#visualize using volcano plot 
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
##################################################Abundance and Prevalence filtering
# Abundance filtering
ps_abundance <- filter_taxa(ps, function(x) sum(x) > 10, TRUE)

# Prevalence filtering
ps_filtered <- filter_taxa(ps_abundance, function(x) sum(x > 0) / length(x) > 0.2, TRUE)

#####################################################DESEQ2 but with abundance and prevalence filtering 
#Converting phyloseq object to DESeq2 object using First_word to differentiate them
dds <- phyloseq_to_deseq2(ps_filtered, ~ First_Word)  

#Perform DESeq2 analysis on my DeSeq object 
dds <- DESeq(dds)

#Viewing the results of the differential abundance testing
res <- results(dds, alpha = 0.05)  #Setting alpha to 0.05 for adjusted p-value cutoff (i tried 0.1 instead and it didnt change anything)
summary(res)

#Extracting significantly differentially abundant taxa (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05), ]
not_sig_res <- res[which(res$padj > 0.05), ]

sig_res <- sig_res[order(sig_res$log2FoldChange, decreasing = TRUE), ]
not_sig_res <- not_sig_res[order(not_sig_res$log2FoldChange, decreasing = TRUE), ]


sig_res
not_sig_res
#found 3 differentially abundant taxa this time!!
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   ASV20   23.2021        8.87907   2.75923   3.21795 1.29108e-03 1.63537e-02
# ASV12   49.7551      -10.10731   2.20319  -4.58758 4.48423e-06 8.52004e-05
# ASV9    45.7928      -24.41916   2.37980 -10.26102 1.05586e-24 4.01226e-23

#WOOOOW

#Convert significant results to a dataframe
output <- as.data.frame(sig_res)
View(output) #only 1 loool 

output2 <- as.data.frame(not_sig_res)
View(output2)

#Visualize differentially abundant taxa with filtering 
# Plotting the log2 fold changes of significantly differentially abundant taxa
# Create an MA plot
plotMA(res, main = "DESeq2 - Differentially Abundant Taxa", ylim = c(-5, 5))


#Accessing and plotting specific taxa
top_taxa <- rownames(sig_res)[1:50]  # Top 50 differentially abundant taxa
print(top_taxa) #only 1 still 

#Check how many significant taxa are available
num_sig_taxa <- nrow(sig_res)
print(paste("Number of significant taxa:", num_sig_taxa))

#Extract the top significant taxa based on availability
top_taxa <- rownames(sig_res)[1:min(50, num_sig_taxa)]

if (length(top_taxa) > 0) {
  plotCounts(dds, gene = top_taxa[1], intgroup = "First_Word")
} else {
  print("No significant taxa to visualize.")
}

dds$Location <- as.factor(dds$First_Word)

plotCounts(dds, gene = top_taxa[1], intgroup = "First_Word")

# Create a volcano plot using ggplot2 
library(ggplot2)

#Convert DESeq2 results to a dataframe for visualization
res_df <- as.data.frame(res)
res_df

#Add a column to indicate significance based on adjusted p-value because its multiple testing 
res_df$significance <- ifelse(res_df$padj <= 0.05, "Significant", "Not Significant")

#visualize using volcano plot 
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

#adundance and prevalence filtering helped! 
##################################
#THE END!

