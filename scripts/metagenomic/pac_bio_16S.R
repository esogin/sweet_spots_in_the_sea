#Full length 16S Analysis using DADA2
#EM Sogin
#Update March 2020
#Description: Analysis of pacbio data with dada2 package

#Set up working space
rm(list=ls())

#libraries
library(ggplot2)
library(ape)
library(dada2)
library(phyloseq)
library(ShortRead)
library(Biostrings)

#Others
path<-"/home/maggie/Documents/Projects/MS1_Seagrass/FullLenth_16S/Analysis/"
path_to_data<-file.path(path, 'Data','CCS_10Passes','fastq')
path.rds <- "Results/RDS/"
path.out<-'Results/'
files<-list.files(path_to_data, pattern="fastq", full.names = T)
GM3F<-'AGAGTTTGATCMTGGC'
GM4R<-"TACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())

##-----------------
##Process data with dada2 package

# 1. Remove primers
nops<-file.path(path_to_data, 'noprimers',basename(files))
for(i in seq_along(files)) {
  fn <- files[[i]]; nop <- nops[[i]]
  dada2:::removePrimers(fn, nop, primer.fwd=GM3F, primer.rev=dada2:::rc(GM4R), orient=TRUE, verbose=T)
}

# 2. Inspect sequence length distribution
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
summary(lens)

#3. Filter data to control for expected sequence lenghts (at least 1000 nt) and quality control.
#https://academic.oup.com/bioinformatics/article/31/21/3476/194979 for more info on EE filtering of data
filts <- file.path(path_to_data, "noprimers", "filtered", basename(files))
track <- filterAndTrim(nops, filts, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE,maxEE=2)
track

# Plot the quality of the basepairs after filtering for each sample
#looks like the quality kicks out after 1600 bp, good place to top the sequence analysis
plotQualityProfile(filts)
plotQualityProfile(filts[1])

#4. Learn error rates
err <- learnErrors(filts, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
plotErrors(err)

saveRDS(err, file.path(path.rds, "errors.rds"))
readRDS(file.path(path.rds, "errors.rds"))

#5. Dereplicate fastq files and run dada2
drp <- derepFastq(filts)
dd <- dada(drp, err=err, multithread=TRUE, BAND_SIZE=32, pool=T)
saveRDS(dd, file = "Results/dada_ASV_full_data.rds")
dd<-readRDS( "Results/dada_ASV_full_data.rds")
dd.rare<-readRDS('Results/dada_ASV.rds')

#6. Make sequence table 
st <- makeSequenceTable(dd); dim(st)
str(st)
rowSums(st)

#7. Assign taxonomy
tax <- assignTaxonomy(st, "~/tax/GTDB_bac-arc_ssu_r86.fa.gz", multithread=TRUE) 

#8. Check for chimeras
bim2 <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim2)
sum(st[,bim2])/ sum(st)

saveRDS(st, 'Results/sequence_table.rds')
saveRDS(tax,'Results/sediment_tax_gtbtk.rds')

#9. Get count table
taxa<-as.data.frame(taxa)
taxa$ASV<-paste("ASV",seq(1:nrow(taxa)), sep="_")

asv_tab <- t(st)
names<-taxa[match(rownames(asv_tab), rownames(taxa)),'ASV']
row.names(asv_tab) <-names
head(asv_tab)
write.table(asv_tab, "Results/ASVs_counts.txt", sep="\t", quote=F, col.names=NA)

#10. Other tasks for getting data into working format
# incorperate count data and sample names into fasta headers
asv_df<-data.frame(asvs=rownames(asv_tab), asv_tab)
asv_long<-reshape2::melt(asv_df)
asvs<-unique(asv_long$asvs)

headers.fa<-data.frame()
for (i in 1:length(asvs)){
  v1<-asvs[i]
  subset<-asv_long[asv_long$asvs==asvs[i] & asv_long$value > 0,]
  s<-paste(subset$variable,'_size=',subset$value,sep="")
  v2<-paste(t(matrix(s)), collapse = ";")
  headers.fa<-rbind(headers.fa, data.frame(v1, v2))
}
headers.fa$headers2<-paste('>',headers.fa$v1,';',headers.fa$v2, sep='')

#match sequences and asv ideas with new fasta headers
taxa$seqs<-rownames(taxa)
fasta.df<-headers.fa
fasta.df$seqs<-taxa[match(fasta.df$v1, taxa$ASV),'seqs']

#Make fasta file
asv_fasta<-c(rbind(fasta.df$headers2, fasta.df$seqs))
write(asv_fasta, "Results/ASVs.fa")


#11. Make phyloseq object
#OTU table 
head(asv_tab)
otus<-otu_table(asv_tab, taxa_are_rows = T)

#sample table
samples<-data.frame(samples=colnames(asv_tab), location=rep(c('Out', 'Edge','In'), c(3,3,3)))
samps<-sample_data(samples)
rownames(samps)<-samples$samples

#Taxa Table
taxa.df<-data.frame(taxa[,colnames(taxa) %in% c("Kingdom" ,"Phylum",  "Class",   "Order",   "Family",  "Genus" )])
rownames(taxa.df)<-taxa$ASV
tax.mat<-as.matrix(taxa.df)
tx<-tax_table(tax.mat)

#merge
ps<-merge_phyloseq(samps, otus, tx)
save(list=c('ps','taxa'), file = 'Results/16s_phyloseq.RData')

#END