#Kallisto_counts.R
#EM Sogin 
#March 2020
#Normalize kallisto raw counts to gene length 

#Combine raw counts
library(tidyr)
files<-list.files(pattern='.tsv',recursive=T)

transcript_count<-data.frame()
for (i in 1:length(files)){
    print(files[i])
	df<-read.table(files[i],header=T)
	transcript_count<-rbind(transcript_count, data.frame(library=files[i],df))
}

# TPM values aren't TPM per bin, but rather TPM per library, need to adjust to per bin 
transcript_count$library1<-gsub('_output/abundance.tsv','',transcript_count$library)

# rename library names
lookup<-data.frame(library=c("4308_A","4308_B","4308_C","4308_D","4308_E","4308_F","4308_G","4308_H","4308_I"),lib2=c("A_Out","B_Out","C_Out","D_Edge","E_Edge","F_Edge","G_In","H_In","I_In"))
transcript_count$library2<-lookup[match(transcript_count$library1, lookup$library),'lib2']

# get mags
transcript_count$MAG<-gsub("_.*", "",transcript_count$target_id)
tmp<-transcript_count[,!colnames(transcript_count) %in% c('tpm', 'library', 'library1', 'MAG','length', 'eff_length')]
count_table<-tmp %>% spread(key=library2,value=est_counts)
write.csv(count_table,'kallisto_est_counts.csv')

# read in data
kcounts<-read.csv('kallisto_est_counts.csv',header=T)
coords<-read.table(file="gene_coords.gff")
head(coords)
colnames(coords)<-c('target_id','type','CDS','Start','Stop','V6','V7','V8','V9','genelength_kb')

# manipulate data
# gene length per kbp
kcounts$genelength_kb<-coords[match(kcounts$target_id, coords$target_id),'genelength_kb']
head(kcounts)
range(kcounts$genelength_kb)

#Calculate TPM per mag properally, avoids the full library normalization from kallisto
kcounts$mag<-gsub('_.*','',kcounts$target_id)
mags<-unique(kcounts$mag)
TPM_counts<-data.frame()
lib<-c("A_Out","B_Out","C_Out","D_Edge","E_Edge","F_Edge","G_In","H_In","I_In")
for(i in 1:length(mags)){
    tmp<-kcounts[kcounts$mag %in% mags[i],]
    FPK <- tmp[,colnames(tmp) %in% lib] / tmp[,'genelength_kb']
    TPM <- sweep(FPK,2, colSums(FPK)/1000000, "/")
    temp2<-data.frame(geneid=tmp$target_id, mag=tmp$mag , TPM)
    TPM_counts<-rbind(TPM_counts, temp2)
}
       
head(TPM_counts)
save('TPM_counts', file='TPM.Kallisto.March2020.RData')