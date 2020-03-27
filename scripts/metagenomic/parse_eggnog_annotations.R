# Rscript to parse emapper annotations
# Description

# Rscript will import emapper results in current directory and combine into a single dataframe
# it will then loop through each bin and count the COG categories per bin and provide a long dataframe
# with columns: cog category, frequence of cog call, bin id
# In addition it will provide dataframes that provide total number of genes in emapper files
# as well as total number of cog categories, which can be used for normalization 
# all data will be saved as an Rdata file for import into an r environment for further analysis. 
# Complete datatable with merged emapper results will be saved as "combined_gene_annotations.csv" in working directory

## Loop through files to get emapper results
## Set up working space 
rm(list=ls())


files<-list.files(pattern='emapper.annotations_simplified')
bin_file<-gsub(pattern = '_orfs.*','.fa',files)
bin_names<-gsub('.fa','',bin_file)

#colnumn names - need to addjust depending on how the data were cut from original DF
# Select emapper columns of interest 
# 1.  query_name, 2. seed eggNOG ortholog, 3. seed ortholog evalue, 7. Gene Ontology terms
# 8. EC number, 9. KEGG_Ko, 10. KEGG_Pathway, 11. KEGG_Module, 16. CAZy, 21. COG Functional Category
# see https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2 for other columns that you could use
names<-c('query_name','seed eggNOG ortholog', 'seed ortholog evalue','Gene Ontology terms','EC number','KEGG_ko',' KEGG_Pathway','KEGG_Module','CAZy','COG Functional Category')
keep<-c('query_name','KEGG_ko','COG Functional Category','CAZy', 'EC number') # select columns to keep in analysis for COGs,Kos CAZy and EC analysis

## loop through each datafile in directory and import data into R as a new dataframe
combined_gene_annotations<-data.frame()
for (i in 1:length(files)){
  #Import dataframe
  df<-read.table(files[i], sep='\t',fill = T)
  head(df)
  colnames(df)<-names
  #Manipulate dataframe 
  tmp<-df
  tmp$bin_emapper_name<-bin_file[i]
  combined_gene_annotations<-rbind(combined_gene_annotations,data.frame(tmp))
}
write.csv(combined_gene_annotations, 'combined_gene_annotations.csv')

## Get out COG categories per bin dataframe
data<-combined_gene_annotations #rename for ease of use

#save a table of total number of genes per bin in eggnog. you can also import this information if you have it saved else where
num_orfs<-data.frame(table(data$bin_emapper_name))

data$COG.Functional.Category<-as.vector(data$COG.Functional.Category)
COG.df<-data[!is.na(data$COG.Functional.Category),] #remove NA values
COG.df<-data[!data$COG.Functional.Category=="",] #remove empty values 

#get frequency of cog counts per bin table 
#get frequency of cog counts per bin table 
cog<-data.frame()
for (i in 1:length(bin_file)){
  tmp<-COG.df[COG.df$bin_emapper_name==bin_file[i],]
  cog.tmp<-data.frame(table(unlist(strsplit(x = tmp$COG.Functional.Category, split = character(0))))) # splits character vector by character and gets frequency count per category
  cog<-rbind(cog, data.frame(cog.tmp, bin=bin_file[i]))
}

total_cogs_per_bin<-aggregate(cog$Freq, by=list(cog$bin), FUN=sum)
save.image(file='emapper_annotations.RData')
#End