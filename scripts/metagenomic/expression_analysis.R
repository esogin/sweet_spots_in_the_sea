# Expression_analysis.R
## EM Sogin
## Update March 2020

# Set up working space
rm(list=ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(packcircles)
library(ggraph)
library(igraph)

#load in data
gene_annots<-read.table('data/Seagrass_MAGs_all_annotations.tsv',sep='\t', header=T,quote="",fill = T)

bin.metadata<-read.csv('data/bin_mdata.csv')
bin.metadata$binid<-gsub('.contigs.fa','',bin.metadata$file_name)
bin.metadata$binid<-gsub('_sub','.sub',bin.metadata$binid)
head(bin.metadata)

#load normalized count data
load('data/TPM.Kallisto.Jan2020.RData') # metabat.241 didn't have any mapping reads in lib C

##---------------------------------------------------------------------------------
## Manipulate data
##---------------------------------------------------------------------------------

#for kallisto, add in annotations
TPM_counts$bin.location<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'location']
TPM_counts$taxonomy<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'GTDB_taxonomy']
TPM_counts$EggNOG<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'EggNOG']
TPM_counts$EggNOG_function<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'EggNOG_function']
TPM_counts$Uniprot<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'Uniprot']
TPM_counts$CAZY_HMM<-gene_annots[match(TPM_counts$geneid, gene_annots$Geneid),'CAZY_HMM']
head(TPM_counts)

TPM_counts$genome.size<-bin.metadata[match(TPM_counts$mag, bin.metadata$binid),'Genome.size..bp.']

## save as a new dataframe
TPM_counts.df<-data.frame(TPM_counts)

# Convert to long format
df.long<-TPM_counts.df %>% 
		gather(key='library',value='TPM', -geneid, -mag, -bin.location, -taxonomy, -EggNOG, -EggNOG_function, -Uniprot, -CAZY_HMM, -genome.size)

# add in library location 
df.long$lib.location<-gsub('[ABCDEFGHI]_','',df.long$library)
head(df.long);tail(df.long)

##---------------------------------------------------------------------------------
## Target genes
##---------------------------------------------------------------------------------

#Transport genes
#ABC Transport - expensive energy transport system
abc<-unique(gene_annots[c(grep('sucrose ABC transporter',gene_annots$Uniprot),
                          grep('COG0395',gene_annots$EggNOG),
                          grep('sugar ABC transporter',gene_annots$Uniprot)),'Geneid'])

#cscB - sucrose permease 
fucP<-unique(gene_annots[grep('COG0738',gene_annots$EggNOG),'Geneid']) #fucose transporter
cscB<-unique(gene_annots[grep('ose permease',gene_annots$RAST),'Geneid'])
hyp<-gene_annots[grep('Hypothetical sugar permease',gene_annots$RAST),'Geneid']
permease<-unique(c(fucP,cscB,hyp))

#scrY - sucrose porin 
porin<-unique(gene_annots[c(grep('LamB',gene_annots$Pfam),
                            grep('COG3659',gene_annots$EggNOG)),'Geneid'])

#scrA PTS sucrose specific
PTS<-unique(gene_annots[c(gene_annots$EggNOG %in% c('COG1263','COG1264','ENOG4105C5Y'),
                   grep('PTS sugar',gene_annots$Uniprot),
                   grep('EIIB',gene_annots$Pfam),
                   grep('EIIC',gene_annots$Pfam)),'Geneid'])

#Sugar degredation genes
#scrP - also known as sucP, Sucrose phosphorylase
#ams -  amylosucerase
phosphorylase<-gene_annots[c(grep('Sucrose phosphorylase',gene_annots$RAST),
                    grep('crose phosphorylase',gene_annots$Uniprot),
                    grep('GH13_4',gene_annots$CAZY_HMM),
                    grep('mylosucrase',gene_annots$Uniprot)),]
phosphorylase<-unique(phosphorylase$Geneid)

#scrB - Sucrose phosphotransfer (coupled to scrA reaction), sucrose 6-phosphate hydrolase
scrB<-gene_annots[c(grep('Sucrose-6-phosphate hydrolase ',
                    gene_annots$RAST),grep('GH32',gene_annots$CAZY_HMM), 
                    grep('ucrose-6-phosphate',gene_annots$Uniprot),
                    grep('Glycosyl hydrolases family 32',gene_annots$Pfam)),]
scrB<-scrB[!scrB$EggNOG_function=='Sulfatase',]
scrB<-unique(scrB$Geneid)                  

hydrolase<-gene_annots[grep('COG1621',gene_annots$EggNOG),'Geneid']

recA<-unique(gene_annots[c(grep('recA',gene_annots$Pfam),grep('RecA',gene_annots$RAST), grep('COG0468',gene_annots$EggNOG)),'Geneid'])

sulfatase<-gene_annots[grep('COG3119',gene_annots$EggNOG),'Geneid']

targets<-data.frame(geneid=c(recA,abc,permease,porin,PTS,phosphorylase,hydrolase,sulfatase), 
                    gene_name=rep(c('recA','abc','permease','porin','PTS','phosphorylase','hydrolase','sulfatase'),
                                  c(length(recA),length(abc),length(permease),length(porin),
                                    length(PTS),length(phosphorylase),length(hydrolase),length(sulfatase))),
                    cat=rep(c('recA','transport','transport','transport','transport','degredation','degredation','sulfatase'),
                            c(length(recA),length(abc),length(permease),length(porin),
                              length(PTS),length(phosphorylase),length(hydrolase),length(sulfatase))))
head(targets)


## Visualization
df.tgs<-df.long[df.long$geneid %in% targets$geneid,]
df.tgs$gene_name<-targets[match(df.tgs$geneid, targets$geneid),'gene_name']
df.tgs$cat<-targets[match(df.tgs$geneid, targets$geneid),'cat']


summary.tgs<-df.tgs %>% 
			group_by(mag,lib.location,gene_name,cat) %>% 
			summarise(TPMs_per_gene=sum(TPM))

RecA<-summary.tgs[summary.tgs$gene_name=='recA',]

summary.tgs$lib.location<-factor(summary.tgs$lib.location, ordered=T, levels=c('In','Edge','Out'))

summary.tgs<-data.frame(summary.tgs)

summary.tgs$cat<-factor(summary.tgs$cat, ordered = T, levels=c('transport','degredation','recA','sulfatase'))
summary.tgs.plot<-summary.tgs[summary.tgs$TPMs_per_gene>0,]
ggplot(summary.tgs.plot, aes(x=gene_name, y=log10(TPMs_per_gene), fill=lib.location)) + 
  geom_hline(yintercept=log10(2181.4))+
  geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.75))+
  geom_point(size=2,alpha=0.15,position=position_jitterdodge(jitter.width =0.2,dodge.width =  0.75)) + 
  facet_grid(~cat, scales='free',space='free') + scale_fill_manual(values=l2p) + theme_new + ylab('log10(Total TPM)') + xlab('')


##---------------------------------------------------------------------------------
## CAZYme Analysis 
##---------------------------------------------------------------------------------
czy.sub<-df.long[grep('[[:alnum:]]',df.long$CAZY_HMM),]
czy.sub$czy_besthit<-gsub(' \\[.*','',czy.sub$CAZY_HMM) # gets best hit
cazy<-data.frame(geneid=czy.sub$geneid, mag=czy.sub$mag,tax=czy.sub$taxonomy, bin.location=czy.sub$bin.location,
                 cazy_hit=czy.sub$czy_besthit,genome_size=czy.sub$genome.size,TPM=czy.sub$TPM, library=czy.sub$library, lib.location=czy.sub$lib.location)

#cazy general category
cazy$general_hit<-gsub("[0-9].*",'',cazy$cazy_hit)

# GHs - carbon degrading genes
cazy.gh<-cazy[cazy$general_hit=='GH',]
head(cazy.gh)

# add in substrate info (simple vs. complex)
substrates<-read.csv('data/cazymes_annotations.csv',header=T)
cazy.gh$substrate2<-substrates[match(cazy.gh$cazy_hit,substrates$Family),'Cat3']
cazy.gh$substrate2<-as.vector(cazy.gh$substrate2)
cazy.gh[which(cazy.gh$substrate2=='plant compound'),'substrate2']<-'complex' 
cazy.gh[which(cazy.gh$substrate=='sucrose'),'substrate2']<-'sucrose'

# Get total TPM counts per genome per substrate category per library 
ss<-cazy.gh %>% 
	group_by(mag,substrate2, tax, genome_size,lib.location) %>% 
	summarise(total_TPM=sum(TPM))
ss$lib.location<-factor(ss$lib.location, ordered=T, levels=c('In','Edge','Out'))
ss<-ss[complete.cases(ss$TPM_per_RecA),] #Remove NA and INF
ss<-ss[!ss$TPM_per_RecA=='Inf',]

#Calculate ratio of TPM Sucrse vs. TPM other
TPM.sucrose<-data.frame(ss[ss$substrate2=='sucrose',])
TPM.sucrose<-data.frame(TPM.sucrose %>% group_by(mag, tax, lib.location) %>% summarize(TPMsuc=sum(total_TPM)))

TPM.other<-ss[!ss$substrate2=='sucrose',]
TPM.other<-data.frame(TPM.other %>% group_by(mag, tax, lib.location) %>% summarize(TPMother=sum(total_TPM)))

TPM.sucrose$other<-TPM.other[match(paste(TPM.sucrose$mag, TPM.sucrose$lib.location),paste(TPM.other$mag, TPM.other$lib.location)),'TPMother']

#calculate the ratio between sucrose TPM / bin to other GHs per bin
TPM.sucrose$suc_other_ratio<-TPM.sucrose$TPMsuc/TPM.sucrose$other
plot.ty<-TPM.sucrose[!TPM.sucrose$suc_other_ratio=='Inf',]
plot.ty<-plot.ty[!is.na(plot.ty$suc_other_ratio),]

subpoints<-plot.ty[plot.ty$suc_other_ratio > 1,]
anit_subpoints<-plot.ty[plot.ty$suc_other_ratio < 1,]

bx<-ggplot() + 
  geom_hline(yintercept = log10(1),linetype=2,color='gray') + 
  geom_boxplot(data=plot.ty,outlier.shape = NA,fill='gray90',aes(x=lib.location, y=log10(suc_other_ratio))) +
  geom_jitter(data=anit_subpoints,width=0.1,height=0,alpha=0.5,size=3,color='black',aes(x=lib.location, y=log10(suc_other_ratio))) + 
  geom_jitter(shape=21,data=subpoints,width=0.1,height=0,size=3,aes(x=lib.location, y=log10(suc_other_ratio),fill=tax))+
  scale_fill_manual(values=c('brown1','deepskyblue','gold')) +
  guides(color=F,size=F,fill=F) + 
  theme_new + ylab('log10(TPM sucrose per TPM other)') + xlab('')
bx


## Circle plotting 
## Visualize Accumulative TPM 
ss.forCicles<-ss[ss$total_TPM>0,]

ss.forCicles$path1<-paste('Sed',ss.forCicles$lib.location,sep='.')
ss.forCicles$path2<-paste('Sed',ss.forCicles$lib.location, ss.forCicles$mag,sep='.')
ss.forCicles$path3<-paste('Sed',ss.forCicles$lib.location, ss.forCicles$mag,ss.forCicles$substrate2,sep='.')

#Edges
e0<-data.frame(node1='Sed',node2=unique(ss.forCicles$path1))
e1<-data.frame(node1=ss.forCicles$path1, node2=ss.forCicles$path2)
e1<-e1[!duplicated(e1),]
e2<-data.frame(node1=ss.forCicles$path2,node2=ss.forCicles$path3)

edges<-rbind(e0,e1,e2)
edges$node1<-as.vector(edges$node1)
edges$node2<-as.vector(edges$node2)

#Nodes
nodes<-data.frame(name=unique(c(edges$node1, edges$node2)))
nodes$size<-ss.forCicles[match(nodes$name,ss.forCicles$path3),'total_TPM']
nodes[is.na(nodes$size),'size']<-1

nodes$shortName<-'bin'
nodes[grep('simple',nodes$name),'shortName']<-'simple'
nodes[grep('complex',nodes$name),'shortName']<-'complex'
nodes[grep('sucrose',nodes$name),'shortName']<-'sucrose'

nodes$labels<-''
nodes[nodes$name=='Sed.Edge','labels']<-'Edge'
nodes[nodes$name=='Sed.In','labels']<-'In'
nodes[nodes$name=='Sed.Out','labels']<-'Out'
nodes[nodes$name %in% c('Sed.In.metabat.241','Sed.Edge.metabat.241','Sed.Out.metabat.241'),'labels']<-'241'
nodes[nodes$name %in% c('Sed.In.metabat.142','Sed.Edge.metabat.142','Sed.Out.metabat.142'),'labels']<-'142'
nodes[nodes$name %in% c('Sed.In.metabat.62','Sed.Edge.metabat.62','Sed.Out.metabat.62'),'labels']<-'62'

# Make Graph 
g <- graph_from_data_frame(
  edges,
  directed = TRUE,vertices = nodes
)

set.seed(10)
bubbles<-ggraph(g, layout='circlepack',weight=size) + 
  geom_node_circle(aes(fill=shortName),alpha=0.75)+
  scale_fill_manual(values=c('white','azure4','azure2','darkblue')) + 
  geom_node_text(aes(label=labels), repel = T) + 
  theme_void() + ggtitle('Total TPM') + guides(alpha=F)
bubbles

#END