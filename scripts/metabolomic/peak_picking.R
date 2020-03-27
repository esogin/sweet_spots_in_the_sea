#Peak Picking Script
#Pick peaks from mzXML data files for untargeted data analysis

##--------------------------------------------------------------------------------------------------
## SET UP WORKING SPACE
##--------------------------------------------------------------------------------------------------

rm(list=ls())
dir<-"path/to/data/analysis/folder"
setwd(file.path(dir, 'Scripts'))

library(xcms)
library(CAMERA)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
source('ExtraFunctions.R')
library(ggrepel)
library(car)
library(agricolae)

## Metabolite database 
setwd(file.path(dir, 'Data'))  

## Metadata
mdata<-read.csv('mdata_salty_samples.csv')
mdata<-mdata[which(mdata$Keep=='keep'),]

##--------------------------------------------------------------------------------
## PEAK PICKING, RETENTION TIME GROUPING & CORRECTION
##--------------------------------------------------------------------------------
setwd(file.path(dir, 'path/to/data'))  

files<-list.files(pattern='.mzXML', recursive = T, full.names=T) 

xs <- xcmsSet(files, method = "matchedFilter", fwhm = 8.4, snthresh = 1,step= 0.25, steps= 2,sigma = 3.56718192627824, max= 500, mzdiff= 1,index= FALSE)
xset1 <- group(xs,method  = "density", bw=2, mzwid= 1, minfrac = 0.5,  minsamp = 1,max = 500) ## Initial peak grouping
xset2<-retcor(xset1, smooth='linear')
plotrt(xset2, leg = FALSE, densplit = TRUE)
xset2 <- group(xset2,method  = "density", bw=2, mzwid= 1, minfrac = 0.5,  minsamp = 1,max = 500)
xset<-fillPeaks(xset2)

## Group peaks in to pseudo-spectra
an<-xsAnnotate(xset)
xsF<-groupFWHM(an, perfwhm=3)

## Select ions
peaks<-getPeaklist(xsF)
peaks<-peaks[order(peaks$pcgroup),]
peaks[is.na(peaks)]<-0


##--------------------------------------------------------------------------------------------------
## DATA MANIPULATION
##--------------------------------------------------------------------------------------------------
 
peaks2<-peaks
## Normalize Peak Data To Ribitol & Remove Samples with Poor Derivitization
Ribitol<-peaks2[peaks2$mz > 217 & peaks2$mz < 218 & peaks2$rt < 1006 & peaks2$rt > 1000,] ## select quantification ion for ribitol 

## Plot Ribitol
ribs<-melt(Ribitol, measure.vars = grep('X', colnames(Ribitol)))
ggplot(ribs, aes(x=variable, y=value)) + geom_bar(stat='identity', width=0.1) + geom_point() + geom_hline(yintercept=1e6)

## Remove samples from dataset with ribitol areas less then 1E6
remove<-as.vector(ribs[which(ribs$value < 1E6), 'variable'])

peaks2<-peaks2[,!colnames(peaks2) %in% remove]
Ribitol<-peaks2[peaks2$mz > 217 & peaks2$mz < 218 & peaks2$rt < 1006 & peaks2$rt > 1000 ,] ## normalize data
r1<-as.vector(t(Ribitol[1,grep('X', colnames(Ribitol))]))
Data<-data.frame(peaks2[grep('X', colnames(peaks2))])
Data.norm<-sweep(Data, MARGIN = 2, STATS = r1, FUN = '/')
df<-data.frame(mz=peaks2$mz, rt=peaks2$rt, pcgroup=peaks2$pcgroup,Data.norm)
dim(df)

head(df)

## Remove Ions Below 146
df<-df[!df$mz < 146,]
dim(df)

## Remove Cholestane and Ribitol PC Groups
head(df)

df<-df[!df$pcgroup %in% c(1,25,5,221),]


## Select one peak per pcgroup 
head(df[,1:10])

df.all.peaks<-df
#df<-df.all.peaks
## Ion select for untargeted analysis

## Rules: 
## 1) mz > 150
#df<-df[df$mz > 150,]
pksLMz<-df[df$mz > 150,]

## 2) Subset by frequency 
pksLMz$count<-rowSums(pksLMz[,grep('X', colnames(pksLMz))] !=0) # count number of times entry appears in dataframe 
pksFreq<-data.frame(pksLMz %>% group_by(pcgroup) %>% top_n(1,count)) #provides only one entry for most metabolites 
pksFreq

## 3) Subset by mean 
pksFreq$means<-rowMeans(pksFreq[,grep('X', colnames(pksFreq))])
pksMean<-data.frame(pksFreq %>% group_by(pcgroup) %>% top_n(1,means))
head(pksMean)
range(pksMean$means)

##Remove very low concentrated peaks 
pksMean<-pksMean[pksMean$means > 0.001,]

## Subset peaks dataframe 
pksReduced<-data.frame(mz=pksMean$mz, rt=pksMean$rt, pcgroup=pksMean$pcgroup,pksMean[,grep('X', colnames(pksMean))])
head(pksReduced)

