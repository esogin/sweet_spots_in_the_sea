#Non-targeted data analysis script
## output is taken from the peak picking script. 
## UNTARGETED APPROACH
## VOLCANO PLOT 
df<-pksReduced

rownames(df)<-paste('PCGRP',df$pcgroup, 'mz', round(df$mz,2), 'rt', round(df$rt, 1), sep='_')
ion_lookup<-df[,colnames(df) %in% c('mz','rt',"pcgroup")]
samples<-colnames(df)[grep('X',colnames(df))]
s<-gsub('X','',samples)
df.volcano<-data.frame(TransectLocation=mdata[match(s, mdata$GCMS_ID),'TransectLocation'],Sediment_depth_cm=mdata[match(s, mdata$GCMS_ID),'Sediment_depth_cm'], t(df[,colnames(df) %in% samples]))

head(df.volcano[,1:10])

## Generate anova models 
mod<-lapply(df.volcano[,grep('PCGRP',colnames(df.volcano))], function(x) anova(lm(log(x+0.001) ~ df.volcano$TransectLocation*df.volcano$Sediment_depth_cm)))

pvals<-data.frame(pvalue=unlist(lapply(mod, function(x) x$`Pr(>F)`[1])))
pvals$Ion<-rownames(pvals)
head(pvals)
range(pvals$pvalue)

## Calculate fold changes 
IN<-colMeans(df.volcano[df.volcano$TransectLocation=='In' & df.volcano$Sediment_depth_cm==-10,grep('PC',colnames(df.volcano))])
EDGE<-colMeans(df.volcano[df.volcano$TransectLocation=='Edge' & df.volcano$Sediment_depth_cm==-10,grep('PC',colnames(df.volcano))])
fcdf<-data.frame(IN, EDGE)
fcdf$FC<-fcdf$IN/fcdf$EDGE
fcdf$Ion<-rownames(fcdf)

## Combine
v.merge<-merge(pvals, fcdf, by='Ion')
v.merge$p.adjust<-p.adjust(v.merge$pvalue, method='BH')
range(v.merge$p.adjust, na.rm=T)
range(v.merge$pvalue, na.rm=T)

v.merge$threshold<-as.factor(v.merge$p.adjust < 0.1 & abs(log2(v.merge$FC)) > 2)
v.merge$pcgroup<-gsub('.mz.*','',v.merge$Ion)
v.merge$pcgroup<-gsub('PCGRP.','',v.merge$pcgroup)
v.merge<-v.merge[complete.cases(v.merge),]

## Plot Volcano
ggplot(data=v.merge, aes(x=log2(FC), y =-log10(p.adjust))) + geom_point(data=v.merge[v.merge$threshold=='FALSE',],aes(x=log2(FC), y =-log10(p.adjust)),colour='gray')  + geom_point(data=v.merge[v.merge$threshold=='TRUE',],aes(x=log2(FC), y =-log10(p.adjust),colour=pcgroup)) + geom_vline(xintercept = 0, linetype=2) + geom_hline(yintercept = -log10(0.1), linetype=2, color='coral') + geom_vline(xintercept = -2, linetype=2, color='lightblue') + geom_vline(xintercept = 2, linetype=2, color='lightblue') 

## USE MASS HUNTER TO ID IONS THAT ARE WITHIN THREHOSLD
v.merge.sigions<-v.merge[v.merge$threshold=='TRUE',]
v.merge.sigions$rt<-ion_lookup[match(v.merge.sigions$Ion, rownames(ion_lookup)),'rt']
v.merge.sigions$pcgroup<-ion_lookup[match(v.merge.sigions$Ion, rownames(ion_lookup)),'pcgroup']
metabolites<-ddply(v.merge.sigions, .(pcgroup), summarise, rt=mean(rt))

Imp.mets<-data.frame(v.merge.sigions)
Imp.mets$rt<-metabolites[match(Imp.mets$pcgroup, metabolites$pcgroup),'rt']
#write.csv(Imp.mets, file = 'Important_Ions_IN_EDGE_10CM_Metabolites.csv')
write.csv(Imp.mets, file = 'Important_Ions_IN_EDGE_10CM_Metabolites_p0.1.csv')