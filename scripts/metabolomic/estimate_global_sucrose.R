# Estimation of global sucrose.R
```{r}
#collect all data from seagrass meadows 
df<-read.csv('Results/pw_sg_species_sugar_concentations.csv')#other sg species
df.2<-read.csv('Results/pw_sg_species_sugar_concentations_KIEL.csv')#observations in kiel
df.3<-read.csv('Results/sugar concentrations/in_edge_out_targeted_sugar_concentations.csv') #in-edge-out
df.3<-df.3[df.3$location=='In',]
df.4<-read.csv('Results/sugar concentrations/diel_pw_sucrose_concentations.csv')
df.5<-read.csv('Results/sugar concentrations/season_targeted_sugar_concentations.csv') #season
df.5<-df.5[!df.5$month=='October',] #remove redundant samples 

df.all<-data.frame(sampleID=c(as.vector(df$Data.File), as.vector(df.2$Data.File),as.vector(df.3$Data.File),as.vector(df.4$Data.File),as.vector(df.5$Data.File)),
                     compound=c(as.vector(df$Compound), as.vector(df.2$Compound),as.vector(df.3$Compound),as.vector(df.4$Compound), as.vector(df.5$Compound)),
                     Concentration_uM=c(df$Concentration_uM, df.2$Concentration_uM,df.3$Concentration_uM,df.4$Concentration_uM, df.5$Concentration_uM))
head(df.all)

sucrose<-df.all[df.all$compound=='Sucrose, octakis(trimethylsily) ether Results','Concentration_uM']

# Only select values with in detection limits of 0.03 uM to 200 uM
range(sucrose)
sucrose[sucrose > 200]<-200
forEst<-sucrose[sucrose<201 & sucrose > 0.02] #conservative estimate because intensities go above calibration scale
range(forEst)

##....................................................
## How much g of sucrose under SG meadows?
#
# g sucrose / L pw * Lpw = g sucrose
##....................................................

##1. Convert sucrose uM concentration to g/L
mol_weight_suc_g_mol<-342.297 #g per mol
mean_sucrose<-mean(forEst) # umol/L
(sucrose_g_L<-mol_weight_suc_g_mol*mean_sucrose*1E-06) #sucrose in grams per L

#2. Calculate L of porewater 

#Values from Fourqurean et al. 2012 NatGeo DOI:10.1038/NGEO1477
seagrass_surface_area_low<-300000 #in km^2
seagrass_surface_area_hi<-600000 ##in km^2
porosity_median<-0.616

sed_depth<-0.0003 #in km, max sed depth of our study
## Calculate km3 of sediment void space
## porosity*Volume sediment = Void space
(Vspace_km3_low<-porosity_median*seagrass_surface_area_low*sed_depth)
(Vspace_km3_hi<-porosity_median*seagrass_surface_area_hi*sed_depth)
Vspace_meadow<-conference_room*porosity_median*sed_depth

## Convert km3 to L of porewater
## 1 m^3 = 1000 L 
## 10^9 m^3 * 1000 = 10^12 

L_hi<-Vspace_km3_hi*1E12
L_low<-Vspace_km3_low*1E12
L_meadow<-Vspace_meadow*1E12

(low_g<-L_low*sucrose_g_L)
(high_g<-L_hi*sucrose_g_L)
(meadow_g<-L_meadow*sucrose_g_L)


# convert to T g of sucrose
low_g*1E-12
high_g*1E-12

