###PERMUTATION
setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

## MNI ##           ## MAX analyses begins line 220##
options(stringsAsFactors = FALSE)
modRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Mod.csv", quote="\"")
modRfemur[31,1] = "BM 550" #renames individual 550 as BM 550
modRfemur[30,1] = "BM 570" #renames individual 570 as BM 570
subRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Sub.csv", quote="\"")
subRfemur<-subRfemur[-14,] #remove APS-2
subRfemur<-subRfemur[-13,] #remove TAO-66-22
subRfemur<-subRfemur[-3,] #remove TAO-66-28, only dated R femurs left

subRfemur_temp<-subRfemur[,-10]
subRfemur_temp<-subRfemur_temp[,-c(2:8)] #removed all measurements except FBE and DML2
subRfemur_temp<-subRfemur_temp[-c(7:12),]
subRfemur_temp<-subRfemur_temp[-4,]
subRfemur_temp<-subRfemur_temp[-1,] #removed all individuals with NAs for all measurements

modRfemur_temp<-modRfemur[,-10]
modRfemur_temp<-modRfemur_temp[,-c(2:8)] #removed all measurements except FBE and DML2
modRfemur_temp<-modRfemur_temp[-28,]
modRfemur_temp<-modRfemur_temp[-26,]
modRfemur_temp<-modRfemur_temp[-25,]
modRfemur_temp<-modRfemur_temp[-23,]
modRfemur_temp<-modRfemur_temp[-12,]
modRfemur_temp<-modRfemur_temp[-10,]
modRfemur_temp<-modRfemur_temp[-6,]
modRfemur_temp<-modRfemur_temp[-2,] #removed all individuals with NAs for all measurements, 23 remaining

#subsample mods to same number of individuals as subs
modRfemur_indiv_sample<-sample(nrow(modRfemur_temp), 4, replace = FALSE)
modRfemur_indiv<-modRfemur_temp[modRfemur_indiv_sample,]
#have to make modern subset look like subfossil
modRfemur_indiv[1,2] = NA
modRfemur_indiv[4,2] = NA

library(EnvStats)
geoMeanFBE<-geoMean(modRfemur_indiv[,2], na.rm = TRUE)
geoMeanDML2<-geoMean(modRfemur_indiv[,3], na.rm = TRUE)
geoMeans<-c("geoMeans", geoMeanFBE, geoMeanDML2)
modRfemur_indiv<-rbind(modRfemur_indiv, geoMeans)


FBE_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,2]))/(as.numeric(modRfemur_indiv[5,2]))-1))
  FBE_foldchange_mod[a]<-c(z)
}

FBE_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,2]))/(as.numeric(modRfemur_indiv[5,2]))-1))
  FBE_foldchange_sub[a]<-c(z)
} 

DML2_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,3]))/(as.numeric(modRfemur_indiv[5,3]))-1))
  DML2_foldchange_mod[a]<-c(z)
}

DML2_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,3]))/(as.numeric(modRfemur_indiv[5,3]))-1))
  DML2_foldchange_sub[a]<-c(z)
} 

Rfem_foldchange_mod<-cbind(modRfemur_indiv[,1], FBE_foldchange_mod, DML2_foldchange_mod)
Rfem_foldchange_mod<-Rfem_foldchange_mod[-5,] #removes calculated Medians/geoMeans row
Rfem_foldchange_sub<-cbind(subRfemur_temp[,1], FBE_foldchange_sub, DML2_foldchange_sub)
Rfem_foldchange<-rbind(Rfem_foldchange_mod, Rfem_foldchange_sub)
Rfem_pops<-c(rep("Modern", 4), rep("Subfossil", 4))
Rfem_foldchange<-cbind(Rfem_pops, Rfem_foldchange)
Rfem_foldchange<-data.frame(Rfem_foldchange)
names(Rfem_foldchange) = c("Population", "Individual", "FBE", "DML2")

Rfemur<-Rfem_foldchange[,-2] #remove Individual
Rfemur<-Rfemur[,-1] #remove Population

Rfemur_averages<-c()
for (a in 1:nrow(Rfemur)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Rfemur[a,]), na.rm = TRUE) #arithmetic mean
  #library(EnvStats)
  #z<-geoMean(as.numeric(Rfemur[a,]), na.rm = TRUE) #geometric mean
  Rfemur_averages[a]<-c(z)
}
Rfemur<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
Rfemur<-data.frame(Rfemur)

for(n in 1:9999){
  #subsample mods to same number of individuals as subs
  modRfemur_indiv_sample<-sample(nrow(modRfemur_temp), 4, replace = FALSE)
  modRfemur_indiv<-modRfemur_temp[modRfemur_indiv_sample,]
  #have to make modern subset look like subfossil
  modRfemur_indiv[1,2] = NA
  modRfemur_indiv[4,2] = NA
  library(EnvStats)
  geoMeanFBE<-geoMean(modRfemur_indiv[,2], na.rm = TRUE)
  geoMeanDML2<-geoMean(modRfemur_indiv[,3], na.rm = TRUE)
  geoMeans<-c("geoMeans", geoMeanFBE, geoMeanDML2)
  modRfemur_indiv<-rbind(modRfemur_indiv, geoMeans)
  
  FBE_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,2]))/(as.numeric(modRfemur_indiv[5,2]))-1))
    FBE_foldchange_mod[a]<-c(z)
  }
  
  FBE_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,2]))/(as.numeric(modRfemur_indiv[5,2]))-1))
    FBE_foldchange_sub[a]<-c(z)
  } 
  
  DML2_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,3]))/(as.numeric(modRfemur_indiv[5,3]))-1))
    DML2_foldchange_mod[a]<-c(z)
  }
  
  DML2_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,3]))/(as.numeric(modRfemur_indiv[5,3]))-1))
    DML2_foldchange_sub[a]<-c(z)
  } 
  
  Rfem_foldchange_mod<-cbind(modRfemur_indiv[,1], FBE_foldchange_mod, DML2_foldchange_mod)
  Rfem_foldchange_mod<-Rfem_foldchange_mod[-5,]
  Rfem_foldchange_sub<-cbind(subRfemur_temp[,1], FBE_foldchange_sub, DML2_foldchange_sub)
  Rfem_foldchange<-rbind(Rfem_foldchange_mod, Rfem_foldchange_sub)
  Rfem_pops<-c(rep("Modern", 4), rep("Subfossil", 4))
  Rfem_foldchange<-cbind(Rfem_pops, Rfem_foldchange)
  Rfem_foldchange<-data.frame(Rfem_foldchange)
  names(Rfem_foldchange) = c("Population", "Individual", "FBE", "DML2")
  
  Rfemur_permute<-Rfem_foldchange[,-2] #remove Individual
  Rfemur_permute<-Rfemur_permute[,-1] #remove Population
  
  Rfemur_averages<-c()
  for (a in 1:nrow(Rfemur_permute)){
    if ((a %% 1000) == 0){print (a)}
    z<-mean(as.numeric(Rfemur_permute[a,]), na.rm = TRUE) #arithmetic mean
    Rfemur_averages[a]<-c(z)
  }
  Rfemur_permute<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
  Rfemur<-rbind(Rfemur, Rfemur_permute)
  if ((n %% 1000) == 0){print (n)}
}

names(Rfemur) = c("Population", "Individual", "RFemAverages")
##################################################
options(stringsAsFactors = FALSE)
Rfemur<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/RFemurMNIPermute.txt", quote="\"")

sub_permute<-subset(Rfemur, Population == "Subfossil")
mod_permute<-subset(Rfemur, Population == "Modern")

average_RFemur_subfossil<-mean(as.numeric(sub_permute[,3]), na.rm = TRUE) #0.1032037
#average_MNI_subfossil from geoMeans script = 0.1029081

quantile(as.numeric(mod_permute[,3]))
#           0%           25%           50%           75%          100% 
#-1.956277e-01 -3.107681e-02 -6.252257e-05  3.073210e-02  2.275600e-01
quantile(as.numeric(mod_permute[,3]), c(0.05, 0.95))
#           5%         95% 
#  -0.08542002  0.09930832

#hist(Rfem_pVal_permute, main = "t-test p-values - Mod vs. Sub R Femur Agg Scores", ylim = c(0, 5000))
hist(as.numeric(mod_permute[,3]), col=rgb(1,0,0,0.5), xlim=c(-0.3, 0.3), ylim=c(0, 12500), main = "R Fem Mod 3D Permute", xlab="Aggregate Scores")
#hist(as.numeric(sub_permute[,3]), col=rgb(0,0,1,0.5), add=T)
#points(0.07179621, 0, pch=19, col="#cec3ae") #subfossil observed value
points(0.1032037, 0, pch=19, col="#cec3ae") #MNI subfossil observed value from permute script
points(0.1029081, 0, pch=19, col="#cec3ae") #MNI subfossil observed value from geoMeans script
#points(0.0198988281883758, 0, pch=19, col="#88ccee") #TAO-66-26 AS
#points(0.020560842806777, 0, pch=19, col="#cc6677") #TAO-66-29 AS
#points(0.126896252573764, 0, pch=19, col="#ddcc77") #TAO-66-32 AS
#points(0.244276570389818, 0, pch=19, col="#3b7834") #TAO-66-33 AS
#points(0.00169789, 0, pch=19, col="red") #modern observed value
points(0.09930832, 0, pch=19, col="black") #95th percentile


#average_MNI_subfossil = 0.1029081
#average_modern = 0.008782976
#average_MNI_modern = 0.00169789
#TAO-66-26 AS = 0.0198988281883758
#TAO-66-29 AS = 0.020560842806777
#TAO-66-32 AS = 0.126896252573764
#TAO-66-33 AS = 0.244276570389818

prob_great_R<-which(as.numeric(mod_permute[,3] ) >= 0.1032037) #1820/40000 = 0.0455
prob_less_R<-which(as.numeric(mod_permute[,3] ) < 0.1032037) #38180/40000 = 0.9545
prob_NA_R<-which(mod_permute[,3] == "NaN") #0/40000 = 0.06703333

prob_great_R<-which(as.numeric(mod_permute[,3] ) >= 0.1029081) #geoMeans: 1842/40000 = 0.04605
prob_less_R<-which(as.numeric(mod_permute[,3] ) < 0.1029081) #geoMeans: 38158/40000 = 0.95395
prob_NA_R<-which(mod_permute[,3] == "NaN") #geoMeans: 0/40000 = 0.06703333

a<-mean(as.numeric(mod_permute[,3]))
s<-sd(as.numeric(mod_permute[,3]))
n<-nrow(mod_permute)
error<-qnorm(0.975)*s/sqrt(n)
left<-a-error #0.000998
right<-a+error #0.002057

write.table(Rfemur, "RFemurMNIPermute.txt", sep = "\t", row.names = TRUE, col.names = TRUE)





setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
## MAX ##
modRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Mod.csv", quote="\"")
modRfemur[31,1] = "BM 550" #renames individual 550 as BM 550
modRfemur[30,1] = "BM 570" #renames individual 570 as BM 570
subRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Sub.csv", quote="\"")
subRfemur<-subRfemur[-14,] #remove APS-2
subRfemur<-subRfemur[-13,] #remove TAO-66-22
subRfemur<-subRfemur[-3,] #remove TAO-66-28, only dated R femurs left

subRfemur_temp<-subRfemur[,-10]
subRfemur_temp<-subRfemur_temp[,-c(2:5)] #removed all columns without any subfossil measurements
subRfemur_temp<-subRfemur_temp[-11,]
subRfemur_temp<-subRfemur_temp[-c(7:9),]
subRfemur_temp<-subRfemur_temp[-4,] #removed all individuals with NAs for all measurements

modRfemur_temp<-modRfemur[,-10]
modRfemur_temp<-modRfemur_temp[,-c(2:5)] #removed all columns without any subfossil measurements
modRfemur_temp<-modRfemur_temp[-26,]
modRfemur_temp<-modRfemur_temp[-25,]
modRfemur_temp<-modRfemur_temp[-23,]
modRfemur_temp<-modRfemur_temp[-12,]
modRfemur_temp<-modRfemur_temp[-10,]
modRfemur_temp<-modRfemur_temp[-2,] #removed all individuals with NAs for all measurements


modLfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLFemur_Mod.csv", quote="\"")
modLfemur[31,1] = "BM 550" #renames individual 550 as BM 550
modLfemur[30,1] = "BM 570" #renames individual 570 as BM 570
subLfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLFemur_Sub.csv", quote="\"")
subLfemur<-subLfemur[-14,] #remove APS-2
subLfemur<-subLfemur[-13,] #remove TAO-66-22
subLfemur<-subLfemur[-11,] #remove TAO-66-21
subLfemur<-subLfemur[-7,] #remove TAO-66-33
subLfemur<-subLfemur[-6,] #remove TAO-66-32
subLfemur<-subLfemur[-4,] #remove TAO-66-29
subLfemur<-subLfemur[-3,] #remove TAO-66-28
subLfemur<-subLfemur[-2,] #remove TAO-66-26
subLfemur<-subLfemur[-1,] #remove TAO-66-25, only dated L femurs left

subLfemur_temp<-subLfemur[,-c(2:4)] #removed all columns without any subfossil measurements

modLfemur_temp<-modLfemur[,-c(2:4)] #removed all columns without any subfossil measurements
modLfemur_temp<-modLfemur_temp[-19,]
modLfemur_temp<-modLfemur_temp[-c(12:14),]
modLfemur_temp<-modLfemur_temp[-10,]
modLfemur_temp<-modLfemur_temp[-6,]
modLfemur_temp<-modLfemur_temp[-4,]
modLfemur_temp<-modLfemur_temp[-2,] #removed all individuals with NAs for all measurements


modRhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRHumerus_Mod.csv", quote="\"")
modRhum[31,1] = "BM 550" #renames individual 550 as BM 550
modRhum[30,1] = "BM 570" #renames individual 570 as BM 570
subRhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRHumerus_Sub.csv", quote="\"")
subRhum<-subRhum[-2,] #remove TAO-66-1

subRhum_temp<-subRhum[,-7]
subRhum_temp<-subRhum_temp[,-c(2:4)] #removed all columns without any subfossil measurements

modRhum_temp<-modRhum[,-7]
modRhum_temp<-modRhum_temp[,-c(2:4)] #removed all columns without any subfossil measurements
modRhum_temp<-modRhum_temp[-26,]
modRhum_temp<-modRhum_temp[-25,]
modRhum_temp<-modRhum_temp[-23,]
modRhum_temp<-modRhum_temp[-c(16:21),]
modRhum_temp<-modRhum_temp[-13,]
modRhum_temp<-modRhum_temp[-12,] #removed all individuals with NAs for all measurements


#creating the MAX subfossil data set
subRfemur_set<-data.frame(subRfemur_temp)
names(subRfemur_set)<-c("Individual", "FFH_3DRF", "FHW_3DRF", "FHSA_3DRF", "FBE_3DRF", "DML2_3DRF")

subLfemur_set<-data.frame(subLfemur_temp)
names(subLfemur_set)<-c("Individual", "DSA_3DLF", "FHH_3DLF", "FHW_3DLF", "FHSA_3DLF", "FBE_3DLF", "DML1_3DLF", "DML2_3DLF")

subRhum_set<-data.frame(subRhum_temp)
names(subRhum_set)<-c("Individual", "VHD_3DRH", "HHW_3DRH", "HHSA_3DRH")

sub_femurs<-merge(subRfemur_set, subLfemur_set, by="Individual", all = TRUE)
MAX_subs<-merge(sub_femurs, subRhum_set, by="Individual", all = TRUE)


#subsample mods to same number of individuals as subs
modRfemur_indiv_sample<-sample(nrow(modRfemur_temp), 6, replace = FALSE)
modRfemur_indiv<-modRfemur_temp[modRfemur_indiv_sample,]

modLfemur_indiv_sample<-sample(nrow(modLfemur_temp), 5, replace = FALSE)
modLfemur_indiv<-modLfemur_temp[modLfemur_indiv_sample,]

modRhum_indiv_sample<-sample(nrow(modRhum_temp), 1, replace = FALSE)
modRhum_indiv<-modRhum_temp[modRhum_indiv_sample,]

#have to make modern subset look like subfossil
modRfemur_indiv[1,c(5:6)] = NA
modRfemur_indiv[2,c(2:5)] = NA
modRfemur_indiv[3,c(2:4)] = NA
modRfemur_indiv[4,c(2:4)] = NA
modRfemur_indiv[5,c(2:5)] = NA
modRfemur_indiv[6,c(5:6)] = NA

modLfemur_indiv[1,c(2:7)] = NA
modLfemur_indiv[2,c(2:6)] = NA
modLfemur_indiv[3,c(2, 6:8)] = NA
modLfemur_indiv[4,c(3:5)] = NA
modLfemur_indiv[5,c(2, 6:8)] = NA

#calculate geoMeans
library(EnvStats)
geoMeanFHH<-geoMean(modRfemur_indiv[,2], na.rm = TRUE)
geoMeanFHW<-geoMean(modRfemur_indiv[,3], na.rm = TRUE)
geoMeanFHSA<-geoMean(modRfemur_indiv[,4], na.rm = TRUE)
geoMeanFBE<-geoMean(modRfemur_indiv[,5], na.rm = TRUE)
geoMeanDML2<-geoMean(modRfemur_indiv[,6], na.rm = TRUE)
geoMeans<-c("geoMeans", geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML2)
modRfemur_indiv<-rbind(modRfemur_indiv, geoMeans)

library(EnvStats)
geoMeanDSA<-geoMean(modLfemur_indiv[,2], na.rm = TRUE)
geoMeanFHH<-geoMean(modLfemur_indiv[,3], na.rm = TRUE)
geoMeanFHW<-geoMean(modLfemur_indiv[,4], na.rm = TRUE)
geoMeanFHSA<-geoMean(modLfemur_indiv[,5], na.rm = TRUE)
geoMeanFBE<-geoMean(modLfemur_indiv[,6], na.rm = TRUE)
geoMeanDML1<-geoMean(modLfemur_indiv[,7], na.rm = TRUE)
geoMeanDML2<-geoMean(modLfemur_indiv[,8], na.rm = TRUE)
geoMeans<-c("geoMeans", geoMeanDSA, geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML1, geoMeanDML2)
modLfemur_indiv<-rbind(modLfemur_indiv, geoMeans)

library(EnvStats)
geoMeanVHD<-geoMean(modRhum_indiv[,2], na.rm = TRUE)
geoMeanHHW<-geoMean(modRhum_indiv[,3], na.rm = TRUE)
geoMeanHHSA<-geoMean(modRhum_indiv[,4], na.rm = TRUE)
geoMeans<-c("geoMeans", geoMeanVHD, geoMeanHHW, geoMeanHHSA)
modRhum_indiv<-rbind(modRhum_indiv, geoMeans)


#calculate fold changes
FHH_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,2]))/(as.numeric(modRfemur_indiv[7,2]))-1))
  FHH_foldchange_mod[a]<-c(z)
}

FHH_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,2]))/(as.numeric(modRfemur_indiv[7,2]))-1))
  FHH_foldchange_sub[a]<-c(z)
} 

FHW_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,3]))/(as.numeric(modRfemur_indiv[7,3]))-1))
  FHW_foldchange_mod[a]<-c(z)
}

FHW_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,3]))/(as.numeric(modRfemur_indiv[7,3]))-1))
  FHW_foldchange_sub[a]<-c(z)
} 

FHSA_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,4]))/(as.numeric(modRfemur_indiv[7,4]))-1))
  FHSA_foldchange_mod[a]<-c(z)
}

FHSA_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,4]))/(as.numeric(modRfemur_indiv[7,4]))-1))
  FHSA_foldchange_sub[a]<-c(z)
} 

FBE_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,5]))/(as.numeric(modRfemur_indiv[7,5]))-1))
  FBE_foldchange_mod[a]<-c(z)
}

FBE_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,5]))/(as.numeric(modRfemur_indiv[7,5]))-1))
  FBE_foldchange_sub[a]<-c(z)
} 

DML2_foldchange_mod<-c()
for (a in 1:nrow(modRfemur_indiv)){
  z<-(((as.numeric(modRfemur_indiv[a,6]))/(as.numeric(modRfemur_indiv[7,6]))-1))
  DML2_foldchange_mod[a]<-c(z)
}

DML2_foldchange_sub<-c()
for (a in 1:nrow(subRfemur_temp)){
  z<-(((as.numeric(subRfemur_temp[a,6]))/(as.numeric(modRfemur_indiv[7,6]))-1))
  DML2_foldchange_sub[a]<-c(z)
} 

Rfem_foldchange_mod<-cbind(modRfemur_indiv[,1], FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML2_foldchange_mod)
Rfem_foldchange_mod<-Rfem_foldchange_mod[-7,] #removes calculated Medians/geoMeans row
Rfem_foldchange_sub<-cbind(subRfemur_temp[,1], FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML2_foldchange_sub)
Rfem_foldchange<-rbind(Rfem_foldchange_mod, Rfem_foldchange_sub)
Rfem_pops<-c(rep("Modern", 6), rep("Subfossil", 6))
Rfem_foldchange<-cbind(Rfem_pops, Rfem_foldchange)
Rfem_foldchange<-data.frame(Rfem_foldchange)
names(Rfem_foldchange) = c("Population", "Individual", "FHH", "FHW", "FHSA", "FBE", "DML2")


DSA_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,2]))/(as.numeric(modLfemur_indiv[6,2]))-1))
  DSA_foldchange_mod[a]<-c(z)
}

DSA_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,2]))/(as.numeric(modLfemur_indiv[6,2]))-1))
  DSA_foldchange_sub[a]<-c(z)
}

FHH_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,3]))/(as.numeric(modLfemur_indiv[6,3]))-1))
  FHH_foldchange_mod[a]<-c(z)
}

FHH_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,3]))/(as.numeric(modLfemur_indiv[6,3]))-1))
  FHH_foldchange_sub[a]<-c(z)
} 

FHW_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,4]))/(as.numeric(modLfemur_indiv[6,4]))-1))
  FHW_foldchange_mod[a]<-c(z)
}

FHW_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,4]))/(as.numeric(modLfemur_indiv[6,4]))-1))
  FHW_foldchange_sub[a]<-c(z)
} 

FHSA_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,5]))/(as.numeric(modLfemur_indiv[6,5]))-1))
  FHSA_foldchange_mod[a]<-c(z)
}

FHSA_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,5]))/(as.numeric(modLfemur_indiv[6,5]))-1))
  FHSA_foldchange_sub[a]<-c(z)
} 

FBE_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,6]))/(as.numeric(modLfemur_indiv[6,6]))-1))
  FBE_foldchange_mod[a]<-c(z)
}

FBE_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,6]))/(as.numeric(modLfemur_indiv[6,6]))-1))
  FBE_foldchange_sub[a]<-c(z)
} 

DML1_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,7]))/(as.numeric(modLfemur_indiv[6,7]))-1))
  DML1_foldchange_mod[a]<-c(z)
}

DML1_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,7]))/(as.numeric(modLfemur_indiv[6,7]))-1))
  DML1_foldchange_sub[a]<-c(z)
} 

DML2_foldchange_mod<-c()
for (a in 1:nrow(modLfemur_indiv)){
  z<-(((as.numeric(modLfemur_indiv[a,8]))/(as.numeric(modLfemur_indiv[6,8]))-1))
  DML2_foldchange_mod[a]<-c(z)
}

DML2_foldchange_sub<-c()
for (a in 1:nrow(subLfemur_temp)){
  z<-(((as.numeric(subLfemur_temp[a,8]))/(as.numeric(modLfemur_indiv[6,8]))-1))
  DML2_foldchange_sub[a]<-c(z)
} 

Lfem_foldchange_mod<-cbind(modLfemur_indiv[,1], DSA_foldchange_mod, FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML1_foldchange_mod, DML2_foldchange_mod)
Lfem_foldchange_mod<-Lfem_foldchange_mod[-6,] #removes calculated Medians/geoMeans row
Lfem_foldchange_sub<-cbind(subLfemur_temp[,1], DSA_foldchange_sub, FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML1_foldchange_sub, DML2_foldchange_sub)
Lfem_foldchange<-rbind(Lfem_foldchange_mod, Lfem_foldchange_sub)
Lfem_pops<-c(rep("Modern", 5), rep("Subfossil", 5))
Lfem_foldchange<-cbind(Lfem_pops, Lfem_foldchange)
Lfem_foldchange<-data.frame(Lfem_foldchange)
names(Lfem_foldchange) = c("Population", "Individual", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2")


VHD_foldchange_mod<-c()
for (a in 1:nrow(modRhum_indiv)){
  z<-(((as.numeric(modRhum_indiv[a,2]))/(as.numeric(modRhum_indiv[2,2]))-1))
  VHD_foldchange_mod[a]<-c(z)
}

VHD_foldchange_sub<-c()
for (a in 1:nrow(subRhum_temp)){
  z<-(((as.numeric(subRhum_temp[a,2]))/(as.numeric(modRhum_indiv[2,2]))-1))
  VHD_foldchange_sub[a]<-c(z)
}

HHW_foldchange_mod<-c()
for (a in 1:nrow(modRhum_indiv)){
  z<-(((as.numeric(modRhum_indiv[a,3]))/(as.numeric(modRhum_indiv[2,3]))-1))
  HHW_foldchange_mod[a]<-c(z)
}

HHW_foldchange_sub<-c()
for (a in 1:nrow(subRhum_temp)){
  z<-(((as.numeric(subRhum_temp[a,3]))/(as.numeric(modRhum_indiv[2,3]))-1))
  HHW_foldchange_sub[a]<-c(z)
} 

HHSA_foldchange_mod<-c()
for (a in 1:nrow(modRhum_indiv)){
  z<-(((as.numeric(modRhum_indiv[a,4]))/(as.numeric(modRhum_indiv[2,4]))-1))
  HHSA_foldchange_mod[a]<-c(z)
}

HHSA_foldchange_sub<-c()
for (a in 1:nrow(subRhum_temp)){
  z<-(((as.numeric(subRhum_temp[a,4]))/(as.numeric(modRhum_indiv[2,4]))-1))
  HHSA_foldchange_sub[a]<-c(z)
}

Rhum_foldchange_mod<-cbind(modRhum_indiv[,1], VHD_foldchange_mod, HHW_foldchange_mod, HHSA_foldchange_mod)
Rhum_foldchange_mod<-Rhum_foldchange_mod[-2,] #removes calculated Medians/geoMeans row
Rhum_foldchange_sub<-cbind(subRhum_temp[,1], VHD_foldchange_sub, HHW_foldchange_sub, HHSA_foldchange_sub)
Rhum_foldchange<-rbind(Rhum_foldchange_mod, Rhum_foldchange_sub)
Rhum_pops<-c(rep("Modern", 1), rep("Subfossil", 1))
Rhum_foldchange<-cbind(Rhum_pops, Rhum_foldchange)
Rhum_foldchange<-data.frame(Rhum_foldchange)
names(Rhum_foldchange) = c("Population", "Individual", "VHD", "HHW", "HHSA")


#calculate averages
Rfemur<-Rfem_foldchange[,-2] #remove Individual
Rfemur<-Rfemur[,-1] #remove Population

Rfemur_averages<-c()
for (a in 1:nrow(Rfemur)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Rfemur[a,]), na.rm = TRUE) #arithmetic mean
  Rfemur_averages[a]<-c(z)
}
Rfemur<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
Rfemur<-data.frame(Rfemur)
names(Rfemur) = c("Population", "Individual", "Averages")


Lfemur<-Lfem_foldchange[,-2] #remove Individual
Lfemur<-Lfemur[,-1] #remove Population

Lfemur_averages<-c()
for (a in 1:nrow(Lfemur)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Lfemur[a,]), na.rm = TRUE) #arithmetic mean
  #library(EnvStats)
  #z<-geoMean(as.numeric(Rfemur[a,]), na.rm = TRUE) #geometric mean
  Lfemur_averages[a]<-c(z)
}
Lfemur<-cbind(Lfem_foldchange[,1], Lfem_foldchange[,2], Lfemur_averages)
Lfemur<-data.frame(Lfemur)
names(Lfemur) = c("Population", "Individual", "Averages")


Rhumerus<-Rhum_foldchange[,-2] #remove Individual
Rhumerus<-Rhumerus[,-1] #remove Population

Rhumerus_averages<-c()
for (a in 1:nrow(Rhumerus)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Rhumerus[a,]), na.rm = TRUE) #arithmetic mean
  #library(EnvStats)
  #z<-geoMean(as.numeric(Rfemur[a,]), na.rm = TRUE) #geometric mean
  Rhumerus_averages[a]<-c(z)
}
Rhumerus<-cbind(Rhum_foldchange[,1], Rhum_foldchange[,2], Rhumerus_averages)
Rhumerus<-data.frame(Rhumerus)
names(Rhumerus) = c("Population", "Individual", "Averages")


#creating the MAX modern data subset
MAX_mods<-rbind(Rfemur, Lfemur, Rhumerus)


for(n in 1:9999){
  modRfemur_indiv_sample<-sample(nrow(modRfemur_temp), 6, replace = FALSE)
  modRfemur_indiv<-modRfemur_temp[modRfemur_indiv_sample,]
  
  modLfemur_indiv_sample<-sample(nrow(modLfemur_temp), 5, replace = FALSE)
  modLfemur_indiv<-modLfemur_temp[modLfemur_indiv_sample,]
  
  modRhum_indiv_sample<-sample(nrow(modRhum_temp), 1, replace = FALSE)
  modRhum_indiv<-modRhum_temp[modRhum_indiv_sample,]
  
  #have to make modern subset look like subfossil
  modRfemur_indiv[1,c(5:6)] = NA
  modRfemur_indiv[2,c(2:5)] = NA
  modRfemur_indiv[3,c(2:4)] = NA
  modRfemur_indiv[4,c(2:4)] = NA
  modRfemur_indiv[5,c(2:5)] = NA
  modRfemur_indiv[6,c(5:6)] = NA
  
  modLfemur_indiv[1,c(2:7)] = NA
  modLfemur_indiv[2,c(2:6)] = NA
  modLfemur_indiv[3,c(2, 6:8)] = NA
  modLfemur_indiv[4,c(3:5)] = NA
  modLfemur_indiv[5,c(2, 6:8)] = NA
  
  #calculate geoMeans
  library(EnvStats)
  geoMeanFHH<-geoMean(modRfemur_indiv[,2], na.rm = TRUE)
  geoMeanFHW<-geoMean(modRfemur_indiv[,3], na.rm = TRUE)
  geoMeanFHSA<-geoMean(modRfemur_indiv[,4], na.rm = TRUE)
  geoMeanFBE<-geoMean(modRfemur_indiv[,5], na.rm = TRUE)
  geoMeanDML2<-geoMean(modRfemur_indiv[,6], na.rm = TRUE)
  geoMeans<-c("geoMeans", geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML2)
  modRfemur_indiv<-rbind(modRfemur_indiv, geoMeans)
  
  library(EnvStats)
  geoMeanDSA<-geoMean(modLfemur_indiv[,2], na.rm = TRUE)
  geoMeanFHH<-geoMean(modLfemur_indiv[,3], na.rm = TRUE)
  geoMeanFHW<-geoMean(modLfemur_indiv[,4], na.rm = TRUE)
  geoMeanFHSA<-geoMean(modLfemur_indiv[,5], na.rm = TRUE)
  geoMeanFBE<-geoMean(modLfemur_indiv[,6], na.rm = TRUE)
  geoMeanDML1<-geoMean(modLfemur_indiv[,7], na.rm = TRUE)
  geoMeanDML2<-geoMean(modLfemur_indiv[,8], na.rm = TRUE)
  geoMeans<-c("geoMeans", geoMeanDSA, geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML1, geoMeanDML2)
  modLfemur_indiv<-rbind(modLfemur_indiv, geoMeans)
  
  library(EnvStats)
  geoMeanVHD<-geoMean(modRhum_indiv[,2], na.rm = TRUE)
  geoMeanHHW<-geoMean(modRhum_indiv[,3], na.rm = TRUE)
  geoMeanHHSA<-geoMean(modRhum_indiv[,4], na.rm = TRUE)
  geoMeans<-c("geoMeans", geoMeanVHD, geoMeanHHW, geoMeanHHSA)
  modRhum_indiv<-rbind(modRhum_indiv, geoMeans)
  
  
  #calculate fold changes
  FHH_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,2]))/(as.numeric(modRfemur_indiv[7,2]))-1))
    FHH_foldchange_mod[a]<-c(z)
  }
  
  FHH_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,2]))/(as.numeric(modRfemur_indiv[7,2]))-1))
    FHH_foldchange_sub[a]<-c(z)
  } 
  
  FHW_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,3]))/(as.numeric(modRfemur_indiv[7,3]))-1))
    FHW_foldchange_mod[a]<-c(z)
  }
  
  FHW_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,3]))/(as.numeric(modRfemur_indiv[7,3]))-1))
    FHW_foldchange_sub[a]<-c(z)
  } 
  
  FHSA_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,4]))/(as.numeric(modRfemur_indiv[7,4]))-1))
    FHSA_foldchange_mod[a]<-c(z)
  }
  
  FHSA_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,4]))/(as.numeric(modRfemur_indiv[7,4]))-1))
    FHSA_foldchange_sub[a]<-c(z)
  } 
  
  FBE_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,5]))/(as.numeric(modRfemur_indiv[7,5]))-1))
    FBE_foldchange_mod[a]<-c(z)
  }
  
  FBE_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,5]))/(as.numeric(modRfemur_indiv[7,5]))-1))
    FBE_foldchange_sub[a]<-c(z)
  } 
  
  DML2_foldchange_mod<-c()
  for (a in 1:nrow(modRfemur_indiv)){
    z<-(((as.numeric(modRfemur_indiv[a,6]))/(as.numeric(modRfemur_indiv[7,6]))-1))
    DML2_foldchange_mod[a]<-c(z)
  }
  
  DML2_foldchange_sub<-c()
  for (a in 1:nrow(subRfemur_temp)){
    z<-(((as.numeric(subRfemur_temp[a,6]))/(as.numeric(modRfemur_indiv[7,6]))-1))
    DML2_foldchange_sub[a]<-c(z)
  } 
  
  Rfem_foldchange_mod<-cbind(modRfemur_indiv[,1], FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML2_foldchange_mod)
  Rfem_foldchange_mod<-Rfem_foldchange_mod[-7,] #removes calculated Medians/geoMeans row
  Rfem_foldchange_sub<-cbind(subRfemur_temp[,1], FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML2_foldchange_sub)
  Rfem_foldchange<-rbind(Rfem_foldchange_mod, Rfem_foldchange_sub)
  Rfem_pops<-c(rep("Modern", 6), rep("Subfossil", 6))
  Rfem_foldchange<-cbind(Rfem_pops, Rfem_foldchange)
  Rfem_foldchange<-data.frame(Rfem_foldchange)
  names(Rfem_foldchange) = c("Population", "Individual", "FHH", "FHW", "FHSA", "FBE", "DML2")
  
  
  DSA_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,2]))/(as.numeric(modLfemur_indiv[6,2]))-1))
    DSA_foldchange_mod[a]<-c(z)
  }
  
  DSA_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,2]))/(as.numeric(modLfemur_indiv[6,2]))-1))
    DSA_foldchange_sub[a]<-c(z)
  }
  
  FHH_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,3]))/(as.numeric(modLfemur_indiv[6,3]))-1))
    FHH_foldchange_mod[a]<-c(z)
  }
  
  FHH_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,3]))/(as.numeric(modLfemur_indiv[6,3]))-1))
    FHH_foldchange_sub[a]<-c(z)
  } 
  
  FHW_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,4]))/(as.numeric(modLfemur_indiv[6,4]))-1))
    FHW_foldchange_mod[a]<-c(z)
  }
  
  FHW_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,4]))/(as.numeric(modLfemur_indiv[6,4]))-1))
    FHW_foldchange_sub[a]<-c(z)
  } 
  
  FHSA_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,5]))/(as.numeric(modLfemur_indiv[6,5]))-1))
    FHSA_foldchange_mod[a]<-c(z)
  }
  
  FHSA_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,5]))/(as.numeric(modLfemur_indiv[6,5]))-1))
    FHSA_foldchange_sub[a]<-c(z)
  } 
  
  FBE_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,6]))/(as.numeric(modLfemur_indiv[6,6]))-1))
    FBE_foldchange_mod[a]<-c(z)
  }
  
  FBE_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,6]))/(as.numeric(modLfemur_indiv[6,6]))-1))
    FBE_foldchange_sub[a]<-c(z)
  } 
  
  DML1_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,7]))/(as.numeric(modLfemur_indiv[6,7]))-1))
    DML1_foldchange_mod[a]<-c(z)
  }
  
  DML1_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,7]))/(as.numeric(modLfemur_indiv[6,7]))-1))
    DML1_foldchange_sub[a]<-c(z)
  } 
  
  DML2_foldchange_mod<-c()
  for (a in 1:nrow(modLfemur_indiv)){
    z<-(((as.numeric(modLfemur_indiv[a,8]))/(as.numeric(modLfemur_indiv[6,8]))-1))
    DML2_foldchange_mod[a]<-c(z)
  }
  
  DML2_foldchange_sub<-c()
  for (a in 1:nrow(subLfemur_temp)){
    z<-(((as.numeric(subLfemur_temp[a,8]))/(as.numeric(modLfemur_indiv[6,8]))-1))
    DML2_foldchange_sub[a]<-c(z)
  } 
  
  Lfem_foldchange_mod<-cbind(modLfemur_indiv[,1], DSA_foldchange_mod, FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML1_foldchange_mod, DML2_foldchange_mod)
  Lfem_foldchange_mod<-Lfem_foldchange_mod[-6,] #removes calculated Medians/geoMeans row
  Lfem_foldchange_sub<-cbind(subLfemur_temp[,1], DSA_foldchange_sub, FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML1_foldchange_sub, DML2_foldchange_sub)
  Lfem_foldchange<-rbind(Lfem_foldchange_mod, Lfem_foldchange_sub)
  Lfem_pops<-c(rep("Modern", 5), rep("Subfossil", 5))
  Lfem_foldchange<-cbind(Lfem_pops, Lfem_foldchange)
  Lfem_foldchange<-data.frame(Lfem_foldchange)
  names(Lfem_foldchange) = c("Population", "Individual", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2")
  
  
  VHD_foldchange_mod<-c()
  for (a in 1:nrow(modRhum_indiv)){
    z<-(((as.numeric(modRhum_indiv[a,2]))/(as.numeric(modRhum_indiv[2,2]))-1))
    VHD_foldchange_mod[a]<-c(z)
  }
  
  VHD_foldchange_sub<-c()
  for (a in 1:nrow(subRhum_temp)){
    z<-(((as.numeric(subRhum_temp[a,2]))/(as.numeric(modRhum_indiv[2,2]))-1))
    VHD_foldchange_sub[a]<-c(z)
  }
  
  HHW_foldchange_mod<-c()
  for (a in 1:nrow(modRhum_indiv)){
    z<-(((as.numeric(modRhum_indiv[a,3]))/(as.numeric(modRhum_indiv[2,3]))-1))
    HHW_foldchange_mod[a]<-c(z)
  }
  
  HHW_foldchange_sub<-c()
  for (a in 1:nrow(subRhum_temp)){
    z<-(((as.numeric(subRhum_temp[a,3]))/(as.numeric(modRhum_indiv[2,3]))-1))
    HHW_foldchange_sub[a]<-c(z)
  } 
  
  HHSA_foldchange_mod<-c()
  for (a in 1:nrow(modRhum_indiv)){
    z<-(((as.numeric(modRhum_indiv[a,4]))/(as.numeric(modRhum_indiv[2,4]))-1))
    HHSA_foldchange_mod[a]<-c(z)
  }
  
  HHSA_foldchange_sub<-c()
  for (a in 1:nrow(subRhum_temp)){
    z<-(((as.numeric(subRhum_temp[a,4]))/(as.numeric(modRhum_indiv[2,4]))-1))
    HHSA_foldchange_sub[a]<-c(z)
  }
  
  Rhum_foldchange_mod<-cbind(modRhum_indiv[,1], VHD_foldchange_mod, HHW_foldchange_mod, HHSA_foldchange_mod)
  Rhum_foldchange_mod<-Rhum_foldchange_mod[-2,] #removes calculated Medians/geoMeans row
  Rhum_foldchange_sub<-cbind(subRhum_temp[,1], VHD_foldchange_sub, HHW_foldchange_sub, HHSA_foldchange_sub)
  Rhum_foldchange<-rbind(Rhum_foldchange_mod, Rhum_foldchange_sub)
  Rhum_pops<-c(rep("Modern", 1), rep("Subfossil", 1))
  Rhum_foldchange<-cbind(Rhum_pops, Rhum_foldchange)
  Rhum_foldchange<-data.frame(Rhum_foldchange)
  names(Rhum_foldchange) = c("Population", "Individual", "VHD", "HHW", "HHSA")
  
  
  #calculate averages
  Rfemur<-Rfem_foldchange[,-2] #remove Individual
  Rfemur<-Rfemur[,-1] #remove Population
  
  Rfemur_averages<-c()
  for (a in 1:nrow(Rfemur)){
    if ((a %% 1) == 0){print (a)}
    z<-mean(as.numeric(Rfemur[a,]), na.rm = TRUE) #arithmetic mean
    Rfemur_averages[a]<-c(z)
  }
  Rfemur<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
  Rfemur<-data.frame(Rfemur)
  names(Rfemur) = c("Population", "Individual", "Averages")
  
  
  Lfemur<-Lfem_foldchange[,-2] #remove Individual
  Lfemur<-Lfemur[,-1] #remove Population
  
  Lfemur_averages<-c()
  for (a in 1:nrow(Lfemur)){
    if ((a %% 1) == 0){print (a)}
    z<-mean(as.numeric(Lfemur[a,]), na.rm = TRUE) #arithmetic mean
    #library(EnvStats)
    #z<-geoMean(as.numeric(Rfemur[a,]), na.rm = TRUE) #geometric mean
    Lfemur_averages[a]<-c(z)
  }
  Lfemur<-cbind(Lfem_foldchange[,1], Lfem_foldchange[,2], Lfemur_averages)
  Lfemur<-data.frame(Lfemur)
  names(Lfemur) = c("Population", "Individual", "Averages")
  
  
  Rhumerus<-Rhum_foldchange[,-2] #remove Individual
  Rhumerus<-Rhumerus[,-1] #remove Population
  
  Rhumerus_averages<-c()
  for (a in 1:nrow(Rhumerus)){
    if ((a %% 1) == 0){print (a)}
    z<-mean(as.numeric(Rhumerus[a,]), na.rm = TRUE) #arithmetic mean
    #library(EnvStats)
    #z<-geoMean(as.numeric(Rfemur[a,]), na.rm = TRUE) #geometric mean
    Rhumerus_averages[a]<-c(z)
  }
  Rhumerus<-cbind(Rhum_foldchange[,1], Rhum_foldchange[,2], Rhumerus_averages)
  Rhumerus<-data.frame(Rhumerus)
  names(Rhumerus) = c("Population", "Individual", "Averages")

  MAX_mods_permute<-rbind(Rfemur, Lfemur, Rhumerus)
  MAX_mods<-rbind(MAX_mods, MAX_mods_permute)
  if ((n %% 1000) == 0){print (n)}
}


write.table(MAX_mods, "MAXPermute.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

##################################################
options(stringsAsFactors = FALSE)
MAX<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/MAXPermute.txt", quote="\"")

sub_permute<-subset(MAX, Population == "Subfossil")
mod_permute<-subset(MAX, Population == "Modern")

average_MAX_subfossil<-mean(as.numeric(sub_permute[,3]), na.rm = TRUE) #0.09091167
#average_MAX_subfossil from geoMeans script = #0.08924061

quantile(mod_permute[,3], na.rm = TRUE)
#           0%           25%           50%           75%          100% 
#-0.19147061     -0.02557775    0.00000000    0.02612512    0.22660732
quantile((mod_permute[,3]), na.rm = TRUE, c(0.05, 0.95))
#           5%         95% 
#  -0.07365117  0.08095605

#hist(Rfem_pVal_permute, main = "t-test p-values - Mod vs. Sub R Femur Agg Scores", ylim = c(0, 5000))
hist(mod_permute[,3], col=rgb(1,0,0,0.5), xlim=c(-0.3, 0.3), ylim=c(0, 35000), main = "MAX Mod 3D Permute", xlab="Aggregate Scores")
#hist(as.numeric(sub_permute[,3]), col=rgb(0,0,1,0.5), add=T)
#points(0.07179621, 0, pch=19, col="#cec3ae") #subfossil observed value
points(0.09091167, 0, pch=19, col="#cec3ae") #MAX subfossil observed value from permuted script
points(0.08924061, 0, pch=19, col="#cec3ae") #MAX subfossil observed value from geoMeans script
#points(0.0198988281883758, 0, pch=19, col="#88ccee") #TAO-66-26 AS
#points(0.020560842806777, 0, pch=19, col="#cc6677") #TAO-66-29 AS
#points(0.126896252573764, 0, pch=19, col="#ddcc77") #TAO-66-32 AS
#points(0.244276570389818, 0, pch=19, col="#3b7834") #TAO-66-33 AS
#points(0.00169789, 0, pch=19, col="red") #modern observed value
points(0.08095605, 0, pch=19, col="black") #95th percentile


#average_MNI_subfossil = 0.1029081
#average_modern = 0.008782976
#average_MNI_modern = 0.00169789
#TAO-66-26 AS = 0.0198988281883758
#TAO-66-29 AS = 0.020560842806777
#TAO-66-32 AS = 0.126896252573764
#TAO-66-33 AS = 0.244276570389818

prob_great_R<-which(mod_permute[,3] >= 0.09091167) #4060/120000 = 0.03383333
prob_less_R<-which(mod_permute[,3] < 0.09091167) #106886/120000 = 0.8907167
prob_NA_R<-which(mod_permute[,3] == "NaN") #9054/40000 = 0.22635

prob_great_R<-which(mod_permute[,3] >= 0.08924061) #geoMeans script: 4423/120000 = 0.03685833
prob_less_R<-which(mod_permute[,3] < 0.08924061) #geoMeans script: 106523/120000 = 0.8876917
prob_NA_R<-which(mod_permute[,3] == "NaN") #geoMeans script: 9054/40000 = 0.22635

a<-mean(mod_permute[,3], na.rm = TRUE)
s<-sd(mod_permute[,3], na.rm = TRUE)
n<-nrow(mod_permute)
error<-qnorm(0.975)*s/sqrt(n)
left<-a-error #0.001031902
right<-a+error #0.001550119
