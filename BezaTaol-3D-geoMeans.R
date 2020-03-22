#geoMeans
setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##R femur##           ##L femur starts line 230##
modRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Mod.csv", quote="\"")
modRfemur[31,1] = "BM 550" #renames individual 550 as BM 550
modRfemur[30,1] = "BM 570" #renames individual 570 as BM 570
subRfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRFemur_Sub.csv", quote="\"")
subRfemur<-subRfemur[-14,] #remove APS-2
subRfemur<-subRfemur[-13,] #remove TAO-66-22
subRfemur<-subRfemur[-3,] #remove TAO-66-28, only dated R femurs left

shapiro.test(modRfemur[,2]) #W = 0.95746, p-value = 0.4945 - MFL distribution probably normal
shapiro.test(modRfemur[,3]) #W = 0.93441, p-value = 0.1362 - MFD distribution probably normal
shapiro.test(modRfemur[,4]) #W = 0.96896, p-value = 0.6641 - FMD distribution probably normal
shapiro.test(modRfemur[,5]) #W = 0.94042, p-value = 0.1835 - DSA distribution probably normal
shapiro.test(modRfemur[,6]) #W = 0.93164, p-value = 0.1061 - FHH distribution probably normal
shapiro.test(modRfemur[,7]) #W = 0.94076, p-value = 0.1696 - FHW distribution probably normal
shapiro.test(modRfemur[,8]) #W = 0.94693, p-value = 0.2323 - FHSA distribution probably normal
shapiro.test(modRfemur[,9]) #W = 0.95149, p-value = 0.3142 - FBE distribution probably normal
shapiro.test(modRfemur[,10]) #W = 0.97297, p-value = 0.7595 - DML1 distribution probably normal
shapiro.test(modRfemur[,11]) #W = 0.95038, p-value = 0.298 - DML2 distribution probably normal

library(EnvStats)
geoMeanMFL<-geoMean(modRfemur[,2], na.rm = TRUE) #177.255
geoMeanMFD<-geoMean(modRfemur[,3], na.rm = TRUE) #8.725
geoMeanFMD<-geoMean(modRfemur[,4], na.rm = TRUE) #9.212
geoMeanDSA<-geoMean(modRfemur[,5], na.rm = TRUE) #623.651
geoMeanFHH<-geoMean(modRfemur[,6], na.rm = TRUE) #12.365
geoMeanFHW<-geoMean(modRfemur[,7], na.rm = TRUE) #12.216
geoMeanFHSA<-geoMean(modRfemur[,8], na.rm = TRUE) #347.593
geoMeanFBE<-geoMean(modRfemur[,9], na.rm = TRUE) #19.837
geoMeanDML1<-geoMean(modRfemur[,10], na.rm = TRUE) #4.987
geoMeanDML2<-geoMean(modRfemur[,11], na.rm = TRUE) #6.863

geoMeans<-c("geoMeans", geoMeanMFL, geoMeanMFD, geoMeanFMD, geoMeanDSA, geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML1, geoMeanDML2)
modRfemur<-rbind(modRfemur, geoMeans)

MFL_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,2]))/(as.numeric(modRfemur[32,2]))-1))
  MFL_foldchange_mod[a]<-c(z)
}

MFL_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,2]))/(as.numeric(modRfemur[32,2]))-1))
  MFL_foldchange_sub[a]<-c(z)
}

MFD_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,3]))/(as.numeric(modRfemur[32,3]))-1))
  MFD_foldchange_mod[a]<-c(z)
}

MFD_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,3]))/(as.numeric(modRfemur[32,3]))-1))
  MFD_foldchange_sub[a]<-c(z)
}

FMD_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,4]))/(as.numeric(modRfemur[32,4]))-1))
  FMD_foldchange_mod[a]<-c(z)
}

FMD_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,4]))/(as.numeric(modRfemur[32,4]))-1))
  FMD_foldchange_sub[a]<-c(z)
}

DSA_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,5]))/(as.numeric(modRfemur[32,5]))-1))
  DSA_foldchange_mod[a]<-c(z)
}

DSA_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,5]))/(as.numeric(modRfemur[32,5]))-1))
  DSA_foldchange_sub[a]<-c(z)
}

FHH_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,6]))/(as.numeric(modRfemur[32,6]))-1))
  FHH_foldchange_mod[a]<-c(z)
}

FHH_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,6]))/(as.numeric(modRfemur[32,6]))-1))
  FHH_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHH_foldchange_mod), as.numeric(FHH_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.001011305   mean of y = 0.072391179
#Test Statistic:                  t = -2.336255
#Test Statistic Parameter:        df = 1.212238
#P-value:                         0.2209746
#95% Confidence Interval:         LCL = -0.3306850   UCL =  0.1879252

FHW_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,7]))/(as.numeric(modRfemur[32,7]))-1))
  FHW_foldchange_mod[a]<-c(z)
}

FHW_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,7]))/(as.numeric(modRfemur[32,7]))-1))
  FHW_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHW_foldchange_mod), as.numeric(FHW_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.0009107089   mean of y = 0.0772468165
#Test Statistic:                  t = -5.054087
#Test Statistic Parameter:        df = 2.264688
#P-value:                         0.02845659
#95% Confidence Interval:         LCL = -0.13456023   UCL =  -0.01811199

FHSA_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,8]))/(as.numeric(modRfemur[32,8]))-1))
  FHSA_foldchange_mod[a]<-c(z)
}

FHSA_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,8]))/(as.numeric(modRfemur[32,8]))-1))
  FHSA_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHSA_foldchange_mod), as.numeric(FHSA_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.00481067   mean of y = 0.12402956
#Test Statistic:                  t = -1.071176
#Test Statistic Parameter:        df = 1.070205
#P-value:                         0.4687673
#95% Confidence Interval:         LCL = -1.331945   UCL =  1.093507

FBE_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,9]))/(as.numeric(modRfemur[32,9]))-1))
  FBE_foldchange_mod[a]<-c(z)
}

FBE_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,9]))/(as.numeric(modRfemur[32,9]))-1))
  FBE_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FBE_foldchange_mod), as.numeric(FBE_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.0006972242   mean of y = 0.0182833927
#Test Statistic:                  t = -0.3317484
#Test Statistic Parameter:        df = 1.045292
#P-value:                         0.7942047
#95% Confidence Interval:         LCL = -0.6259493   UCL =  0.5907769

DML1_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,10]))/(as.numeric(modRfemur[32,10]))-1))
  DML1_foldchange_mod[a]<-c(z)
}

DML1_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,10]))/(as.numeric(modRfemur[32,10]))-1))
  DML1_foldchange_sub[a]<-c(z)
}

DML2_foldchange_mod<-c()
for (a in 1:nrow(modRfemur)){
  z<-(((as.numeric(modRfemur[a,11]))/(as.numeric(modRfemur[32,11]))-1))
  DML2_foldchange_mod[a]<-c(z)
}

DML2_foldchange_sub<-c()
for (a in 1:nrow(subRfemur)){
  z<-(((as.numeric(subRfemur[a,11]))/(as.numeric(modRfemur[32,11]))-1))
  DML2_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(DML2_foldchange_mod), as.numeric(DML2_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.002557063   mean of y = 0.130630701
#Test Statistic:                  t = -2.416837
#Test Statistic Parameter:        df = 3.54291
#P-value:                         0.08127933
#95% Confidence Interval:         LCL = -0.28300089   UCL =  0.02685361


Rfem_foldchange_mod<-cbind(modRfemur[,1], MFL_foldchange_mod, MFD_foldchange_mod, FMD_foldchange_mod, DSA_foldchange_mod, FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML1_foldchange_mod, DML2_foldchange_mod)
Rfem_foldchange_sub<-cbind(subRfemur[,1], MFL_foldchange_sub, MFD_foldchange_sub, FMD_foldchange_sub, DSA_foldchange_sub, FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML1_foldchange_sub, DML2_foldchange_sub)
Rfem_foldchange<-rbind(Rfem_foldchange_mod, Rfem_foldchange_sub)
Rfem_foldchange<-Rfem_foldchange[-32,] #remove "geoMeans" row
Rfem_pops<-c(rep("Modern", 31), rep("Subfossil", 11))
Rfem_foldchange<-cbind(Rfem_pops, Rfem_foldchange)
Rfem_foldchange<-data.frame(Rfem_foldchange)
names(Rfem_foldchange) = c("Population", "Individual", "MFL", "MFD", "FMD", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2")

library(reshape)
Rfem_foldchange.df<-melt(Rfem_foldchange, id.vars = c("Population", "Individual"), measure.vars = c("MFL", "MFD", "FMD", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2"))
Rfem_foldchange.df<-rename(Rfem_foldchange.df, c(X1="Population", X2="measurement"))
library(ggplot2)
noAvizo<-Rfem_foldchange.df[c(1:126, 169:252, 295:420),]
ggplot(Rfem_foldchange.df, aes(x=Population, y=as.numeric(value), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Fold-Change (based on modern population geoMean)") +
  ggtitle("All Right Femoral Measurements") + geom_jitter(width = 0.1, na.rm = TRUE) + facet_wrap(~variable)

#ggsave("~/Desktop/Rfem.pdf")

write.table(Rfem_foldchange, "3D_Rfem_foldchange_geoMeans.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#



setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##L femur##           ##R humerus starts line 455##
modLfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLFemur_Mod.csv", quote="\"")
modLfemur[31,1] = "BM 550" #renames individual 550 as BM 550
modLfemur[30,1] = "BM 570" #renames individual 570 as BM 570
subLfemur<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLFemur_Sub.csv", quote="\"")
subLfemur<-subLfemur[-14,] #remove APS-2
subLfemur<-subLfemur[-13,] #remove TAO-66-22
subLfemur<-subLfemur[-3,] #remove TAO-66-28, only dated L femurs left

shapiro.test(modLfemur[,2]) #W = 0.91113, p-value = 0.07764 - MFL distribution probably normal
shapiro.test(modLfemur[,3]) #W = 0.93596, p-value = 0.1325 - MFD distribution probably normal
shapiro.test(modLfemur[,4]) #W = 0.92022, p-value = 0.05906 - FMD distribution probably normal
shapiro.test(modLfemur[,5]) #W = 0.93445, p-value = 0.169 - DSA distribution probably normal
shapiro.test(modLfemur[,6]) #W = 0.93824, p-value = 0.2221 - FHH distribution probably normal
shapiro.test(modLfemur[,7]) #W = 0.93863, p-value = 0.2259 - FHW distribution probably normal
shapiro.test(modLfemur[,8]) #W = 0.95133, p-value = 0.3877 - FHSA distribution probably normal
shapiro.test(modLfemur[,9]) #W = 0.96878, p-value = 0.7058 - FBE distribution probably normal
shapiro.test(modLfemur[,10]) #W = 0.94131, p-value = 0.2313 - DML1 distribution probably normal
shapiro.test(modLfemur[,11]) #W = 0.94915, p-value = 0.3283 - DML2 distribution probably normal

library(EnvStats)
geoMeanMFL<-geoMean(modLfemur[,2], na.rm = TRUE) #176.825
geoMeanMFD<-geoMean(modLfemur[,3], na.rm = TRUE) #8.549
geoMeanFMD<-geoMean(modLfemur[,4], na.rm = TRUE) #9.293
geoMeanDSA<-geoMean(modLfemur[,5], na.rm = TRUE) #565.141
geoMeanFHH<-geoMean(modLfemur[,6], na.rm = TRUE) #12.408
geoMeanFHW<-geoMean(modLfemur[,7], na.rm = TRUE) #12.357
geoMeanFHSA<-geoMean(modLfemur[,8], na.rm = TRUE) #347.279
geoMeanFBE<-geoMean(modLfemur[,9], na.rm = TRUE) #19.772
geoMeanDML1<-geoMean(modLfemur[,10], na.rm = TRUE) #4.521
geoMeanDML2<-geoMean(modLfemur[,11], na.rm = TRUE) #6.756

geoMeans<-c("geoMeans", geoMeanMFL, geoMeanMFD, geoMeanFMD, geoMeanDSA, geoMeanFHH, geoMeanFHW, geoMeanFHSA, geoMeanFBE, geoMeanDML1, geoMeanDML2)
modLfemur<-rbind(modLfemur, geoMeans)

MFL_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,2]))/(as.numeric(modLfemur[32,2]))-1))
  MFL_foldchange_mod[a]<-c(z)
}

MFL_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,2]))/(as.numeric(modLfemur[32,2]))-1))
  MFL_foldchange_sub[a]<-c(z)
}

MFD_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,3]))/(as.numeric(modLfemur[32,3]))-1))
  MFD_foldchange_mod[a]<-c(z)
}

MFD_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,3]))/(as.numeric(modLfemur[32,3]))-1))
  MFD_foldchange_sub[a]<-c(z)
}

FMD_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,4]))/(as.numeric(modLfemur[32,4]))-1))
  FMD_foldchange_mod[a]<-c(z)
}

FMD_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,4]))/(as.numeric(modLfemur[32,4]))-1))
  FMD_foldchange_sub[a]<-c(z)
}

DSA_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,5]))/(as.numeric(modLfemur[32,5]))-1))
  DSA_foldchange_mod[a]<-c(z)
}

DSA_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,5]))/(as.numeric(modLfemur[32,5]))-1))
  DSA_foldchange_sub[a]<-c(z)
}

#t.test(as.numeric(DSA_foldchange_mod), as.numeric(DSA_foldchange_sub))

FHH_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,6]))/(as.numeric(modLfemur[32,6]))-1))
  FHH_foldchange_mod[a]<-c(z)
}

FHH_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,6]))/(as.numeric(modLfemur[32,6]))-1))
  FHH_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHH_foldchange_mod), as.numeric(FHH_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.001115396   mean of y = 0.048922739
#Test Statistic:                  t = -4.410532
#Test Statistic Parameter:        df = 20.92123
#P-value:                         0.0002453617
#95% Confidence Interval:         LCL = -0.07035419   UCL =  -0.02526049

FHW_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,7]))/(as.numeric(modLfemur[32,7]))-1))
  FHW_foldchange_mod[a]<-c(z)
}

FHW_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,7]))/(as.numeric(modLfemur[32,7]))-1))
  FHW_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHW_foldchange_mod), as.numeric(FHW_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.001146035   mean of y = 0.076301679
#Test Statistic:                  t = -2.361401
#Test Statistic Parameter:        df = 1.275455
#P-value:                         0.2094836
#95% Confidence Interval:         LCL = -0.3217819   UCL =  0.1714706 

FHSA_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,8]))/(as.numeric(modLfemur[32,8]))-1))
  FHSA_foldchange_mod[a]<-c(z)
}

FHSA_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,8]))/(as.numeric(modLfemur[32,8]))-1))
  FHSA_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(FHSA_foldchange_mod), as.numeric(FHSA_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.005399407   mean of y = 0.062518539
#Test Statistic:                  t = -0.9006759
#Test Statistic Parameter:        df = 1.340897
#P-value:                         0.5000872
#95% Confidence Interval:         LCL = -0.5090145   UCL =  0.3947762

FBE_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,9]))/(as.numeric(modLfemur[32,9]))-1))
  FBE_foldchange_mod[a]<-c(z)
}

FBE_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,9]))/(as.numeric(modLfemur[32,9]))-1))
  FBE_foldchange_sub[a]<-c(z)
}

#t.test(as.numeric(FBE_foldchange_mod), as.numeric(FBE_foldchange_sub))


DML1_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,10]))/(as.numeric(modLfemur[32,10]))-1))
  DML1_foldchange_mod[a]<-c(z)
}

DML1_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,10]))/(as.numeric(modLfemur[32,10]))-1))
  DML1_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(DML1_foldchange_mod), as.numeric(DML1_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.002632843   mean of y = 0.165583992
#Test Statistic:                  t = -8.919964
#Test Statistic Parameter:        df = 12.06183
#P-value:                         1.167918e-06
#95% Confidence Interval:         LCL = -0.2027314   UCL =  -0.1231709

DML2_foldchange_mod<-c()
for (a in 1:nrow(modLfemur)){
  z<-(((as.numeric(modLfemur[a,11]))/(as.numeric(modLfemur[32,11]))-1))
  DML2_foldchange_mod[a]<-c(z)
}

DML2_foldchange_sub<-c()
for (a in 1:nrow(subLfemur)){
  z<-(((as.numeric(subLfemur[a,11]))/(as.numeric(modLfemur[32,11]))-1))
  DML2_foldchange_sub[a]<-c(z)
}

t.test(as.numeric(DML2_foldchange_mod), as.numeric(DML2_foldchange_sub))
#Estimated Parameter(s):          mean of x = 0.002626718   mean of y = -0.003380161
#Test Statistic:                  t = 0.06747748
#Test Statistic Parameter:        df = 2.133477
#P-value:                         0.9519974
#95% Confidence Interval:         LCL = -0.3549458   UCL =  0.3669596 


Lfem_foldchange_mod<-cbind(modLfemur[,1], MFL_foldchange_mod, MFD_foldchange_mod, FMD_foldchange_mod, DSA_foldchange_mod, FHH_foldchange_mod, FHW_foldchange_mod, FHSA_foldchange_mod, FBE_foldchange_mod, DML1_foldchange_mod, DML2_foldchange_mod)
Lfem_foldchange_sub<-cbind(subLfemur[,1], MFL_foldchange_sub, MFD_foldchange_sub, FMD_foldchange_sub, DSA_foldchange_sub, FHH_foldchange_sub, FHW_foldchange_sub, FHSA_foldchange_sub, FBE_foldchange_sub, DML1_foldchange_sub, DML2_foldchange_sub)
Lfem_foldchange<-rbind(Lfem_foldchange_mod, Lfem_foldchange_sub)
Lfem_foldchange<-Lfem_foldchange[-32,] #remove "geoMeans" row
Lfem_pops<-c(rep("Modern", 31), rep("Subfossil", 11))
Lfem_foldchange<-cbind(Lfem_pops, Lfem_foldchange)
Lfem_foldchange<-data.frame(Lfem_foldchange)
names(Lfem_foldchange) = c("Population", "Individual", "MFL", "MFD", "FMD", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2")

library(reshape)
Lfem_foldchange.df<-melt(Lfem_foldchange, id.vars = c("Population", "Individual"), measure.vars = c("MFL", "MFD", "FMD", "DSA", "FHH", "FHW", "FHSA", "FBE", "DML1", "DML2"))
Lfem_foldchange.df<-rename(Lfem_foldchange.df, c(X1="Population", X2="measurement"))
library(ggplot2)
noAvizo<-Lfem_foldchange.df[c(1:126, 169:252, 295:420),]
ggplot(Lfem_foldchange.df, aes(x=Population, y=as.numeric(value), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Fold-Change (based on modern population geoMean)") +
  ggtitle("All Left Femoral Measurements") + geom_jitter(width = 0.1, na.rm = TRUE) + facet_wrap(~variable)

write.table(Lfem_foldchange, "3D_Lfem_foldchange_geoMeans.txt", sep = "\t", row.names = TRUE, col.names = TRUE)




setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##R humerus##           ##L humerus starts line 595## 
modRhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRHumerus_Mod.csv", quote="\"")
modRhum[31,1] = "BM 550" #renames individual 550 as BM 550
modRhum[30,1] = "BM 570" #renames individual 570 as BM 570
subRhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DRHumerus_Sub.csv", quote="\"")

shapiro.test(modRhum[,2]) #W = 0.96293, p-value = 0.6591 - MHL distribution probably normal
shapiro.test(modRhum[,3]) #W = 0.97647, p-value = 0.8382 - MHD1 distribution probably normal
shapiro.test(modRhum[,4]) #W = 0.96555, p-value = 0.5839 - MHD2 distribution probably normal
shapiro.test(modRhum[,5]) #W = 0.93748, p-value = 0.2148 - VHD distribution probably normal
shapiro.test(modRhum[,6]) #W = 0.92199, p-value = 0.1082 - HHW distribution probably normal
shapiro.test(modRhum[,7]) #W = 0.90744, p-value = 0.04199 - BBH distribution probably NOT normal
shapiro.test(modRhum[,8]) #W = 0.97709, p-value = 0.8913 - HHSA distribution probably normal

library(EnvStats)
geoMeanMHL<-geoMean(modRhum[,2], na.rm = TRUE) #91.228
geoMeanMHD1<-geoMean(modRhum[,3], na.rm = TRUE) #7.420
geoMeanMHD2<-geoMean(modRhum[,4], na.rm = TRUE) #6.116
geoMeanVHD<-geoMean(modRhum[,5], na.rm = TRUE) #11.292
geoMeanHHW<-geoMean(modRhum[,6], na.rm = TRUE) #9.207
geoMeanBBH<-geoMean(modRhum[,7], na.rm = TRUE) #20.922
geoMeanHHSA<-geoMean(modRhum[,8], na.rm = TRUE) #173.149

geoMeans<-c("geoMeans", geoMeanMHL, geoMeanMHD1, geoMeanMHD2, geoMeanVHD, geoMeanHHW, geoMeanBBH, geoMeanHHSA)
modRhum<-rbind(modRhum, geoMeans)

MHL_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,2]))/(as.numeric(modRhum[32,2]))-1))
  MHL_foldchange_mod[a]<-c(z)
}

MHL_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,2]))/(as.numeric(modRhum[32,2]))-1))
  MHL_foldchange_sub[a]<-c(z)
}

MHD1_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,3]))/(as.numeric(modRhum[32,3]))-1))
  MHD1_foldchange_mod[a]<-c(z)
}

MHD1_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,3]))/(as.numeric(modRhum[32,3]))-1))
  MHD1_foldchange_sub[a]<-c(z)
}

MHD2_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,4]))/(as.numeric(modRhum[32,4]))-1))
  MHD2_foldchange_mod[a]<-c(z)
}

MHD2_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,4]))/(as.numeric(modRhum[32,4]))-1))
  MHD2_foldchange_sub[a]<-c(z)
}

VHD_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,5]))/(as.numeric(modRhum[32,5]))-1))
  VHD_foldchange_mod[a]<-c(z)
}

VHD_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,5]))/(as.numeric(modRhum[32,5]))-1))
  VHD_foldchange_sub[a]<-c(z)
}

HHW_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,6]))/(as.numeric(modRhum[32,6]))-1))
  HHW_foldchange_mod[a]<-c(z)
}

HHW_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,6]))/(as.numeric(modRhum[32,6]))-1))
  HHW_foldchange_sub[a]<-c(z)
}

BBH_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,7]))/(as.numeric(modRhum[32,7]))-1))
  BBH_foldchange_mod[a]<-c(z)
}

BBH_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,7]))/(as.numeric(modRhum[32,7]))-1))
  BBH_foldchange_sub[a]<-c(z)
}

HHSA_foldchange_mod<-c()
for (a in 1:nrow(modRhum)){
  z<-(((as.numeric(modRhum[a,8]))/(as.numeric(modRhum[32,8]))-1))
  HHSA_foldchange_mod[a]<-c(z)
}

HHSA_foldchange_sub<-c()
for (a in 1:nrow(subRhum)){
  z<-(((as.numeric(subRhum[a,8]))/(as.numeric(modRhum[32,8]))-1))
  HHSA_foldchange_sub[a]<-c(z)
}


Rhum_foldchange_mod<-cbind(modRhum[,1], MHL_foldchange_mod, MHD1_foldchange_mod, MHD2_foldchange_mod, VHD_foldchange_mod, HHW_foldchange_mod, BBH_foldchange_mod, HHSA_foldchange_mod)
Rhum_foldchange_sub<-cbind(subRhum[,1], MHL_foldchange_sub, MHD1_foldchange_sub, MHD2_foldchange_sub, VHD_foldchange_sub, HHW_foldchange_sub, BBH_foldchange_sub, HHSA_foldchange_sub)
Rhum_foldchange<-rbind(Rhum_foldchange_mod, Rhum_foldchange_sub)
Rhum_foldchange<-Rhum_foldchange[-32,] #remove "geoMeans" row
Rhum_pops<-c(rep("Modern", 31), rep("Subfossil", 2))
Rhum_foldchange<-cbind(Rhum_pops, Rhum_foldchange)
Rhum_foldchange<-data.frame(Rhum_foldchange)
names(Rhum_foldchange) = c("Population", "Individual", "MHL", "MHD1", "MHD2", "VHD", "HHW", "BBH", "HHSA")

library(reshape)
Rhum_foldchange.df<-melt(Rhum_foldchange, id.vars = c("Population", "Individual"), measure.vars = c("MHL", "MHD1", "MHD2", "VHD", "HHW", "BBH", "HHSA"))
Rhum_foldchange.df<-rename(Rhum_foldchange.df, c(X1="Population", X2="measurement"))
library(ggplot2)
noAvizo<-Rhum_foldchange.df[1:198,]
ggplot(Rhum_foldchange.df, aes(x=Population, y=as.numeric(value), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Fold-Change (based on modern population geoMean)") +
  ggtitle("All Right Humeral Measurements") + geom_jitter(width = 0.1, na.rm = TRUE) + facet_wrap(~variable)

write.table(Rhum_foldchange, "3D_Rhum_foldchange_geoMeans.txt", sep = "\t", row.names = TRUE, col.names = TRUE)




setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##L humerus##           ##All elements, all sides starts line 455##
modLhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLHumerus_Mod.csv", quote="\"")
modLhum[31,1] = "BM 550" #renames individual 550 as BM 550
modLhum[30,1] = "BM 570" #renames individual 570 as BM 570
subLhum<-read.csv("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3DLHumerus_Sub.csv", quote="\"")

shapiro.test(modLhum[,2]) #W = 0.9613, p-value = 0.5982 - MHL distribution probably normal
shapiro.test(modLhum[,3]) #W = 0.93428, p-value = 0.1866 - MHD1 distribution probably normal
shapiro.test(modLhum[,4]) #W = 0.94034, p-value = 0.2433 - MHD2 distribution probably normal
shapiro.test(modLhum[,5]) #W = 0.97974, p-value = 0.9307 - VHD distribution probably normal
shapiro.test(modLhum[,6]) #W = 0.95468, p-value = 0.4438 - HHW distribution probably normal
shapiro.test(modLhum[,7]) #W = 0.8815, p-value = 0.01883 - BBH distribution probably NOT normal
shapiro.test(modLhum[,8]) #W = 0.96038, p-value = 0.5515 - HHSA distribution probably normal

library(EnvStats)
geoMeanMHL<-geoMean(modLhum[,2], na.rm = TRUE) #90.788
geoMeanMHD1<-geoMean(modLhum[,3], na.rm = TRUE) #7.441
geoMeanMHD2<-geoMean(modLhum[,4], na.rm = TRUE) #6.112
geoMeanVHD<-geoMean(modLhum[,5], na.rm = TRUE) #11.412
geoMeanHHW<-geoMean(modLhum[,6], na.rm = TRUE) #9.193
geoMeanBBH<-geoMean(modLhum[,7], na.rm = TRUE) #20.896
geoMeanHHSA<-geoMean(modLhum[,8], na.rm = TRUE) #177.670

geoMeans<-c("geoMeans", geoMeanMHL, geoMeanMHD1, geoMeanMHD2, geoMeanVHD, geoMeanHHW, geoMeanBBH, geoMeanHHSA)
modLhum<-rbind(modLhum, geoMeans)

MHL_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,2]))/(as.numeric(modLhum[32,2]))-1))
  MHL_foldchange_mod[a]<-c(z)
}

MHL_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,2]))/(as.numeric(modLhum[32,2]))-1))
  MHL_foldchange_sub[a]<-c(z)
}

MHD1_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,3]))/(as.numeric(modLhum[32,3]))-1))
  MHD1_foldchange_mod[a]<-c(z)
}

MHD1_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,3]))/(as.numeric(modLhum[32,3]))-1))
  MHD1_foldchange_sub[a]<-c(z)
}

MHD2_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,4]))/(as.numeric(modLhum[32,4]))-1))
  MHD2_foldchange_mod[a]<-c(z)
}

MHD2_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,4]))/(as.numeric(modLhum[32,4]))-1))
  MHD2_foldchange_sub[a]<-c(z)
}

VHD_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,5]))/(as.numeric(modLhum[32,5]))-1))
  VHD_foldchange_mod[a]<-c(z)
}

VHD_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,5]))/(as.numeric(modLhum[32,5]))-1))
  VHD_foldchange_sub[a]<-c(z)
}

HHW_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,6]))/(as.numeric(modLhum[32,6]))-1))
  HHW_foldchange_mod[a]<-c(z)
}

HHW_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,6]))/(as.numeric(modLhum[32,6]))-1))
  HHW_foldchange_sub[a]<-c(z)
}

BBH_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,7]))/(as.numeric(modLhum[32,7]))-1))
  BBH_foldchange_mod[a]<-c(z)
}

BBH_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,7]))/(as.numeric(modLhum[32,7]))-1))
  BBH_foldchange_sub[a]<-c(z)
}

HHSA_foldchange_mod<-c()
for (a in 1:nrow(modLhum)){
  z<-(((as.numeric(modLhum[a,8]))/(as.numeric(modLhum[32,8]))-1))
  HHSA_foldchange_mod[a]<-c(z)
}

HHSA_foldchange_sub<-c()
for (a in 1:nrow(subLhum)){
  z<-(((as.numeric(subLhum[a,8]))/(as.numeric(modLhum[32,8]))-1))
  HHSA_foldchange_sub[a]<-c(z)
}


Lhum_foldchange_mod<-cbind(modLhum[,1], MHL_foldchange_mod, MHD1_foldchange_mod, MHD2_foldchange_mod, VHD_foldchange_mod, HHW_foldchange_mod, BBH_foldchange_mod, HHSA_foldchange_mod)
Lhum_foldchange_sub<-cbind(subLhum[,1], MHL_foldchange_sub, MHD1_foldchange_sub, MHD2_foldchange_sub, VHD_foldchange_sub, HHW_foldchange_sub, BBH_foldchange_sub, HHSA_foldchange_sub)
Lhum_foldchange<-Lhum_foldchange_mod
Lhum_foldchange<-Lhum_foldchange[-32,] #remove "geoMeans" row
Lhum_pops<-c(rep("Modern", 31))
Lhum_foldchange<-cbind(Lhum_pops, Lhum_foldchange)
Lhum_foldchange<-data.frame(Lhum_foldchange)
names(Lhum_foldchange) = c("Population", "Individual", "MHL", "MHD1", "MHD2", "VHD", "HHW", "BBH", "HHSA")

library(reshape)
Lhum_foldchange.df<-melt(Lhum_foldchange, id.vars = c("Population", "Individual"), measure.vars = c("MHL", "MHD1", "MHD2", "VHD", "HHW", "BBH", "HHSA"))
Lhum_foldchange.df<-rename(Lhum_foldchange.df, c(X1="Population", X2="measurement"))
library(ggplot2)
noAvizo<-Lhum_foldchange.df[1:198,]
ggplot(Lhum_foldchange.df, aes(x=Population, y=as.numeric(value), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Fold-Change (based on modern population geoMean)") +
  ggtitle("All Left Humeral Measurements") + geom_jitter(width = 0.1, na.rm = TRUE) + facet_wrap(~variable)

write.table(Lhum_foldchange, "3D_Lhum_foldchange_geoMeans.txt", sep = "\t", row.names = TRUE, col.names = TRUE)




setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##All elements, all sides (MAX)##           ##R femur MNI starts line 855##
Rfem_foldchange<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3D_Rfem_foldchange_geoMeans.txt", quote="\"")
Lfem_foldchange<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3D_Lfem_foldchange_geoMeans.txt", quote="\"")
Rhum_foldchange<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3D_Rhum_foldchange_geoMeans.txt", quote="\"")
Lhum_foldchange<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3D_Lhum_foldchange_geoMeans.txt", quote="\"")

#Rfemur<-Rfem_foldchange[,-9] #remove FHSA
#Rfemur<-Rfemur[,-6] #remove DSA
Rfemur<-Rfem_foldchange[,-2] #remove Individual
Rfemur<-Rfemur[,-1] #remove Population

#Lfemur<-Lfem_foldchange[,-9] #remove FHSA
#Lfemur<-Lfemur[,-6] #remove DSA
Lfemur<-Lfem_foldchange[,-2] #remove Individual
Lfemur<-Lfemur[,-1] #remove Population

#Rhumerus<-Rhum_foldchange[,-9] #remove HHSA
Rhumerus<-Rhum_foldchange[,-2] #remove Individual
Rhumerus<-Rhumerus[,-1] #remove Population

#Lhumerus<-Lhum_foldchange[,-9] #remove HHSA
Lhumerus<-Lhum_foldchange[,-2] #remove Individual
Lhumerus<-Lhumerus[,-1] #remove Population

Rfemur_averages<-c()
for (a in 1:nrow(Rfemur)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Rfemur[a,]), na.rm = TRUE)
  Rfemur_averages[a]<-c(z)
}
Rfemur<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
Rfemur<-data.frame(Rfemur)
names(Rfemur) = c("Population", "Individual", "RFemAverages")

Lfemur_averages<-c()
for (a in 1:nrow(Lfemur)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Lfemur[a,]), na.rm = TRUE)
  Lfemur_averages[a]<-c(z)
}
Lfemur<-cbind(Lfem_foldchange[,1], Lfem_foldchange[,2], Lfemur_averages)
Lfemur<-data.frame(Lfemur)
names(Lfemur) = c("Population", "Individual", "LFemAverages")

Rhumerus_averages<-c()
for (a in 1:nrow(Rhumerus)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Rhumerus[a,]), na.rm = TRUE)
  Rhumerus_averages[a]<-c(z)
}
hum_pop_ind<-cbind(Rhum_foldchange[,1], Rhum_foldchange[,2])
Rhumerus<-cbind(hum_pop_ind, Rhumerus_averages)
Rhumerus<-data.frame(Rhumerus)
names(Rhumerus) = c("Population", "Individual", "RHumAverages")

Lhumerus_averages<-c()
for (a in 1:nrow(Lhumerus)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(Lhumerus[a,]), na.rm = TRUE)
  Lhumerus_averages[a]<-c(z)
}
hum_pop_ind<-cbind(Lhum_foldchange[,1], Lhum_foldchange[,2])
Lhumerus<-cbind(hum_pop_ind, Lhumerus_averages)
Lhumerus<-data.frame(Lhumerus)
names(Lhumerus) = c("Population", "Individual", "LHumAverages")

femur<-merge(Rfemur, Lfemur, by="Individual")
femur<-femur[,-4] #remove extra "Mod/Sub" column
names(femur) = c("Individual", "Population", "RFemAverages", "LFemAverages")
humerus<-merge(Rhumerus, Lhumerus, by="Individual")
humerus<-humerus[,-4] #remove extra "Mod/Sub" column
sub_hum<-c("TAO-66-2", "Subfossil", 0.281129030759131, "NaN")
humerus<-rbind(humerus, sub_hum)
names(humerus) = c("Individual", "Population", "RHumAverages", "LHumAverages")

all<-merge(femur, humerus, by=c("Individual", "Population"), all=TRUE)
all_temp<-all[,-2]
all_temp<-all_temp[,-1]
averages<-c()
for (a in 1:nrow(all_temp)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(all_temp[a,]), na.rm = TRUE)
  averages[a]<-c(z)
}
averages<-cbind(all[,1], all[,2], averages)
averages<-data.frame(averages)
names(averages) = c("Individual", "Population", "Averages")

library(ggplot2)
ggplot(averages, aes(x=Population, y=as.numeric(Averages), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Aggregate Score (based on geoMean fold-changes)") +
  ggtitle("Aggregate Score of Maximum Number of Individuals") + geom_jitter(width = 0.1, na.rm = TRUE)
#ggsave(all_agg, file = "~/Desktop/Aggregate Score.pdf")

subfossil<-subset(averages, Population == "Subfossil")
modern<-subset(averages, Population == "Modern")
t.test(as.numeric(modern[,3]), as.numeric(subfossil[,3]))
#Test Name:                       Welch Two Sample t-test
#Estimated Parameter(s):          mean of x = 0.008782976
#                                 mean of y = 0.089240611
#Data:                            as.numeric(modern[, 3]) and as.numeric(subfossil[, 3])
#Test Statistic:                  t = -2.310231
#Test Statistic Parameter:        df = 12.29329
#P-value:                         0.03898436
#95% Confidence Interval:         LCL = -0.156138125
#                                 UCL = -0.004777145

average_subfossil<-mean(as.numeric(subfossil[,3]), na.rm = TRUE) #0.08924061
sd_subfossil<-sd(as.numeric(subfossil[,3]), na.rm = TRUE) #0.1173007
average_modern<-mean(as.numeric(modern[,3]), na.rm = TRUE) #0.008782976
sd_modern<-sd(as.numeric(modern[,3]), na.rm = TRUE) #0.04532682



setwd("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/")

options(stringsAsFactors = FALSE)

###stats###
##All elements, all sides (MAX)##           ##R femur MNI starts line 455##
Rfem_foldchange<-read.table("~/Dropbox/PennState/PerryLab/2014-18_Lemur/2016-2017_Beza_Taolambiby_Propithecus/Processing_Stats/3D/3D_Rfem_foldchange_geoMeans.txt", quote="\"")

#Rfemur<-Rfem_foldchange[,-9] #remove FHSA
#Rfemur<-Rfemur[,-6] #remove DSA
Rfemur<-Rfem_foldchange[,-2] #remove Individual
Rfemur<-Rfemur[,-1] #remove Population
MNI<-Rfemur[,-c(1, 2, 3, 4, 5, 6, 7)]
MNI<-MNI[,-2]#remove all measurements but FBE and DML2

Rfemur_averages<-c()
for (a in 1:nrow(MNI)){
  if ((a %% 1) == 0){print (a)}
  z<-mean(as.numeric(MNI[a,]), na.rm = TRUE)
  Rfemur_averages[a]<-c(z)
}
MNI_averages<-cbind(Rfem_foldchange[,1], Rfem_foldchange[,2], Rfemur_averages)
MNI_averages<-data.frame(MNI_averages)
names(MNI_averages) = c("Population", "Individual", "Averages")

library(ggplot2)
ggplot(MNI_averages, aes(x=Population, y=as.numeric(Averages), fill=Population)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  geom_boxplot(na.rm = TRUE) + xlab("Population") + ylab("Aggregate Score (based on geoMean fold-changes)") +
  ggtitle("MNI Subfossil Individuals (TAO-66-26, 29, 32, 33)") + geom_jitter(width = 0.1, na.rm = TRUE)

MNI_subfossil<-subset(MNI_averages, Population == "Subfossil")
MNI_modern<-subset(MNI_averages, Population == "Modern")
t.test(as.numeric(MNI_modern[,3]), as.numeric(MNI_subfossil[,3]))
#Test Name:                       Welch Two Sample t-test
#Estimated Parameter(s):          mean of x = 0.001697889
#                                 mean of y = 0.102908123
#Data:                            as.numeric(big_modern[, 3]) and as.numeric(big_subfossil[, 3])
#Test Statistic:                  t = -1.8592
#Test Statistic Parameter:        df = 3.237
#P-value:                         0.1532
#95% Confidence Interval:         LCL = -0.26749761
#                                 UCL =  0.06507714

average_MNI_subfossil<-mean(as.numeric(MNI_subfossil[,3]), na.rm = TRUE) #0.1029081
sd_MNI_subfossil<-sd(as.numeric(MNI_subfossil[,3]), na.rm = TRUE) #0.1068209
average_MNI_modern<-mean(as.numeric(MNI_modern[,3]), na.rm = TRUE) #0.00169789
sd_MNI_modern<-sd(as.numeric(MNI_modern[,3]), na.rm = TRUE) #0.05049076



header<-c("Element", "Mod%Diff", "ModStdDev", "Sub%Diff", "SubStdDev")
Rfemur<-c("Rfemur", -0.03, 3.02, 1.66, 0.58)
Lfemur<-c("Lfemur", 0.80, 1.74, 1.54, 0.66)
Rhum<-c("Rhumerus", -2.02, 7.04, 0.41, 1.67)
Lhum<-c("Lhumerus", -1.48, 6.47, NA, NA)

averages<-rbind(header, Rfemur, Lfemur, Rhum, Lhum)

t.test(as.numeric(averages[,2]), as.numeric(averages[,4], na.rm = TRUE))
#Test Name:                       Welch Two Sample t-test
#Estimated Parameter(s):          mean of x = -0.682500
#                                 mean of y = 1.203333
#Data:                            as.numeric(big_modern[, 3]) and as.numeric(big_subfossil[, 3])
#Test Statistic:                  t = -2.4778
#Test Statistic Parameter:        df = 4.6884
#P-value:                         0.05925
#95% Confidence Interval:         LCL = -3.8820274
#                                 UCL =  0.1103607



Rfemur_mod<-c(0.65, 1.59, 0.26, 0.00, 0.26, 0.53, -0.01, 0.81, -0.10, -0.25, 2.27, 2.73, -0.90, -0.82, 0.68, 2.67, 0.16, 1.24, 0.29, 0.99, -17.69, 0.29, 0.58, 1.44, 1.56)
Rfemur_sub<-c(2.02, 1.84, 1.12)
t.test(as.numeric(Rfemur_mod), as.numeric(Rfemur_sub))
#Test Name:                       Welch Two Sample t-test
#Estimated Parameter(s):          mean of x = -0.0308
#                                 mean of y = 1.6600 
#Data:                            as.numeric(big_modern[, 3]) and as.numeric(big_subfossil[, 3])
#Test Statistic:                  t = -2.0929
#Test Statistic Parameter:        df = 25.457
#P-value:                         0.04648
#95% Confidence Interval:         LCL = -3.3531467
#                                 UCL =  -0.0284533

Lfemur_mod<-c(1.88, 0.98, 0.25, -0.21, -0.17, 2.08, 1.66, 0.95, 1.94, 0.26, 0.36, 1.45, 0.61, 1.23, -1.07, 1.42, 0.70, 0.60, 0.55, 0.78, 1.55, -0.30, 0.48, 1.21)
Lfemur_sub<-c(-0.48, 1.61, 3.50)
t.test(as.numeric(Lfemur_mod), as.numeric(Lfemur_sub))
#Welch Two Sample t-test
#t = -0.64081, df = 2.0793, p-value = 0.5851
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -5.559484  4.071984
#sample estimates:
#  mean of x mean of y 
#0.7995833 1.5433333 


Rhum_mod<-c(0.30, -3.33, -1.28, 2.17, -3.09, -2.38, -0.99, -1.02, -2.87, 2.82, 2.88, -5.49, -1.09, -1.70, -2.81, -2.35, -6.85, -5.64, -1.17, -0.81, -7.32, 0.10, -4.49)
Rhum_sub<-c(0.41)
t.test(as.numeric(Rhum_mod), as.numeric(Rhum_sub))
#Error in t.test.default(as.numeric(Rhum_mod), as.numeric(Rhum_sub)) : not enough 'y' observations
