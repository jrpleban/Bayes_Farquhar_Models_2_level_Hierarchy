## Bayesian Estimate of PS paramters from Brassica ACi curves ###
setwd("~/Documents/ACiP4/Scripts")   ### Change 4 server
require(rjags)### Change 4 server

source("CiCc_Temp_Jf_model_H.R") ### Bayesian Modelw/ Temp Dependency & farQ_tempDep limitation


#### Set up for rjags #########
parameters = c("Vcmax25", "Rd25","gammaS25","gm25","Kc25", "Ko25",
"mu.Vcmax25", "mu.Rd25","mu.gammaS25","mu.gm25",
"tau.Vcmax25", "tau.Rd25","tau.gammaS25","tau.gm25",
"EgammaS" ,
"ERd" ,"Egm","EKc", "EKo",  "EVcmax",
"tau")

adaptSteps = 100000             # Number of steps to "tune" the samplers.
burnInSteps = 200000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps=50000        # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per

######################
######################
#####DATA SET UP######
######################
ACdat<-read.delim("~/Documents/ACiP4/Data/ACdata_reduced1_clean.txt")  ### Change 4 server
#useID<-read.delim("~/Documents/ACiP4/Data/IDs_to_use-9-9-17.txt")
#gID<-useID$ID
#IDgeno<-as.character(useID$geno)

ACi<-ACdat#[ACdat$AC.ID %in% gID, ]

####  DATA SET  UP  ############
#### PULL OUT Data for each genotype   ###
ID<-sort(unique(ACi$geno))
#### PULL OUT Data for each genotype   ###
ID<-sort(unique(ACi$geno))
#### ID Data
IDdat<-vector("list", 6)
for(i in 1:6){
    IDdat[[i]]<-as.factor(ACi[ACi$geno==ID[i],]$AC.ID)
}
#### Photo Data
A = vector("list", 6)
for(i in 1:6){
    A[[i]]<-ACi[ACi$geno==ID[i],]$Photo
}
## CP = CO2 concentraion intercelluar space (in partial pressure)
ACi$CP<-ACi$Ci*ACi$Press/1000   ## calcs CO2 concentraion intercelluar space (in partial pressure)
CP = vector("list", 6)   ### pulls for each ID
for(i in 1:6){
    CP[[i]]<-ACi[ACi$geno==ID[i],]$CP
}

### Leaf Temp (C) ###
T = vector("list", 6) ### pulls for each ID
for(i in 1:6){
    T[[i]]<-ACi[ACi$geno==ID[i],]$Tleaf
}


### Leaf Temp (K) ###
K = vector("list", 6)  ### pulls for each ID
for(i in 1:6){
    K[[i]]<-ACi[ACi$geno==ID[i],]$Tleaf+273.15
}


###  O2 in pp from atmospheric pressure (kPa)
### PP(Pa) =  MolFact * pressure of the gas mixture
ACi$O<-(.21)*(ACi$Press*1000)
O = vector("list", 6)   ### pulls for each ID
for(i in 1:6){
    O[[i]]<-ACi[ACi$geno==ID[i],]$O
}

### PARi
Q = vector("list", 6)   ### pulls for each ID
for(i in 1:6){
    Q[[i]]<-ACi[ACi$geno==ID[i],]$PARi
}

### ETR
Jf= vector("list", 6)   ### pulls for each ID
for(i in 1:6){
    Jf[[i]]<-ACi[ACi$geno==ID[i],]$ETR
}


#### Sample size
#N<-36
## Contants ##
Kref=298.15; R=0.008314; Tref=25

datalist1<-list(N=length(A[[1]]), NID=length(unique(IDdat[[1]])), An=A[[1]], CiP=CP[[1]], T=T[[1]], K=K[[1]], Jf=Jf[[1]], O=O[[1]], ID=IDdat[[1]], Kref=298.15, R=0.008314, Tref=25) #301
datalist2<-list(N=length(A[[2]]), NID=length(unique(IDdat[[2]])), An=A[[2]], CiP=CP[[2]], T=T[[2]], K=K[[2]], Jf=Jf[[2]], O=O[[2]], ID=IDdat[[2]], Kref=298.15, R=0.008314, Tref=25) #46
datalist3<-list(N=length(A[[3]]), NID=length(unique(IDdat[[3]])), An=A[[3]], CiP=CP[[3]], T=T[[3]], K=K[[3]], Jf=Jf[[3]], O=O[[3]], ID=IDdat[[3]], Kref=298.15, R=0.008314, Tref=25) #bro
datalist4<-list(N=length(A[[4]]), NID=length(unique(IDdat[[4]])), An=A[[4]], CiP=CP[[4]], T=T[[4]], K=K[[4]], Jf=Jf[[4]], O=O[[4]], ID=IDdat[[4]], Kref=298.15, R=0.008314, Tref=25)  #cab
datalist5<-list(N=length(A[[5]]), NID=length(unique(IDdat[[5]])), An=A[[5]], CiP=CP[[5]], T=T[[5]], K=K[[5]], Jf=Jf[[5]], O=O[[5]], ID=IDdat[[5]], Kref=298.15, R=0.008314, Tref=25)  #oil
datalist6<-list(N=length(A[[6]]), NID=length(unique(IDdat[[6]])), An=A[[6]], CiP=CP[[6]], T=T[[6]], K=K[[6]], Jf=Jf[[6]], O=O[[6]], ID=IDdat[[6]], Kref=298.15, R=0.008314, Tref=25) #tur
#################################
##################################
##### Impliment model in JAGS ####
#################################
### running each curve
source("CiCc_Temp_Jf_model_H.R") ### Bayesian Modelw/ Temp Dependency & farQ_tempDep limitation
print("initialize models")

model1 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist1, n.chains=nChains , n.adapt=adaptSteps)

model2 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist2, n.chains=nChains , n.adapt=adaptSteps)


model3 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist3, n.chains=nChains , n.adapt=adaptSteps)
model4 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist4, n.chains=nChains , n.adapt=adaptSteps)
model5 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist5, n.chains=nChains , n.adapt=adaptSteps)
model6 <- jags.model(textConnection(CiCc_Temp_Jf),
data = datalist6, n.chains=nChains , n.adapt=adaptSteps)

#################################
print("updating")
update(model1, burnInSteps)
update(model2, burnInSteps)
update(model3, burnInSteps)
update(model4, burnInSteps)
update(model5, burnInSteps)
update(model6, burnInSteps)

##########################################
##### mcmc_samples  model 6 genos & IDs  #####
#####   ID order   1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
##########################################
###### SAMPLE all 18 individual curves ###
##########################################
##########################################
print("sampling chains")
##########################################
mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
####### Plot results #####
#plot(mcmc_samples1)
mcmcChain_r301_5 = as.matrix( mcmc_samples1)
chainLength = NROW(mcmcChain_r301_5)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain_r301_5[, "tau" ] )
mcmcChain_r301_5 = as.data.frame(cbind( mcmcChain_r301_5, sigma ))
g1<-gelman.diag(mcmc_samples1)
####
mcmc_samples2<- coda.samples(model2,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_r46_5 = as.matrix( mcmc_samples2)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain_r46_5[, "tau" ] )
mcmcChain_r46_5 = as.data.frame(cbind( mcmcChain_r46_5, sigma ))
g2<-gelman.diag(mcmc_samples2)
######
mcmc_samples3 <- coda.samples(model3,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_bro_5 = as.matrix( mcmc_samples3)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_bro_5[, "tau" ] )
mcmcChain_bro_5 = as.data.frame(cbind( mcmcChain_bro_5, sigma ))
g3<-gelman.diag(mcmc_samples3)
#####
mcmc_samples4 <- coda.samples(model4,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_cab_5 = as.matrix( mcmc_samples4)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_cab_5[, "tau" ] )
mcmcChain_cab_5 = as.data.frame(cbind( mcmcChain_cab_5, sigma ))
g4<-gelman.diag(mcmc_samples4)
######
mcmc_samples5 <- coda.samples(model5,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_oil_5 = as.matrix( mcmc_samples5)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_oil_5[, "tau" ] )
mcmcChain_oil_5 = as.data.frame(cbind( mcmcChain_oil_5, sigma ))
g5<-gelman.diag(mcmc_samples5)
######
mcmc_samples6 <- coda.samples(model6,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_tur_5 = as.matrix( mcmc_samples6)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_tur_5[, "tau" ] )
mcmcChain_tur_5 = as.data.frame(cbind( mcmcChain_tur_5, sigma ))
g6<-gelman.diag(mcmc_samples6)

##########################################

##########################################

######### Gelman potential scale reduction factor summary for all par and multivariate #############
gelmin<-c(min(g1$psrf[,1]), min(g2$psrf[,1]),min(g3$psrf[,1]), min(g4$psrf[,1]),
min(g5$psrf[,1]), min(g6$psrf[,1]))


gelmax<-c(max(g1$psrf[,1]), max(g2$psrf[,1]),max(g3$psrf[,1]), max(g4$psrf[,1]),
max(g5$psrf[,1]), max(g6$psrf[,1]))

gelmulit<-c(g1$mpsrf, g2$mpsrf, g3$mpsrf, g4$mpsrf,g5$mpsrf, g6$mpsrf)

gelsum<-as.data.frame(cbind(gelmin, gelmax,gelmulit))
colnames(gelsum)<-c("gel_min", "gel_max","gel_multi")

##########################################
#########  SAVE Samples ####################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
###
print("writing samples")


setwd("~/Documents/ACiP4/Post_Data_H")
write.table(mcmcChain_r301_5,"mcmcChain_r301_5", sep="\t", col.name=TRUE)
write.table(mcmcChain_r46_5,"mcmcChain_r46_5", sep="\t", col.name=TRUE)
write.table(mcmcChain_bro_5,"mcmcChain_bro_5", sep="\t", col.name=TRUE)
write.table(mcmcChain_cab_5,"mcmcChain_cab_5", sep="\t", col.name=TRUE)
write.table(mcmcChain_oil_5,"mcmcChain_oil_5", sep="\t", col.name=TRUE)
write.table(mcmcChain_tur_5,"mcmcChain_tur_5", sep="\t", col.name=TRUE)
#########################################################################

#################    PARAMETER ESTIMATES      ###########################

#########################################################################

#########################################################################
###  for revision 3 (2_24_18) unbalanced data results in a different number of ind par est ###
###  for simplicity select only 3 indivudals for saving in data structure
common_cols <- colnames(mcmcChain_oil_5)
mcmcChain_r301_5<-subset(mcmcChain_r301_5, select = common_cols)
mcmcChain_r46_5<-subset(mcmcChain_r46_5, select = common_cols)
mcmcChain_bro_5<-subset(mcmcChain_bro_5, select = common_cols)
mcmcChain_cab_5<-subset(mcmcChain_cab_5, select = common_cols)
mcmcChain_tur_5<-subset(mcmcChain_tur_5, select = common_cols)
#############################################


#############################################
EstPars1<-apply(mcmcChain_r301_5, 2, median)
EstPars2<-apply(mcmcChain_r46_5, 2, median)
EstPars3<-apply(mcmcChain_bro_5, 2, median)
EstPars4<-apply(mcmcChain_cab_5, 2, median)
EstPars5<-apply(mcmcChain_oil_5, 2, median)
EstPars6<-apply(mcmcChain_tur_5, 2, median)

ParEst<-rbind(EstPars1,EstPars2,EstPars3,EstPars4,EstPars5,EstPars6)


#########################################################################

#########################################################################

#########################################################################

#########################################################################

#################   Calc DIC with pD criteria  ###########################

#########################################################################

#########################################################################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760

print("sampling for DIC")
dic1 <- dic.samples(model1, DICsteps, "pD")
dic1a <- dic.samples(model1, DICsteps, "popt")
dic_r301_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model2, DICsteps, "pD")
dic1a <- dic.samples(model2, DICsteps, "popt")
dic_r46_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model3, DICsteps, "pD")
dic1a <- dic.samples(model3, DICsteps, "popt")
dic_bro_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model4, DICsteps, "pD")
dic1a <- dic.samples(model4, DICsteps, "popt")
dic_cab_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model5, DICsteps, "pD")
dic1a <- dic.samples(model5, DICsteps, "popt")
dic_oil_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model6, DICsteps, "pD")
dic1a <- dic.samples(model6, DICsteps, "popt")
dic_tur_5<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))


dics<-rbind(dic_r301_5, dic_r46_5, dic_bro_5, dic_cab_5, dic_oil_5,dic_tur_5)
colnames(dics)<-c("pD_dev", "pD_pen", "pD_DIC", "popt_dev", "popt_pen", "popt_DIC")


modelname<-rep("CiCc_Temp_Jf",6)

###################################################

###################################################

#   FINAL TABLE W/ median estimates and DICS ########

###################################################
print("writing PAR Estimates")
ParEstFull<-as.data.frame(cbind(modelname, ID, ParEst,dics, gelsum))
write.table(ParEstFull,"ParEstFull_model_5", sep="\t", col.name=TRUE,row.names = F)

print("Finished")



require(ggmcmc)

pdf("~/Documents/ACiP4/Trace_plots/density_r301_5.pdf",
colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_46_5.pdf",
colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_bro_5.pdf", colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_cab_5.pdf", colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_oil_5.pdf", colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_tur_5.pdf", colormodel='cmyk',width=6.25, height=(180))
ggs_density(ggs(mcmc_samples6))
dev.off()




pdf("~/Documents/ACiP4/Trace_plots/cor_r301_5.pdf",
colormodel='cmyk',width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_46_5.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_bro_5.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_cab_5.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_oil_5.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_tur_5.pdf", colormodel='cmyk',width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples6))
dev.off()


pdf("~/Documents/ACiP4/Trace_plots/Rhat_r301_5.pdf",
colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_46_5.pdf",
colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_bro_5.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_cab_5.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_oil_5.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_tur_5.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples6))
dev.off()



