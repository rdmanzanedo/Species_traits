###load the data from Sci.data sdms models
##set wd
#install.packages('raster', dependencies=T)
library(raster)

##load the rasters for the multiple species
#we are going to use the ensemble sdm, using current probability and potential distribution
abal = raster("Abies_alba_ens-sdms_cur_prob_pot.tif")
acca = raster("Acer_campestre_ens-sdms_cur_prob_pot.tif")
acpl = raster("Acer_platanoides_ens-sdms_cur_prob_pot.tif")
acps = raster("Acer_pseudoplatanus_ens-sdms_cur_prob_pot.tif")
algl = raster("Alnus_glutinosa_ens-sdms_cur_prob_pot.tif")
alin = raster("Alnus_incana_ens-sdms_cur_prob_pot.tif")
bepe = raster("Betula_pendula_ens-sdms_cur_prob_pot.tif")
cabe = raster("Carpinus_betulus_ens-sdms_cur_prob_pot.tif")
fasy = raster("Fagus_sylvatica_ens-sdms_cur_prob_pot.tif")
lade = raster("Larix_decidua_ens-sdms_cur_prob_pot.tif")
piab = raster("Picea_abies_ens-sdms_cur_prob_pot.tif")
pice = raster("Pinus_cembra_ens-sdms_cur_prob_pot.tif")
#pimu = raster("Pinus_mugo_ens-sdms_cur_prob_pot.tif")
pini = raster("Pinus_nigra_ens-sdms_cur_prob_pot.tif")
pisy = raster("Pinus_sylvestris_ens-sdms_cur_prob_pot.tif")
prav = raster("Prunus_avium_ens-sdms_cur_prob_pot.tif")
qepe = raster("Quercus_petraea_ens-sdms_cur_prob_pot.tif")
quro = raster("Quercus_robur_ens-sdms_cur_prob_pot.tif")
soau = raster("Sorbus_aucuparia_ens-sdms_cur_prob_pot.tif")
#soto = raster("Sorbus_torminalis_ens-sdms_cur_prob_pot.tif")
tico = raster("Tilia_cordata_ens-sdms_cur_prob_pot.tif")
ulgl = raster("Ulmus_glabra_ens-sdms_cur_prob_pot.tif")

#two species are not available
par(mfcol=c(4,5))
species = c(abal, acca, acpl, acps, algl,alin,
            bepe, cabe, fasy, lade, piab,pice,
            pini, pisy, prav, qepe, quro, soau,
            tico, ulgl)
for (i in 1:20){
  plot(species[[i]], main=substr(paste(names(species[[i]])),1,12), line=-0.5)
}

#get the climate data
bio_curr_var <- getData('worldclim', download=TRUE,var='bio', res=2.5)

max.temp.warmest = bio_curr_var$bio5
min.temp.coldest = bio_curr_var$bio6
prec.driest.quar = bio_curr_var$bio17

#install.packages('virtualspecies')
library(virtualspecies)

#dataframe for results
results= data.frame('species' = rep(NA,20), 'sampling.size'=rep(NA,20),
                    'tmax90' = rep(NA,20),'tmax95' = rep(NA,20),'tmax99' = rep(NA,20),
                    'tmin10' = rep(NA,20),'tmin05' = rep(NA,20),'tmin01' = rep(NA,20),
                    'dro10' = rep(NA,20),'dro05' = rep(NA,20),'dro01' = rep(NA,20))
par(mfcol=c(1,1))
##loopyloop 
for (i in c(1:20)){
  #set the sampling random term
  set.seed(23)
  
  #make sure the raster is right
  plot(species[[i]], main=substr(paste(names(species[[i]])),1,12), line=-1, axes=F)
  
  ##the sampling size should be dependent on the size, let's use 
  sampling.n = length(species[[i]][!is.na(species[[i]])])/6
  
  #sampling 9999 occurences, wiht probability of sampling directly linked to the suitabitliy (allows us to not care about thresholds)
  points_real= sampleOccurrences(species[[i]]>0, sampling.n, bias='manual', weights=species[[i]])
  
  #put sampled points together
  points_s = data.frame(points_real$sample.points$x, points_real$sample.points$y)
  
  #sample the 3 interest variables
  sam_max.temp.warmest <- raster::extract(max.temp.warmest,points_s)
  sam_min.temp.coldest <- raster::extract(min.temp.coldest,points_s)
  sam_prec.driest.quar <- raster::extract(prec.driest.quar,points_s)
  
  #distributions? are they roughly normal looking
  par(mfcol=c(1,3))
  plot(density(na.omit(sam_max.temp.warmest/10)))
  plot(density(na.omit(sam_min.temp.coldest/10)))
  plot(density(na.omit(sam_prec.driest.quar)))
  par(mfcol=c(1,1))
  
  #quantiles 90,95,99 for them
  max.t = quantile(na.omit(sam_max.temp.warmest)/10, probs=c(0.9,0.95,0.99))
  min.t = quantile(na.omit(sam_min.temp.coldest)/10, probs=c(0.01,0.05,0.1))
  dro.t = quantile(na.omit(sam_prec.driest.quar), probs=c(0.01,0.05,0.1))  
  
  #save the results before continue looping
  results$species[i] = substr(paste(names(species[[i]])),1,12)
  results$sampling.size[i] = sampling.n
  results$tmax90[i] = max.t[1]
  results$tmax95[i] = max.t[2]
  results$tmax99[i] = max.t[3]
  results$tmin10[i] = min.t[3]
  results$tmin05[i] = min.t[2]
  results$tmin01[i] = min.t[1]
  results$dro10[i] = dro.t[3]
  results$dro05[i] = dro.t[2]
  results$dro01[i] = dro.t[1]
}

#we can now do it for number 12 with lower sample size ()
head(results)
results$sp.code = c('abal','acca','acpl','acps','algl','alin','bepe','cabe',
                    'fasy','lade','piab','pice','pini','pisy','prav',
                    'qupe','quro','soau','tico','ulgl')

#load data of kristiina
setwd("D:/Universidad/2022/Traits Kristiina project")
data.kv = read.csv('traits table KV.csv', header=T)
head(data.kv)
head(results)

#merge both
all.data = merge(data.kv, results, by='sp.code')
names(all.data)
all.data2 = all.data[,c(1,3,4,6,8,10,12,14,17:25)]
head(all.data2)
plot(all.data2)

#install.packages('corrplot')
library(corrplot)

cors=cor(all.data2[,c(4,5,7,10,13,15)],use='complete.obs') # evaluate correlations
corrplot(cors, addCoef.col = "grey",number.cex=.6)
plot(prec.driest.quar$bio17)

cors=cor(results[,c(3,6,9)],use='complete.obs') # evaluate correlations
corrplot(cors, addCoef.col = "grey",number.cex=.6)

write.csv()