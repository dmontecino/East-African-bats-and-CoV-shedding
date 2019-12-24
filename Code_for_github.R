# ------------------------------------------------------------------------------------ #
# Script showing the data analysis: descriptive as well as statistical modeling of the #
# project "Reproduction of East-African bats may guide risk mitigation for coronavirus #
# spillover" published in the Journal One health Outlook                               #  
# Author: Diego Montecino-Latorre                                                      #  
# ------------------------------------------------------------------------------------ #

#---load packages needed ---#

library('R2jags')
library("rethinking")
library("sp")
library("lubridate")
library("Hmisc")
library("GGally")
library("rstan")
library("sp")
library("plyr")
library("rstanarm")
library("loo")
library("stringr")

# ---download data (constructed in Bats/Chapter 2/Data analysis/code to construct dataset) ---#

download.file(url="https://ndownloader.figshare.com/files/20411751",
              destfile = "/bats.csv")

# ---load data---#

bats=read.csv("/bats.csv")

# --- some descriptives ---#

n.obs.bats=nrow(bats) # 747 wo Epomophorus
n.sites.bats=length(unique(bats$SiteName)) # 34 sites

#sites
bats$SiteName=factor(bats$SiteName)
sites.bats=as.numeric(bats$SiteName)


#Coronavirus (alpaha or beta)
CoV=ifelse(c(bats$Result.alpha + bats$Result.beta)>0,1,0)

# Age
age.bats=ifelse(bats$AgeClass=="Adult",0,1)


table(bats$Age)
table(bats$new.season)
table(bats$Age, bats$new.season)
table(bats$Family)
table(bats$Family,bats$new.season)
table(bats$Family,bats$Age)

tapply(CoV, list(bats$Age), mean)
tapply(CoV, list(bats$new.season), mean)
tapply(CoV, list(bats$Age, bats$new.season), mean)
tapply(CoV, list(bats$Family, bats$new.season), mean)
tapply(CoV, list(bats$Family, bats$Age), mean)


#################
## STAN MODELS ##
#################

sort(unique(bats$SiteName))

#------#
# DATA #
#------#

#-----------------------------#
# data to model CoV detection #
#-----------------------------#

#------ number of observations ------#
N=nrow(bats) 

#------ Results of CoV in samples .Binary outcome of interest ------#
CoV=ifelse(c(bats$Result.alpha + bats$Result.beta)>0,1,0)

#------ ages ------#
age=ifelse(bats$Age=="Juvenile",1,0)


#------ sampling event vector ------#

sort(unique(bats$SiteName))
sampling=factor(bats$SiteName)
sampling=as.numeric(sampling)

#------ number of sampling events ------#
n.sampling=length(unique(sampling))

#------ Families ------#

# All hipossideros are hipposideros
bats$Family=factor(bats$Family) #family
family=as.integer(bats$Family) #family
n.family=length(unique(family)) #number of unique families in the dataset

# ------------------------------------------------------#




# ------------------------------------------------------#
# data imputation of the repro season when non-inferred #
# ------------------------------------------------------#

lat=c(scale(bats$Latitude, center = T)) #latitude

# --- historical precipitation --- #
# Move the historical precipitation to a category pf dry or wet season
# category added according to pripitation at each site

bats$historic.monthly.prep.cat=ifelse(bats$historic.monthly.prep>100,"wet","dry")

table(bats$historic.monthly.prep.cat, bats$Age)
table(bats$historic.monthly.prep.cat)
table(bats$historic.monthly.prep.cat, bats$new.season)
historical_prep=bats$historic.monthly.prep.cat # historical precipitation

# --- Number of litters per year per species --- #

bats$litters.per.year=NA # litters per year
bats$litters.per.year[bats$barcode.sp%in%c("Neoromicia nanus", "Hipposideros ruber", "Hipposideros cf. caffer", # one
                                           "Nycteris cf. thebaica", "Rhinolophus cf. clivosus",
                                           "Pipistrellus cf. hesperidus", "Triaenops persicus", 
                                           "Eidolon helvum")]=1
bats$litters.per.year[bats$barcode.sp%in%c("Coleura afra", "Mops condylurus", "Taphozous mauritianus", "Rousettus aegyptiacus", "Lissonycteris angolensis")]=2 # two
bats$litters.per.year[bats$barcode.sp=="Chaerephon pumilus"]=3 #three litters per year
unique(bats[is.na(bats$litters.per.year),]$barcode.sp)

litters_per_year=ifelse(bats$litters.per.year>1,1,0)

# --- Day of the year the bats were sampled --- #

dates <- bats$EventDate # Day of the year bats were sampled #
dates <- ymd(dates)
bats$Day.of.the.year=yday(dates)
day_of_the_year=c(scale(bats$Day.of.the.year, center = T))

# --- Reproductive season missing or not --- #

season_cat0=bats$new.season
season_cat=NA
season_cat[season_cat0=="Recently weaned"]=1
season_cat[season_cat0=="Not recently weaned"]=0
table(season_cat)
season_cat_miss=ifelse(is.na(season_cat), 1, 0)
season_cat[is.na(season_cat)]=(-1)

# --- Bats touch in the roost --- #

touch=NA
#touch=ifelse(bats$touch.roost=="Yes",1,0)

touch[bats$barcode.sp=="Chaerephon pumilus"]=0
touch[bats$barcode.sp=="Coleura afra"]=0
touch[bats$barcode.sp=="Hipposideros cf. caffer"]=1
touch[bats$barcode.sp=="Hipposideros ruber"]=1
touch[bats$barcode.sp=="Mops condylurus"]=1
touch[bats$barcode.sp=="Neoromicia nanus"]=0
touch[bats$barcode.sp=="Nycteris cf. thebaica"]=0
touch[bats$barcode.sp=="Pipistrellus cf. hesperidus"]=1
touch[bats$barcode.sp=="Rhinolophus cf. clivosus"]=0
touch[bats$barcode.sp=="Taphozous mauritianus"]=0
touch[bats$barcode.sp=="Triaenops persicus"]=1
touch[bats$barcode.sp=="Eidolon helvum"]=1
touch[bats$barcode.sp=="Lissonycteris angolensis"]=0
touch[bats$barcode.sp=="Rousettus aegyptiacus"]=1

tapply(CoV, list(touch), mean)


## --- FINAL MODEL HAS AGE THE REPRODUCTIVE PERIOD, THE GROUPING 
## --- PER SAMPLING EVENT AND SPECIES AND A 
## --- A COVARIATE FOR E. HELVUM AND T. PERSICUS

# ----------------#
# STAN MODEL CODE #
# ----------------#

modelstring <- "
data{

int N; //number of observations
int L; //number of unique species
int J; //number of unique sampling events


int CoV[N]; // binary outcome


//repro season covariate
int season_cat[N]; // categorical variable when observed
int season_cat_miss[N]; // categorical variable when unobserved (index with 1's and zero's). 1 if unobserved

//age covariate
int age[N];

//age covariate
int eid[N];

// index of sampling events
int index_sampling[N];

// index of species
int index_species[N];

//data to impute repoductive season
real day_of_the_year[N];
real lat[N];
int historical_prep[N];
int litters_per_year[N];
}

parameters{
real  alpha;
real  beta_rec_weaned; // coefficient of the binary outcome as a function of the categorical variable level 2 (dummy 2)
real  beta_age; // coefficient of the binary outcome as a function of the categorical variable level 2 (dummy 2)
real  beta_eid; // coefficient of the binary outcome as a function of the categorical variable level 2 (dummy 2)

real a_imp; // intercept for the imputation model and in the chunk to model the probability of category 2
real b1_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real b2_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real b3_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real b4_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2


real<lower=0> sigma_theta_sampling; // standard deviation across families - specific intercepts
real<lower=0> sigma_theta_species; // standard deviation across families- specific intercepts
vector[L] epsilon;
vector[L] mu;

vector[J] kappa;
vector[J] mu_2;

}

transformed parameters {
vector[L] theta_species;
vector[J] theta_site;

theta_species = mu + sigma_theta_species * epsilon;
theta_site = mu_2 + sigma_theta_sampling * kappa;


}


model{
// priors
alpha ~ normal(0,1.5); // explained above
beta_rec_weaned ~ normal(0,1.5); // explained above
beta_age ~ normal(0,1.5); // explained above
beta_eid ~ normal(0,1.5); // explained above
epsilon~normal(0,1.5);
kappa~normal(0,1.5);

sigma_theta_species~cauchy(0,5); //standard deviation across family- specific intercepts
sigma_theta_sampling~cauchy(0,5); //standard deviation across sampling events - specific intercepts

a_imp ~ normal(0,1.5); // explained above
b1_imp ~ normal(0,1.5); // explained above
b2_imp ~ normal(0,1.5); // explained above
b3_imp ~ normal(0,1.5); // explained above
b4_imp ~ normal(0,1.5); // explained above

for(l in 1:L){
mu[l]~normal(0,1.5);
}

for(j in 1:J){
mu_2[j]~normal(0,1.5);
}

//add Data in the model
for (i in 1:N) {

vector[2] p;
p[1] = a_imp + b1_imp*day_of_the_year[i] + b2_imp*lat[i] + b3_imp*litters_per_year[i] + b4_imp*historical_prep[i]; // modeling the prob 2 as a function of the covariate to model the category
p[2]=1-p[1];


if (season_cat_miss[i] == 0) {
CoV[i] ~  bernoulli_logit(alpha+
                          beta_rec_weaned*season_cat[i]+
                          beta_age*age[i]+
                          beta_eid*eid[i]+
                          theta_site[index_sampling[i]]+
                          theta_species[index_species[i]]);

season_cat[i] ~ bernoulli(softmax(p)[1]);

}

// x missing model the category and posteiorly model the binary outcome as a function of the imputed category
else {

vector[2] log_lik_cats; // vector to hold the log probabilities for each alternate scenario (category 1, 2, or 3)

log_lik_cats[1] = softmax(p)[1] + bernoulli_logit_lpmf( CoV[i] | alpha + beta_rec_weaned + beta_age*age[i] + beta_eid*eid[i] + theta_site[index_sampling[i]] + theta_species[index_species[i]]);  // category 1: recently weaned
log_lik_cats[2] = softmax(p)[1] + bernoulli_logit_lpmf( CoV[i] | alpha + beta_age*age[i] + beta_eid*eid[i] + theta_site[index_sampling[i]] + theta_species[index_species[i]]);  // category 0: not recently weaned


target += log_sum_exp(log_lik_cats); // sum log probabilities across the scenarios (i.e., marginalize over missingness)
}
}
} // close model block


generated quantities { // generate estimates of the imputed category for all observations

vector[N] x_imp_all;
vector[N] x_imp_unobs;
int y_rep[N];

for (i in 1:N) {

vector[2] p;
p[1] = a_imp + b1_imp*day_of_the_year[i] + b2_imp*lat[i] + b3_imp*litters_per_year[i] + b4_imp*historical_prep[i]; // modeling the prob 2 as a function of the covariate to model the category
p[2]=1-p[1];

if (season_cat_miss[i] == 1) {

x_imp_all[i] = bernoulli_rng(softmax(p)[1]); 
x_imp_unobs[i] = bernoulli_rng(softmax(p)[1]); // realization of the actual category
} 

else {

x_imp_all[i] = bernoulli_rng(softmax(p)[1]); // realization of the actual category
x_imp_unobs[i] = season_cat[i]; // when the category has been observed, we know the category
}


y_rep[i] = bernoulli_logit_rng(alpha + beta_rec_weaned * x_imp_unobs[i] + beta_age*age[i] + beta_eid*eid[i] + theta_species[index_species[i]] + theta_site[index_sampling[i]]);

} // close loop
} // close generated quantities  block

"


# ---- Data for the model ---- #

dat <- list(N = N,
            L = n.species,
            J = n.sampling,
            index_species=species,
            index_sampling=sampling,
            CoV=CoV,
            season_cat=season_cat,
            season_cat_miss=season_cat_miss,
            eid=ifelse(bats$barcode.sp%in%c("Eidolon helvum", "Triaenops persicus"),1,0),
            age=age,
            day_of_the_year=day_of_the_year,
            lat=lat,
            historical_prep=historic.monthly.prep.cat,
            litters_per_year=litters_per_year)


# ---- Run the model ---- #

stan.model.season.age.species.sites.eid=stan(model_code = modelstring,  
                                             iter = 4000, warmup = 3000, #thin=1,
                                             chains = , cores = 1,
                                             data = dat,
                                             control=list(adapt_delta=0.995, max_treedepth = 15))


# alpha

HPDI(c(extract(stan.model.season.age.species.sites.eid)$alpha), 0.9)


# coeff recent weaned

HPDI(c(extract(stan.model.season.age.species.sites.eid)$beta_rec_weaned), 0.9)

# coeff age

HPDI(c(extract(stan.model.season.age.species.sites.eid)$beta_age), 0.9)

# coeff eid

HPDI(c(extract(stan.model.season.age.species.sites.eid)$beta_eid), 0.9)


#get the samples per parameter#
for(i in c(1:11, 16,17)){
  if(class(extract(stan.model.season.age.species.sites.eid)[[i]])=="matrix"){
    for(j in 1:ncol(extract(stan.model.season.age.species.sites.eid)[[i]])){
      assign(paste0(names(extract(stan.model.season.age.species.sites.eid))[i],".",j), value = extract(stan.model.season.age.species.sites.eid)[[i]][,j])}}
  else{assign(names(extract(stan.model.season.age.species.sites.eid))[i], value = extract(stan.model.season.age.species.sites.eid)[[i]])}}


# assess correlation among samples 

assess.cor=cor(cbind(
  sapply(names(extract(stan.model.season.age.species.sites.eid)[c(1:11)]), function(x) get(x)),
  do.call(cbind, mget(ls()[grep(pattern = "theta_site.", ls(), fixed = T)])),
  do.call(cbind, mget(ls()[grep(pattern = "theta_species.", ls(), fixed = T)]))))

print(assess.cor, max = 3000)

#OR when recently weaned

exp(HPDI(c(extract(stan.model.season.age.species.sites.eid)$beta_rec_weaned), 0.9))

#OR juvenile

exp(HPDI(c(extract(stan.model.season.age.species.sites.eid)$beta_age), 0.9))


# Detection from the posterior prob. distributions

juvenile.not.weaned=extract(stan.model.season.age.species.sites.eid)$alpha+
  extract(stan.model.season.age.species.sites.eid)$beta_age

inv_logit(HPDI(c(juvenile.not.weaned), 0.9))


juvenile.weaned=extract(stan.model.season.age.species.sites.eid)$alpha+
  extract(stan.model.season.age.species.sites.eid)$beta_age+
  extract(stan.model.season.age.species.sites.eid)$beta_rec

inv_logit(HPDI(c(juvenile.weaned), 0.9))



adult.weaned=c(extract(stan.model.season.age.species.sites.eid)$alpha+
                 extract(stan.model.season.age.species.sites.eid)$beta_rec)

inv_logit(HPDI(c(adult.weaned), 0.9))



adult.not.weaned=c(extract(stan.model.season.age.species.sites.eid)$alpha)


inv_logit(HPDI(c(adult.not.weaned), 0.9))



# Detection based on posterior predicitve distributions and plot. 


# Detections 
# detection adults when not recently weaned

adult.not.weaned=c(extract(stan.model.season.age.species.sites.eid)$alpha)
adult.not.weaned=exp(adult.not.weaned)/(1+exp(adult.not.weaned))

# simulating realizations of detection in 1000 individuals with each of the 1000 prob values
adult.not.weaned=sapply(adult.not.weaned, function(x) rbinom(n = 1000, prob = x, size = 1))

# detection per simulation
adult.not.weaned=apply(adult.not.weaned, MARGIN = 2, mean)
#dens(adult.not.weaned, xlim=c(0,0.1))



# detection juveniles when not recently weaned

juvenile.not.weaned=extract(stan.model.season.age.species.sites.eid)$alpha+
  extract(stan.model.season.age.species.sites.eid)$beta_age

juvenile.not.weaned=exp(juvenile.not.weaned)/(1+exp(juvenile.not.weaned))

# simulating realizations of detection in 1000 individuals with each of the 1000 prob values
juvenile.not.weaned=sapply(juvenile.not.weaned, function(x) rbinom(n = 1000, prob = x, size = 1))

# detection per simulation
juvenile.not.weaned=apply(juvenile.not.weaned, MARGIN = 2, mean)



# detection adults recently weaned

adult.weaned=sample(c(extract(stan.model.season.age.species.sites.eid)$alpha+
                        extract(stan.model.season.age.species.sites.eid)$beta_rec), 1000)

adult.weaned=exp(adult.weaned)/(1+exp(adult.weaned))

# simulating realizations of detection in 1000 individuals with each of the 1000 prob values
adult.weaned=sapply(adult.weaned, function(x) rbinom(n = 1000, prob = x, size = 1))

# detection per simulation
adult.weaned=apply(adult.weaned, MARGIN = 2, mean)


#detection juveniles recently weaned

juvenile.weaned=extract(stan.model.season.age.species.sites.eid)$alpha+
  extract(stan.model.season.age.species.sites.eid)$beta_age+
  extract(stan.model.season.age.species.sites.eid)$beta_rec

juvenile.weaned=exp(juvenile.weaned)/(1+exp(juvenile.weaned))

# simulating realizations of detection in 1000 individuals with each of the 1000 prob values
juvenile.weaned=sapply(juvenile.weaned, function(x) rbinom(n = 1000, prob = x, size = 1))

# detection per simulation
juvenile.weaned=apply(juvenile.weaned, MARGIN = 2, mean)

###################
### FIGURE  3 #####
###################

juvenile.weaned.0.9=HPDI(juvenile.weaned,0.9)
juvenile.weaned.0.51=HPDI(juvenile.weaned,0.51)
juvenile.not.weaned.0.9=HPDI(juvenile.not.weaned,0.9)
juvenile.not.weaned.0.51=HPDI(juvenile.not.weaned,0.51)
adult.weaned.0.9=HPDI(adult.weaned,0.9)
adult.weaned.0.51=HPDI(adult.weaned,0.51)
adult.not.weaned.0.9=HPDI(adult.not.weaned,0.9)
adult.not.weaned.0.51=HPDI(adult.not.weaned,0.51)


juveniles.0.9=cbind(juvenile.weaned.0.9,juvenile.not.weaned.0.9)
juveniles.0.51=cbind(juvenile.weaned.0.51,juvenile.not.weaned.0.51)
adults.0.9=cbind(adult.weaned.0.9, adult.not.weaned.0.9)
adults.0.51=cbind(adult.weaned.0.51, adult.not.weaned.0.51)

colors=rep(c("grey", "lightblue", "pink", "plum3", "yellowgreen", "tan1"),3)
colors2=rep(c("black", "dodgerblue4", "firebrick1", "hotpink4", "forestgreen", "darkorange2"),3)

#labels.species

temp=lapply(strsplit(levels(bats$barcode.sp), split = " "), function(x) 
  paste0(str_sub(x[1], 1,1), ". ", x[-1]))


labels.species=sapply(temp, function(x) if(length(x)>1){
  paste0(
    x[1], " ",
    strsplit(x[2], " ")[[1]][2])}else{x})


#index of pteropodid species
ptero=which(sapply(levels(bats$barcode.sp), function(x) unique(bats[bats$barcode.sp==x,]$Family), USE.NAMES = F)=="Pteropodidae")

# leaving pteropodids a the end
species.output=cbind(extract(stan.model.season.age.species.sites.eid)$theta_species[,-ptero],
                     extract(stan.model.season.age.species.sites.eid)$theta_species[,ptero])


#labels.species

temp=lapply(strsplit(levels(bats$barcode.sp), split = " "), function(x) 
  paste0(str_sub(x[1], 1,1), ". ", x[-1]))


labels.species=sapply(temp, function(x) if(length(x)>1){
  paste0(
    x[1], " ",
    strsplit(x[2], " ")[[1]][2])}else{x})

labels.species=c(labels.species[-ptero], labels.species[ptero])
labels.species[3]="Hipposiderid sp."


# function to mke a chracter vector italic
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))


png('/Fig 3.png',  units = "px", width = 1200, height = 600)
par(mar=c(14, 12, 2, 0), mfrow=c(1,2))

plot (NULL, xlim=c(0.7,2.3), ylim=c(-0.02,1.05), xaxs="i", yaxs="i", axes=F,  xaxt='n',main='', xlab='', ylab="")
axis(side = 2, at = seq(0,1,length.out = 6), labels =seq(0,1,length.out = 6), las=2, cex.axis=2.5)
axis(side = 1, at = 1:2, labels = c("Recent\nWeaning", "Not recent\nweaning"), cex.axis=2.8,  mgp=c(3, 1.6, 0), tck=-0.03, las=2)
mtext(side = 2, text = "Adjusted coronavirus\ndetection", cex=2.8, line=6.7)

for(i in 1:2){
  lines(x = c(i-0.1,i-0.1), y=c(juveniles.0.9[,i]), col=alpha("lightblue", 1), lwd=5)
  lines(x = c(i-0.025-0.1,i+0.025-0.1), y=rep(juveniles.0.9[2,i],2), col=alpha("lightblue", 1), lwd=5) #top whiskers
  lines(x = c(i-0.025-0.1,i+0.025-0.1), y=rep(juveniles.0.9[1,i],2), col=alpha("lightblue", 1), lwd=5) #top whiskers
}

for(i in 1:2){
  lines(x = c(i-0.1,i-0.1), y=c(juveniles.0.51[,i]), col=alpha("deepskyblue3", 1), lwd=10)
  lines(x = c(i-0.05-0.1,i+0.05-0.1), y=rep(juveniles.0.51[2,i],2), col=alpha("deepskyblue3", 1), lwd=10) #top whiskers
  lines(x = c(i-0.05-0.1,i+0.05-0.1), y=rep(juveniles.0.51[1,i],2), col=alpha("deepskyblue3", 1), lwd=10) #top whiskers
}

for(i in 1:2){
  lines(x = c(i+0.1,i+0.1), y=c(adults.0.9[,i]), col=alpha("burlywood1", 1), lwd=5)
  lines(x = c(i+0.1-0.025,i+0.1+0.025), y=rep(adults.0.9[2,i],2), col=alpha("burlywood1", 1), lwd=5) #top whiskers
  lines(x = c(i+0.1-0.025,i+0.1+0.025), y=rep(adults.0.9[1,i],2), col=alpha("burlywood1", 1), lwd=5) #top whiskers
}

for(i in 1:2){
  lines(x = c(i+0.1,i+0.1), y=c(adults.0.51[,i]), col=alpha("darkorange2", 1), lwd=10)
  lines(x = c(i+0.1-0.05,i+0.1+0.05), y=rep(adults.0.51[2,i],2), col=alpha("darkorange2", 1), lwd=10) #top whiskers
  lines(x = c(i+0.1-0.05,i+0.1+0.05), y=rep(adults.0.51[1,i],2), col=alpha("darkorange2", 1), lwd=10) #top whiskers
}

legend(legend = c("Juveniles", "Adults"), fill = c("deepskyblue3", "darkorange2"), x=1.5, y=1.04, cex = 2)

#dev.off()

par(mar=c(14, 8, 2, 0))
plot (NULL, xlim=c(0,2.5), ylim=c(-5.02,5.02), xaxs="i", yaxs="i", axes=F,  xaxt='n',main='', xlab='', ylab="")
axis(side = 2, at = seq(-5,5,length.out = 5), labels =seq(-5,5,length.out = 5), las=2, cex.axis=2.1)
axis(side = 1, at = seq(0.1,2.2, length.out = 13), labels = make.italic(labels.species), cex.axis=2,  mgp=c(3, 1.6, 0), tck=-0.03, las=2)
mtext(side = 2, text = "Species-specific\nintercept" , cex=2.8, line=4)


for(i in 1:n.species){
  temp.0.9=HPDI(species.output[,i], 0.9)
  lines(x = rep(seq(0.1,2.2, length.out = 13)[i],2), y=temp.0.9, col=alpha(colors[i], 1), lwd=5)
  lines(x = c(seq(0.1,2.2, length.out = 13)[i]-0.025,seq(0.1,2.2, length.out = 13)[i]+0.025), y=rep(temp.0.9[2],2), col=alpha(colors[i], 1), lwd=5) #top whiskers
  lines(x = c(seq(0.1,2.2, length.out = 13)[i]-0.025,seq(0.1,2.2, length.out = 13)[i]+0.025), y=rep(temp.0.9[1],2), col=alpha(colors[i], 1), lwd=5) #bottom whiskers
}

for(i in 1:n.species){
  temp.0.51=HPDI(species.output[,i], 0.51)
  lines(x = rep(seq(0.1,2.2, length.out = 13)[i],2), y=temp.0.51, col=alpha(colors2[i], 1), lwd=5)
  lines(x = c(seq(0.1,2.2, length.out = 13)[i]-0.025,seq(0.1,2.2, length.out = 13)[i]+0.025), y=rep(temp.0.51[2],2), col=alpha(colors2[i], 1), lwd=5) #top whiskers
  lines(x = c(seq(0.1,2.2, length.out = 13)[i]-0.025,seq(0.1,2.2, length.out = 13)[i]+0.025), y=rep(temp.0.51[1],2), col=alpha(colors2[i], 1), lwd=5) #bottom whiskers
}

dev.off()


# --- Ratio of detections ---#


#HPDI(juvenile.weaned/juvenile.not.weaned, 0.9)
mean(juvenile.weaned)/mean(juvenile.not.weaned)

#mean(juvenile.weaned)/mean(adult.weaned)


#HPDI(juvenile.weaned/adult.not.weaned, 0.9)


#HPDI(adult.weaned/juvenile.not.weaned, 0.9)


mean(adult.weaned)/mean(adult.not.weaned)

# HPDI(juvenile.weaned/adult.weaned, 0.9)
# HPDI(juvenile.not.weaned/adult.not.weaned, 0.9)


mean(juvenile.weaned)/mean(adult.weaned)
mean(juvenile.not.weaned)/mean(adult.not.weaned)

## TABLE 1 ##

table1=ftable(bats$barcode.sp, bats$Age, bats$new.season)


table2=ftable(bats$barcode.sp, bats$Age, bats$new.season, CoV)

table2=data.frame(Species= c(sapply(attr(table2, "row.vars")[[1]], function(x)rep(x,4), USE.NAMES = F)),
                  Age=rep(c("Adult", "Adult", "Juvenile", "Juvenile"), 13), 
                  Repro=rep(c("Not recently weaned", "Recently weaned"), 26),
                  Pos=table2[,2],
                  tested=rowSums(table2),
                  prop.pos=round(table2[,2]/rowSums(table2), 3))


## -- Imputation model -- ##

summary(stan.model.season.age.species.sites.eid)$summary

HPDI(c(a_imp), 0.9)
HPDI(c(b1_imp), 0.9)
HPDI(c(b2_imp), 0.9)
HPDI(c(b3_imp), 0.9)
HPDI(c(b4_imp), 0.9)

##################################################
# Additional Information 3: Figures imputed data #
##################################################

#-----------------------#
# Additional figure 3.1 #
#-----------------------#

temp.results=extract(stan.model.season.age.species.sites.eid)$x_imp_all

#global
temp.results=extract(stan.model.season.age.species.sites.eid)$x_imp_all[,season_cat_miss==0]
temp.results=t(apply(temp.results, MARGIN = 1, FUN = function(x){season_cat[season_cat_miss==0]-x}))
temp.results.global=prop.table(table(temp.results))

# RW only

temp.results=extract(stan.model.season.age.species.sites.eid)$x_imp_all[,season_cat_miss==0 & season_cat==1]
temp.results=t(apply(temp.results, MARGIN = 1, FUN = function(x){season_cat[season_cat_miss==0 & season_cat==1]-x}))
temp.results.rw=prop.table(table(temp.results))
temp.results.rw=c(NA, temp.results.rw)


# N-RW only

temp.results=extract(stan.model.season.age.species.sites.eid)$x_imp_all[,season_cat_miss==0 & season_cat==0]
temp.results=t(apply(temp.results, MARGIN = 1, FUN = function(x){season_cat[season_cat_miss==0 & season_cat==0]-x}))
temp.results.nrw=prop.table(table(temp.results))
temp.results.nrw=c(temp.results.nrw,NA)


png('/Additional_figure_3_1.png', units = "px", width = 2200, height = 800)
par(mar=c(14,14.5,14,3), mfrow=c(1,3))

rep.season.names=c("Recent weaning", "Not recent weaning")

h=barplot(temp.results.global,  xaxs="i", yaxs="i", axes=F,  xaxt='n',
          main='', xlab='', ylab="", ylim=c(0,1), col="purple")
axis(side = 2, at = seq(0,1,length.out = 6), labels = seq(0,1,length.out = 6), las=2, cex.axis=4.5)
axis(side = 1, at = c(h[,1]), labels =c(-1,0,1) , cex.axis=4.5, mgp=c(3, 5, 0))
mtext(text = 'Global', side = 3, line = 4, cex=5)
mtext(text = 'Proportion', side = 2, line = 10, cex=4)
mtext(text = 'Difference', side =1, line = 12, cex=4)


h=barplot(temp.results.rw,  xaxs="i", yaxs="i", axes=F,  xaxt='n',
          main='', xlab='', ylab="", ylim=c(0,1), col="purple")
axis(side = 2, at = seq(0,1,length.out = 6), labels = seq(0,1,length.out = 6), las=2, cex.axis=4.5)
axis(side = 1, at = c(h[,1]), labels =c(-1,0,1) , cex.axis=4.5, mgp=c(3, 5, 0))
mtext(text = 'Recent weaning', side = 3, line = 4, cex=5)
#mtext(text = 'Proportion', side = 2, line = 10, cex=4)
mtext(text = 'Difference', side =1, line = 12, cex=4)

h=barplot(temp.results.nrw,  xaxs="i", yaxs="i", axes=F,  xaxt='n',
          main='', xlab='', ylab="", ylim=c(0,1), col="purple")
axis(side = 2, at = seq(0,1,length.out = 6), labels = seq(0,1,length.out = 6), las=2, cex.axis=4.5)
axis(side = 1, at = c(h[,1]), labels =c(-1,0,1) , cex.axis=4.5, mgp=c(3, 5, 0))
mtext(text = 'Not recent weaning', side = 3, line = 4, cex=5)
#mtext(text = 'Proportion', side = 2, line = 10, cex=4)
mtext(text = 'Difference', side =1, line = 12, cex=4)


dev.off()


#-----------------------#
#Additional  figure 3.2 #
#-----------------------#

temp.cat=apply(extract(stan.model.season.age.species.sites.eid)$x_imp_unobs, MARGIN = 1, FUN = function(x) sum(x))
hist(temp.cat)
temp.cat=t(apply(extract(stan.model.season.age.species.sites.eid)$x_imp_unobs, MARGIN = 1, FUN = function(x) prop.table(table(x))))

temp.results=t(apply(extract(stan.model.season.age.species.sites.eid)$x_imp_unobs, MARGIN = 1, FUN = function(x) tapply(CoV, list(x), mean)))
dim(temp.results)

apply(temp.results, MARGIN = 2, summary)
prop.table(table(CoV))


# adults
prop.table(table(CoV[bats$Age=="Adult"], bats[bats$Age=="Adult",]$new.season), margin = 2)

# juveniles
prop.table(table(CoV[bats$Age!="Adult"], bats[bats$Age!="Adult",]$new.season), margin = 2)


temp.results.juv=t(apply(extract(stan.model.season.age.species.sites.eid)$x_imp_unobs[,bats$Age!="Adult"], MARGIN = 1, FUN = function(x) tapply(CoV[bats$Age!="Adult"], list(x), mean)))
temp.results.adults=t(apply(extract(stan.model.season.age.species.sites.eid)$x_imp_unobs[,bats$Age=="Adult"], MARGIN = 1, FUN = function(x) tapply(CoV[bats$Age=="Adult"], list(x), mean)))

rep.season.names=rev(c("Recent\nweaning", "Not recent\nweaning"))

png('/Additional_Figure_3_2.png', units = "px", width = 1450, height = 800)

par(mar=c(15, 15, 4, 2), mfrow=c(1,2))


matplot(t(temp.cat), col="black", type="l", ylim=c(0,1), xlim=c(0.65, 2.5), lty=3, xaxs="i", yaxs="i", axes=F,  xaxt='n',
        main='', xlab='', ylab="")
points(x=rep(1, nrow(temp.cat)), y=temp.cat[,1], pch=21, bg="white", col="blue", cex=3)
points(x=rep(2, nrow(temp.cat)), y=temp.cat[,2], pch=21, bg="white", col="red", cex=3)
axis(side = 2, at = seq(0,1,length.out = 6), labels = seq(0,1,length.out = 6), las=2, cex.axis=3)
axis(side = 1, at = c(1:2), labels =rep.season.names , cex.axis=3, mgp=c(3, 2, 0), las=2)
mtext(text = 'Proportion bats\nin posterior dist.', side = 2, line = 8, cex=4)

par(mar=c(15, 15, 4, 2))
boxplot(temp.results.juv, xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="",  border=NA, col=alpha(c("lightblue", "pink"), 1),  ylim=c(-0.01,0.6),  xlim=c(0.65, 2.5))
boxplot(temp.results.adults, add=T, xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="",  border=NA, col=alpha(c("blue", "red"), 1))

boxplot(temp.results,  xaxs="i", yaxs="i", axes=F,  xaxt='n', border=c("black"),
        main='', xlab='', ylab="", col=alpha(c("darkblue", "red4", "darkgreen", "purple"), 0.1), add=T)
axis(side = 2, at = seq(0,0.6,length.out = 4), labels = seq(0,0.6,length.out = 4), las=2, cex.axis=3)
axis(side = 1, at = c(1:2), labels =rep.season.names , cex.axis=3, mgp=c(3, 2, 0), las=2)
mtext(text = 'Period', side = 1, line = 18, cex=4)
mtext(text = 'Coronavirus detection', side = 2, line = 8, cex=4)

dev.off()





##################################################
# Additional Information 4: Model for imputation #
##################################################

# data wo NA in the season to fit the model for the imputation and get the WAIC and LOO

#------#
# DATA #
#------#

# number of observations
N_complete=nrow(bats[season_cat!=(-1),]) 

#lat_complete=c(scale(bats$Latitude[season_cat!=(-1)], center = T)) #latitude
lat_complete=bats$Latitude[season_cat!=(-1)] #latitude

#historical_prep_complete=c(scale(bats$historic.monthly.prep[season_cat!=(-1)], center = T)) #historical precipitation
historical_prep_complete=bats$historic.monthly.prep[season_cat!=(-1)] #historical precipitation


litters_per_year_complete=bats$litters.per.year[season_cat!=(-1)] #his
litters_per_year_complete=ifelse(litters_per_year_complete>1,1,0)

day_of_the_year_complete=bats$Day.of.the.year[season_cat!=(-1)]

season_cat_complete=season_cat[season_cat!=(-1)]


# Stan code for imputation #

modelstring <- "
data{

int N_complete; //number of observations

//repro season covariate
int season_cat_complete[N_complete]; // categorical variable when observed

//data to impute terpoductive season
real day_of_the_year_complete[N_complete];
real lat_complete[N_complete];
real historical_prep_complete[N_complete];
int litters_per_year_complete[N_complete];
}

parameters{
real <lower=-6,upper=6> a_imp; // intercept for the imputation model and in the chunk to model the probability of category 2
real <lower=-6,upper=6> b1_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real <lower=-6,upper=6> b2_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real <lower=-6,upper=6> b3_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
real <lower=-6,upper=6> b4_imp; // coefficent for the imputation model to multiply the covariate predicting the category variable with missing data and in the chunk to model the probability of category 2
}


model{
// priors
a_imp ~ normal(0,1.5); // explained above
b1_imp ~ normal(0,1.5); // explained above
b2_imp ~ normal(0,1.5); // explained above
b3_imp ~ normal(0,1.5); // explained above
b4_imp ~ normal(0,1.5); // explained above


//add Data in the model
for (i in 1:N_complete) {

vector[2] p;
p[1] = a_imp + b1_imp*day_of_the_year_complete[i] + b2_imp*lat_complete[i] + b3_imp*litters_per_year_complete[i] + b4_imp*historical_prep_complete[i]; // modeling the prob 2 as a function of the covariate to model the category
p[2]=1-p[1];

season_cat_complete[i] ~ bernoulli(softmax(p)[1]);
}

} // close model block


generated quantities { // generate estimates of the imputed category for all observations


vector[N_complete] log_lik; 
real dev;
dev = 0;

for (i in 1:N_complete) {
real p;
p = a_imp + b1_imp*day_of_the_year_complete[i] + b2_imp*lat_complete[i] + b3_imp*litters_per_year_complete[i] + b4_imp*historical_prep_complete[i]; // modeling the prob 2 as a function of the covariate to model the category

dev = dev + (-2)*bernoulli_logit_lpmf( season_cat_complete[i] | p);
log_lik[i] = bernoulli_logit_lpmf(season_cat_complete[i] | p);

} // close loop

} // close generated quantities  block

"

#Prepare data

dat <- list(N_complete = N_complete,
            season_cat_complete=season_cat_complete,
            day_of_the_year_complete=c(scale(day_of_the_year_complete,center = T)),
            lat_complete=c(scale(lat_complete,center = T)),
            historical_prep_complete=c(scale(historical_prep_complete,center = T)),
            litters_per_year_complete=litters_per_year_complete)


# Run model 

stan.model.imputation=stan(model_code = modelstring,  
                           iter = 4000, warmup = 3000, thin=1,
                           chains = 4, cores = 4,
                           data = dat,
                           control=list(adapt_delta=0.995, max_treedepth = 15))



# extract relevant information: LOO and WAIC
log_lik1 <- extract_log_lik(stan.model.imputation, merge_chains = F, parameter_name = "log_lik")
rel_n_eff <- relative_eff(exp(log_lik1))
product<-loo(log_lik1, r_eff = rel_n_eff, cores = 2)
imputation.model.waic=waic(log_lik1)
imputation.model.loo=loo(log_lik1, r_eff =rel_n_eff )




##################################################
# Additional Information 5: Model for imputation #
##################################################


#-----------------------#
# Additional figure 5.1 #
#-----------------------#

# General detection #

post.pred.CoV=extract(stan.model.season.age.species.sites.eid)$y_rep
dim(post.pred.CoV)

hist(apply(post.pred.CoV, MARGIN = 1, FUN = mean))
abline(v=mean(CoV))



# Detection per age #

predicted.detection.per.age=apply(post.pred.CoV, MARGIN = 1, FUN = function(x){tapply(x, list(age), mean)})
dim(predicted.detection.per.age)
predicted.detection.per.age=t(predicted.detection.per.age)

#histogram predicted detection in adults
hist(predicted.detection.per.age[,1])
abline(v=tapply(CoV, list(age), mean)[1])

#histogram predicted detection in juveniles
hist(predicted.detection.per.age[,2])
abline(v=tapply(CoV, list(age), mean)[2])



# Detection per reproductive season #

predicted.detection.per.repro.season=apply(post.pred.CoV[,season_cat!=(-1)], MARGIN = 1, FUN = function(x){tapply(x, list(season_cat[season_cat!=(-1)]), mean)})
dim(predicted.detection.per.repro.season)
predicted.detection.per.repro.season=t(predicted.detection.per.repro.season)

#histogram predicted detection birth pulse and pup nursing
hist(predicted.detection.per.repro.season[,1])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)]), mean)[1])

hist(predicted.detection.per.repro.season[,2])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)]), mean)[2])



# Detection per age per reproductive season #

predicted.detection.per.age.per.repro.season=lapply(c(1:4000), function(x){tapply(post.pred.CoV[x,season_cat!=(-1)], list(season_cat[season_cat!=(-1)], age[season_cat!=(-1)]), mean)})

predicted.detection.per.age.per.repro.season=simplify2array(predicted.detection.per.age.per.repro.season)
dim(predicted.detection.per.age.per.repro.season)

# detection in adults Non-recnet weamimg
hist(predicted.detection.per.age.per.repro.season[1,1,])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)], age[season_cat!=(-1)]), mean)[1,1])

# detection in juveniles  Non-recnet weamimg
hist(predicted.detection.per.age.per.repro.season[1,2,])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)], age[season_cat!=(-1)]), mean)[1,2])

# detection in adults during recent weaing
hist(predicted.detection.per.age.per.repro.season[2,1,])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)], age[season_cat!=(-1)]), mean)[2,1])


# detection in juveniles during recent weaning
hist(predicted.detection.per.age.per.repro.season[2,2,])
abline(v=tapply(CoV[season_cat!=(-1)], list(season_cat[season_cat!=(-1)], age[season_cat!=(-1)]), mean)[2,2])





png('/posterior_predictions_and_observed.png',  units = "px", width = 1000, height = 1200)
par(mar=c(10, 12.5, 3, 3), mfrow=c(3,2))

# general detection

d=hist(apply(post.pred.CoV, MARGIN = 1, FUN = mean), plot = FALSE)
hist(apply(post.pred.CoV, MARGIN = 1, FUN = mean),  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, signif(max(d$counts), 3))+50, xlim=c(0,0.8), col="cornflowerblue")
axis(side = 1, at = seq(-0,1,length.out = 6), labels = seq(-0,1,length.out = 6), las=1, cex.axis=3.5, mgp=c(4, 3, 0), lwd.ticks = 3, tck=-0.025)
axis(side = 2, at = seq(0,signif(max(d$counts), 2), length.out = 5), labels = seq(0,signif(max(d$counts), 2), length.out = 5)/100, las=1, cex.axis=3.5, mgp=c(2, 1.5, 0), lwd.ticks = 3, tck=-0.025)
mtext(text ="Frequency", side = 2, line = 9.5, cex=3)
text(x=0.875, y=signif(max(d$counts), 2), "A", cex=4, pos = 1)
abline(v=mean(CoV), lwd=4, col="black")

# Detection per age

par(mar=c(10, 9, 3, 3))
e=hist(predicted.detection.per.age[,1], plot = FALSE)
f=hist(predicted.detection.per.age[,2], plot = FALSE)
hist(predicted.detection.per.age[,1],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, round_any(max(c(e$counts, f$counts)), 100, f = ceiling)), xlim=c(0,0.8), col=scales::alpha("cornflowerblue", 0.3))
axis(side = 1, at = seq(-0,1,length.out = 6), labels = seq(-0,1,length.out = 6), las=1, cex.axis=3.5, mgp=c(4, 3, 0), lwd.ticks = 3, tck=-0.025)
axis(side = 2, at = seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5), 
     labels = seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5)/100, las=1, cex.axis=3.5, mgp=c(2, 1.5, 0), lwd.ticks = 3, tck=-0.025)
#mtext(text ="Frequency", side = 2, line = 9.5, cex=3)
text(x=0.875, y=round_any(max(c(e$counts, f$counts)), 100, f = ceiling), "B", cex=4, pos = 1)
abline(v=tapply(CoV, list(age), mean)[1], lwd=4, col="navyblue")

hist(predicted.detection.per.age[,2],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, signif(max(f$counts), 3)), xlim=c(0,1), col=scales::alpha("lightgoldenrod1", 0.5), add=T)
abline(v=tapply(CoV, list(age), mean)[2], lwd=4, col="coral1")



# Detection per reproductive season

par(mar=c(10, 12.5, 3, 3))
e=hist(predicted.detection.per.repro.season[,1], plot = FALSE)
f=hist(predicted.detection.per.repro.season[,2], plot = FALSE)
#g=hist(predicted.detection.per.repro.season[,3], plot = FALSE)

hist(predicted.detection.per.repro.season[,1],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, round_any(max(c(e$counts, f$counts)), 100, f = ceiling)), xlim=c(0,0.8), col=scales::alpha("cornflowerblue", 0.3))
axis(side = 1, at = seq(-0,1,length.out = 6), labels = seq(-0,1,length.out = 6), las=1, cex.axis=3.5, mgp=c(4, 3, 0), lwd.ticks = 3, tck=-0.025)
axis(side = 2, at = seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5), 
     labels = round(seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5)/100,1), las=1, cex.axis=3.5, mgp=c(2, 1.5, 0), lwd.ticks = 3, tck=-0.025)
mtext(text ="Frequency", side = 2, line = 9.5, cex=3)
text(x=0.875, y=round_any(max(c(e$counts, f$counts)), 100, f = ceiling), "C", cex=4, pos = 1)
abline(v=tapply(CoV, list(season_cat), mean)[2], lwd=4, col="navyblue")

hist(predicted.detection.per.repro.season[,2],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, signif(max(f$counts), 3)), xlim=c(0,0.8), col=scales::alpha("lightgoldenrod1", 0.5), add=T)
abline(v=tapply(CoV, list(season_cat), mean)[3], lwd=4, col="coral1")



# Detection per age per reproductive season 
#juveniles

par(mar=c(10, 9, 3, 3))
e=hist(predicted.detection.per.age.per.repro.season[1,2,], plot = FALSE)
f=hist(predicted.detection.per.age.per.repro.season[2,2,], plot = FALSE)

hist(predicted.detection.per.age.per.repro.season[1,2,],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, round_any(max(c(e$counts, f$counts)), 100, f = ceiling)), xlim=c(0,0.8), col=scales::alpha("cornflowerblue", 0.3))
axis(side = 1, at = seq(-0,1,length.out = 6), labels = seq(-0,1,length.out = 6), las=1, cex.axis=3.5, mgp=c(4, 3, 0), lwd.ticks = 3, tck=-0.025)
axis(side = 2, at = seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5), 
     labels = round(seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5)/100,1), las=1, cex.axis=3.5, mgp=c(2, 1.5, 0), lwd.ticks = 3, tck=-0.025)
text(x=0.875, y=round_any(max(c(e$counts, f$counts)), 100, f = ceiling), "D", cex=4, pos = 1)
abline(v=tapply(CoV, list(season_cat, age), mean)[2,2]+0.01, lwd=4, col="navyblue")

hist(predicted.detection.per.age.per.repro.season[2,2,],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, signif(max(f$counts), 3)), xlim=c(0,0.8), col=scales::alpha("lightgoldenrod1", 0.5), add=T)
abline(v=tapply(CoV, list(season_cat, age), mean)[3,2], lwd=4, col="coral1")
mtext(text ="Coronavirus detection", side = 1, line = 8, cex=3)


#adults

par(mar=c(10, 12.5, 3, 3))
e=hist(predicted.detection.per.age.per.repro.season[1,1,], plot = FALSE)
f=hist(predicted.detection.per.age.per.repro.season[2,1,], plot = FALSE)

hist(predicted.detection.per.age.per.repro.season[1,1,],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, round_any(max(c(e$counts, f$counts)), 100, f = ceiling)), xlim=c(0,0.8), col=scales::alpha("cornflowerblue", 0.3))
axis(side = 1, at = seq(-0,1,length.out = 6), labels = seq(-0,1,length.out = 6), las=1, cex.axis=3.5, mgp=c(4, 3, 0), lwd.ticks = 3, tck=-0.025)
axis(side = 2, at = seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5), 
     labels = round(seq(0,round_any(max(c(e$counts, f$counts)), 100, f = ceiling), length.out = 5)/100,1), las=1, cex.axis=3.5, mgp=c(2, 1.5, 0), lwd.ticks = 3, tck=-0.025)
mtext(text ="Frequency", side = 2, line = 9.5, cex=3)
text(x=0.875, y=round_any(max(c(e$counts, f$counts)), 100, f = ceiling), "E", cex=4, pos = 1)
abline(v=tapply(CoV, list(season_cat, age), mean)[2,1]+0.01, lwd=4, col="navyblue")

hist(predicted.detection.per.age.per.repro.season[2,1,],  xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", ylim=c(0, signif(max(f$counts), 3)), xlim=c(0,0.8), col=scales::alpha("lightgoldenrod1", 0.5), add=T)
abline(v=tapply(CoV, list(season_cat, age), mean)[3,1], lwd=4, col="coral1")

mtext(text ="Coronavirus detection", side = 1, line = 8, cex=3)


dev.off()




#-----------------------#
# Additional figure 5.2 #
#-----------------------#

bats$barcode.sp=factor(bats$barcode.sp)
#labels.species

temp=lapply(strsplit(levels(bats$barcode.sp), split = " "), function(x) 
  paste0(str_sub(x[1], 1,1), ". ", x[-1]))


labels.species=sapply(temp, function(x) if(length(x)>1){
  paste0(
    x[1], " ",
    strsplit(x[2], " ")[[1]][2])}else{x})

labels.species[c(4)]="Hipposiderid sp."

### PLOTS  ####
temp.names=c(names(extract(stan.model.season.age.species.sites.eid))[1:4], paste0(names(extract(stan.model.season.age.species.sites.eid))[16], ".", c(1:13)), paste0(names(extract(stan.model.season.age.species.sites.eid))[17], ".", c(1:30)))

# PPD densities
png('Model_coefficient_PPDs1.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 1:9){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E. helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()


# PPD densities 2
png('/Micro_final_model_coefficient_PPDs2.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 10:18){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E. helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()


# PPD densities 3
png('/Micro_final_model_coefficient_PPDs3.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 19:27){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E. helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()


# PPD densities 4
png('/Micro_final_model_coefficient_PPDs4.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 28:36){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E. helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()



# PPD densities 5

png('/Micro_final_model_coefficient_PPDs5.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 37:45){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E. helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()



# PPD densities 6

png('Micro_final_model_coefficient_PPDs6.png',  units = "px", width = 1200, height = 1200)
par(mar=c(12, 4, 3, 3), mfrow=c(3,3))

for(i in 46:47){
  plot(density(get(temp.names[i])), xaxs="i", yaxs="i", axes=F,  xaxt='n', main='', xlab='', ylab="", xlim=c(-10,10))
  polygon(density(get(temp.names[i])), col="cornflowerblue", border="black")
  axis(side = 1, at = seq(-10,10,length.out = 11), labels = seq(-10,10,length.out = 11), las=1, cex.axis=3.5, mgp=c(5, 3.3, 0), lwd.ticks = 3, tck=-0.025)
  
  mtext(text = c(expression(alpha[0]), "Age", "RW", "E.helvum - T. persicus", make.italic(labels.species), paste("Sampling site", 1:30, sep=' '))[i], side = 1, line = 9, cex=3)
  abline(v=0)
}

dev.off()




##########################################################
# Additional Information 6: Sampling specific intercepts #
##########################################################


# leaving pteropodids a the end
sites.output=cbind(extract(stan.model.season.age.species.sites.eid)$theta_site)

colors=rep(colors,2)
colors2=rep(colors2,2)

png('Supp_6.png',  units = "px", width = 1200, height = 600)

par(mar=c(11, 11, 2, 0), mfrow=c(1,1))

plot (NULL, xlim=c(0,6.4), ylim=c(-4.02,4.02), xaxs="i", yaxs="i", axes=F,  xaxt='n',main='', xlab='', ylab="")
axis(side = 2, at = seq(-4,4,length.out = 5), labels =seq(-4,4,length.out = 5), las=2, cex.axis=2.5)
axis(side = 1, at = seq(0.1,6.1, length.out = 30), labels = paste("SE", c(1:30), sep=' '), cex.axis=2.8,  mgp=c(3, 1.6, 0), tck=-0.03, las=2)
mtext(side = 2, text = "Sampling event\ncoefficient" , cex=2.8, line=4.7)


for(i in 1:30){
  temp.0.9=HPDI(sites.output[,i], 0.9)
  lines(x = rep(c(seq(0.1,6.1, length.out = 30)[i]),2), y=temp.0.9, col=alpha(colors[i], 1), lwd=5)
  lines(x = c(seq(0.1,6.1, length.out = 30)[i]-0.025,seq(0.1,6.1, length.out = 30)[i]+0.025), y=rep(temp.0.9[2],2), col=alpha(colors[i], 1), lwd=5) #top whiskers
  lines(x = c(seq(0.1,6.1, length.out = 30)[i]-0.025,seq(0.1,6.1, length.out = 30)[i]+0.025), y=rep(temp.0.9[1],2), col=alpha(colors[i], 1), lwd=5) #bottom whiskers
}

for(i in 1:30){
  temp.0.51=HPDI(sites.output[,i], 0.51)
  lines(x = rep(c(seq(0.1,6.1, length.out = 30)[i]),2), y=temp.0.51, col=alpha(colors2[i], 1), lwd=5)
  lines(x = c(seq(0.1,6.1, length.out = 30)[i]-0.025,seq(0.1,6.1, length.out = 30)[i]+0.025), y=rep(temp.0.51[2],2), col=alpha(colors2[i], 1), lwd=5) #top whiskers
  lines(x = c(seq(0.1,6.1, length.out = 30)[i]-0.025,seq(0.1,6.1, length.out = 30)[i]+0.025), y=rep(temp.0.51[1],2), col=alpha(colors2[i], 1), lwd=5) #bottom whiskers
}


dev.off()
