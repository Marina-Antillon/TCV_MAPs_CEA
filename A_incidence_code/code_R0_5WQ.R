library("Hmisc")
library('deSolve')
library("MASS")
library("matrixcalc")
library("emdbook")

library("RColorBrewer")
library("ggplot2")

library("epitools")
library("LaplacesDemon") # so it will give me a dirichlet distribution

library("tidyr")
library("dplyr")
library("abind")

library("readxl")
library("R.matlab")

fls = list.files("./functions/", pattern="^[fcn]")
for (i in 1:length(fls)){source(paste0("./functions/", fls[i]))}

# The cea-calculator is now a function, and 
# the file has the citation/link to the original.
# source("./code/functions/icer_calculator-master/icer_calculator.R")

# make directory to deposit output
subdir = "out" 
if (!dir.exists(file.path("../", subdir))){
  dir.create(file.path("../", subdir))
}

#***************************************************************
## Step 0: Bring in data from Matlab ----
#***************************************************************

source("./functions/fcn_ci_elegant.R")

meanpi = function(x){c(mean=mean(x), lci=quantile(x, 0.025), hci=quantile(x, 0.975))}

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
avgage_yale = readMat("./data/burden_dynamic_link_Yale.mat")$avgage.country
avgage_ihme = readMat("./data/burden_dynamic_link_IHME.mat")$avgage.country

agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv") # MWI: 90, NPL: 96

indicators = read_xlsx("./data/DHS_MICS_GlobalAnalysis_updated.xlsx", "data_to_r")
indicators$ISO = countrycode::countrycode(indicators$Country, "country.name", "iso3c")
indicators$ISO[indicators$Country=="Kosovo"] = "KSV"
indicators$ISO[indicators$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

# Country data from MMGH
Countries <- read_excel("./data/Countries.xlsx", sheet = "Country Data w intro date")
colnames(Countries) = gsub(" - ", "_",colnames(Countries))
colnames(Countries) = gsub(" ", "_",colnames(Countries))
colnames(Countries) = gsub("ISO_code", "ISO",colnames(Countries))
colnames(Countries) = gsub("-", "",colnames(Countries))
colnames(Countries) = gsub("\\(", "",colnames(Countries))
colnames(Countries) = gsub(")", "",colnames(Countries))
colnames(Countries) = gsub("WB_Group_June_2021", "WB_group",colnames(Countries))
Countries$NS_intro = sapply(1:183, function(i){max(Countries[i, c("Low_intro_date", "Medium_Intro_date", "High_Intro_date")])})
Countries$NS_intro[Countries$NS_intro==0 & Countries$WB_group %in% c("Low income", "Lower middle income")] = 2030
Countries$NS_intro[Countries$NS_intro==0 & Countries$WB_group == "Upper middle income"] = 2030
#Continent TODO: change EMRO to respective continents
Countries = Countries %>% 
  mutate(continent = ifelse(Countries$WHO %in% c("WPRO", "SEARO"), "Asia", 
                            ifelse(Countries$WHO=="EMRO", "Mideast",
                                   ifelse(Countries$WHO=="AFRO", "Africa", 
                                          ifelse(Countries$WHO=="EURO", "Eurasia", "Americas")))))
Countries$ISO[Countries$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

no_imp_water = read.csv("./data/share-without-improved-water_OWID_JMP.csv") %>% 
  filter(Year==2020)

no_imp_sanitation = read.csv("./data/sanitation_service_level_JMP.csv") %>% 
  filter(Year==2020, Service.level=="Unimproved")

indicators = right_join(indicators, data.frame(ISO=Countries$ISO, WB_group=Countries$WB_group, continent=Countries$continent)) %>% 
  filter(WB_group!="High income") %>% 
  mutate(across(contains("sani_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("water_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("mcv1_wq"), ~as.numeric(.x))) %>% 
  mutate_at(c("mcv1_total", "sani_total", "water_total"), as.numeric)
  
for(k in which(is.na(indicators$sani_total))){
  indicators$sani_total[k] = ifelse(indicators$ISO[k] %in% no_imp_sanitation$ISO3,
                                     100-no_imp_sanitation$Coverage[no_imp_sanitation$ISO3==indicators$ISO[k]],
                                     99.9)
}
indicators$sani_total[indicators$sani_total>99.9] = 99.9

indicators$sani_dif = indicators$sani_wq5-indicators$sani_wq1
indicators$sani_difa = indicators$sani_total-indicators$sani_wq1
indicators$sani_difb = indicators$sani_wq5-indicators$sani_total

indicators = indicators %>% 
  mutate(sani_wq1 = ifelse(is.na(sani_wq1), pmax(sani_total-20, 0.1), sani_wq1),
         sani_wq2 = ifelse(is.na(sani_wq2), pmax(sani_total-10, 0.1), sani_wq2),
         sani_wq3 = ifelse(is.na(sani_wq3), pmin(sani_total, 99.9), sani_wq3),
         sani_wq4 = ifelse(is.na(sani_wq4), pmin(sani_total+10, 99.9), sani_wq4),
         sani_wq5 = ifelse(is.na(sani_wq5), pmin(sani_total+20, 99.9), sani_wq5)) %>% 
  mutate(across(contains("sani_wq"), ~pmax(pmin(.x,99.9),0.1))) 

#*************************************************************
## Step 1: Fit R0 for the different WQ ----
#*************************************************************
# countries for today: Pakistan, Nigeria, Kenya, if time allows, Malawi
# country numbers: 97, 94, 62, 90

I2toCtf = F # option for secondary infections to go to be eligible for chronic carriage: no

# R0 = 3.53 # India; get rep: 0.15 and get observed incidence (and avg age...)
# R0_WQ = c(3.00, 3.25, 3.53, 3.75, 4)
num_wq = 5

birthrates = c(36.6, 23.6, 15.0)

# for (z in 1:10){ # 97, 94, 62, 90
  
  # The number of the run in the config file which is being run.
  # This will be set automatically if the script is run from slurm.
  zz <- strtoi(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=Sys.getenv("LINE_NUMBER", unset="1")))
  z = which(which(ISO$countryiso==indicators$ISO[zz])==CN)
  # Should line up with the 133 countries of the indicators.

  subdir = paste0("../out/", ISO[CN[z],2])
  if (!dir.exists(file.path("./", subdir))){
    dir.create(file.path("./", subdir))
  }
  
  if(agedist[z]==1){
    #These are only these ages: 0-4, 5-9, 10-14, 15-20, 20-25. The rest are
    #1-sum(distribution)
    low_inc = c(0.158914104, 0.14116809, 0.124802931, 0.108397066, 0.091790805)
    py = c(low_inc[1]/5, low_inc[1]/5, low_inc[1]*3/5, sum(low_inc[2:3]), sum(low_inc[4:5]), 1-sum(low_inc))*1e5
  }else if(agedist[z]==2){
    lmid_inc = c(0.108066773, 0.103190497, 0.098970907, 0.094718397, 0.090066389)
    py = c(lmid_inc[1]/5, lmid_inc[1]/5, lmid_inc[1]*3/5, sum(lmid_inc[2:3]), sum(lmid_inc[4:5]), 1-sum(lmid_inc))*1e5
  }else{
    umid_inc = c(0.072224644, 0.068942809, 0.066816143, 0.068359158, 0.080071929)
    py = c(umid_inc[1]/5, umid_inc[1]/5, umid_inc[1]*3/5, sum(umid_inc[2:3]), sum(umid_inc[4:5]), 1-sum(umid_inc))*1e5
  }

  # with WQ ----
  parameters = list(delta=1/4, # rate of recovery
                    alpha=0.01, # 0.01, # prob dying because of typhoid
                    omega=1/104, # duration of immunity. Alternative: -log(1-0.25)/52
                    # see article on typhoid outbreaks in British troops two years in a row.
                    theta = rep(c(rep(0.003, 4), 0.021, 0.021), num_wq), # prob I1 becomes chronic. c(rep(0.003, 4), 0.021, 0.082)
                    theta2 = 0, # prob I2 becomes chronic (this is the focus of one kind of sensitivity analysis)
                    epsilon = 0,
                    # r = 0.24, # contribution of secondary carriers to transmission
                    # rC = 0.24, # contribution of chronic carriers to transmission
                    mub = log(1+birthrates[agedist[z]]/1000)/52, # birthrate, transformed to account for exponential growth.
                    # mu = -log(1-data$mort_per_py)/52,
                    u = rep(c(12/9, 12/15, 12/35, 1/10, 1/10, 0)/52, num_wq),
                    v = rep(c(0, 12/9, 12/15, 12/35, 1/10, 1/10)/52, num_wq),
                    al = 6, # number of age groups is equivalent to the number of rows in "data" table
                    py = py,
                    population=sum(py),
                    num_wq = num_wq,
                    pop_wq=rep(1/num_wq, num_wq))

  aging_in = c(parameters$mub*sum(py), parameters$u[1:(parameters$al-1)]*py[1:(parameters$al-1)])
  aging_out = parameters$u[1:parameters$al]*py
  
  parameters$mu = rep((aging_in - aging_out)/parameters$py, parameters$num_wq)
  # parameters$mu
  # [1] -0.00092951  0.00000000  0.00000000 -0.00008146 -0.00004802  0.00062887
  
  data = list()
  data$py = py
  # data$min_rep = min(TransPar$repsamples[,56])
  # data$max_rep = max(TransPar$repsamples[,56])
  # 
  # data$mean_rep = mean(TransPar$repsamples[,56]) # /TransPar$repsamples[,56])
  # data$sd_rep = sd(TransPar$repsamples[,56]) # /TransPar$repsamples[,56])
  
  # data$reporting = TransPar$repsamples[,z]
  
  data$min_inc = min(TransPar$incsamples[,z]) # /TransPar$repsamples[,z]
  data$max_inc = max(TransPar$incsamples[,z]) # /TransPar$repsamples[,56])
  
  data$mean_inc = mean(TransPar$incsamples[,z]) # /TransPar$repsamples[,56])
  data$sd_inc = sd(TransPar$incsamples[,z]) # /TransPar$repsamples[,56])
  data$sd_inc = ifelse(is.na(data$sd_inc)|data$sd_inc==0, 0.005*data$mean_inc, data$sd_inc) # VCT didn't have uncertainty because IHME didn't estimate it 
  
  # data$repinc_Sigma = cov(data.frame(inc=TransPar$incsamples[,56], rep=TransPar$repsamples[,56]))
  
  # OR 1.56 (1.25 – 1.95) for unsafe waste management
  # put in access to improved water.
  data$oddsrat = rlnorm(10000, meanlog = log(1.56), sdlog = (log(1.56)-log(1.25))/2)
  data$haves = as.numeric(indicators[zz, paste0("sani_wq", 1:5)])
  # india: c(37.5, 57.5, 73.3, 86.9, 95.6)
  
  odds_wq = (rep(1, 10000) %o% data$haves/100) + (data$oddsrat %o% (100-data$haves)/100)
  odds_wq = sweep(odds_wq, 1, rowSums(odds_wq), "/") # how to make it dirichlet... see Briggs.
  # dirichlet: 
  data$dirichletalphas = (mean(odds_wq[,1])*(1-mean(odds_wq[,1]))/var(odds_wq[,1])-1)*colMeans(odds_wq)/sum(colMeans(odds_wq))
  
  data$m1_loc = mean(qlogis(TransPar$m1))
  data$m2_loc = mean(qlogis(TransPar$m2))
  data$m_Sigma = cov(data.frame(m1=qlogis(TransPar$m1), m2=qlogis(TransPar$m2)))
  data$m_Sigma = data$m_Sigma/10
  # the variance is smaller, but so is the covariance between the two.

  data$rC_loc = mean(qlogis(TransPar$rC.sample))
  data$rC_sca = sqrt(var(qlogis(TransPar$rC.sample))*(3/pi^2)) #/10
  
  data$rep_loc = mean(qlogis(pmin(TransPar$repsamples[,z], 0.999)))
  data$rep_sca = sqrt(var(qlogis(pmin(TransPar$repsamples[,z], 0.999)))*(3/pi^2))
  data$rep_sca = ifelse(data$rep_sca<0.20, data$rep_sca, 0.2)
  data$rep_sca = ifelse(data$rep_sca>0, data$rep_sca, 0.01)
  data$reporting = pmin(TransPar$repsamples[,z], 0.999)
  
  data$mdpts_age = c(4.5/12, 0.5*(24+9)/12, 3.5, 10, 20, 45)
  data$mean_age = mean(c(avgage_yale[z], avgage_ihme[z]), na.rm=T)
  data$sd_age = 0.25*abs(avgage_yale[z] - avgage_ihme[z]) # so range = 95% CI
  data$sd_age = ifelse(is.na(data$sd_age)|data$sd_age==0, 0.005*data$mean_age, data$sd_age) # VCT didn't have uncertainty because IHME didn't estimate it 
    
  # parameters$m1 = mean(TransPar$m1)
  # parameters$m2 = mean(TransPar$m2)
  # parameters$rC = mean(TransPar$rC.sample)
    
  R0init = quantile(TransPar$R0samples[,z], 0.025)

  SIR_fit = fit_mod_wq(data, parameters, R0init, R0minit=3)
  # SAVE
  save(SIR_fit, data, file=paste0("../out/", ISO[CN[z],2], "/fit_global.Rdata"))
  
