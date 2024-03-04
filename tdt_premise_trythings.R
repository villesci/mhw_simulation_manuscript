
require(tidyverse)
require(lubridate)
require(tibble)
require(ggplot2)
require(purrr)
require(ggridges)
require(parallel)
require(foreach)
require(data.table)
require(zoo)
require(cowplot)
require(doParallel)
require(rlang)
require(reshape2)
require(heatwaveR)
require(scico)
require(ggpubr)
require(SpatialTools)


options(scipen=999)

source("mhwsim/Thermal_landscape_functions.R")
source("mhwsim/hour_fxns_heatwaveR.R")


merged_temp_tri_alphas<-readRDS("mhwsim/data/merged_temp_tri_alphas.RDs") #reading this file in can take a while




set.seed(123)

#Let's create three species: Two of which have intersecting TDT at 30C and 15 days. The third will intersect the high CT species at 30C at 1 day. We will set CTmax

#first, find potential values of CTmax to work with that intsersect at 15 days at 30C
temp.seq<-seq(35,50,by=1)
z_list<-vector("list",length=length(temp.seq))

for (i in 1:length(temp.seq)){
  slope<-(log10(15*60*24-1)/(30-temp.seq[i]))
  z_list[[i]]<-(-1/slope)
}

z_bound<-unlist(z_list)
sim_vals<-data.frame("ctmax"=temp.seq,"z"=z_bound)

#we select 48C and 44C to start off
highct_highz<-sim_vals[14,]
highct_highz$species<-"highct_highz"
lowct_midz<-sim_vals[10,]
lowct_midz$species<-"lowct_midz"

#thus far, we have selected tdt curves that intersect in thè middle of our parameter field (15 days, 4C above)
#we now are interestd in a species that still has a lower CTmax, but will have a lower z such that it intersects the highCTmax Highz species quickly at 12 hours.

equaltdt <- function(CTmax1, CTmax2, z1, z2, timetodeath) {
  if (is.na(CTmax2)) {
    CTmax2 <- CTmax1 - (log10(timetodeath) * (z1 - z2))
    tempatdeath <- CTmax2 - z2 * log10(timetodeath)
    return(list(CTmax2 = CTmax2, z2 = z2, temp_at_equal_TDT = tempatdeath))
  } else if (is.na(z2)) {
    z2 <- z1 - ((CTmax1 - CTmax2) / log10(timetodeath))
    tempatdeath <- CTmax2 - z2 * log10(timetodeath)
    return(list(z2 = z2, temp_at_equal_TDT = tempatdeath))
  } else if (is.na(timetodeath)) {
    timetodeath<-10^((CTmax1-CTmax2)/(z1-z2))
    tempatdeath <- CTmax2 - z2 * log10(timetodeath)
    
    return(list(timetodeath=timetodeath,tempatdeath=tempatdeath))
  }
}

#acute-chronic
equ_res1a<-equaltdt(CTmax1=highct_highz$ctmax,CTmax2=lowct_midz$ctmax,z1=highct_highz$z,z2=NA,timetodeath=12*60)

lowct_lowz<-data.frame("ctmax"=44,"z"=equ_res1a$z2,"species"="lowct_lowz")
equ_res1b<-equaltdt(CTmax1=lowct_midz$ctmax,CTmax2=highct_highz$ctmax,z1=lowct_midz$z,z2=highct_highz$z,timetodeath=NA)


#acute-mixed
equ_res2<-equaltdt(CTmax1=lowct_lowz$ctmax,CTmax2=lowct_midz$ctmax,z1=lowct_lowz$z,z2=lowct_midz$z,timetodeath=NA)
#chronic-mixed
equ_res3<-equaltdt(CTmax1=highct_highz$ctmax,CTmax2=lowct_lowz$ctmax,z1=highct_highz$z,z2=lowct_lowz$z,timetodeath=NA)


#set ctmax and z parameters
example_TDT<-bind_rows(highct_highz,lowct_midz,lowct_lowz)


temp<-rep(seq(26,46,by=2))

# Create an empty data frame to store the results
sim_pred <- vector("list",length=nrow(example_TDT))
error<-10^(-1)

# Loop over the species
for(m in 1:nrow(example_TDT)){
  # Model of TDT curve
  time <- 10^((example_TDT$ctmax[m] - temp) / example_TDT$z[m])
  
  fit <- lm(log10(time) ~ temp)
  
  
  
  sim_results <- data.frame(species = rep(example_TDT$species[m], length(time)),
                            temp = temp,
                            time = time,
                            geom = "line")
  simulated_data <- data.frame()
  
  for (j in seq_along(temp)) {
    # Generate 10 replicates for each temperature with normal distribution
    replicates <- data.frame(
      temp = rep(temp[j], times = 10),
      time = 10^(rnorm(10, mean(predict(fit, newdata = data.frame(temp = temp[j]))), sqrt(error)))
    )
    # Add to the simulated_data
    replicates$species <- rep(example_TDT$species[m], nrow(replicates))
    replicates$geom="point"
    simulated_data <- bind_rows(simulated_data, replicates)
    
  }
  
  sim_pred[[m]] <- rbind(simulated_data, sim_results)
}

sim_pred_bound<-bind_rows(sim_pred)

simulated_line<-sim_pred_bound%>%filter(geom=="line")
simulated_point<-sim_pred_bound%>%filter(geom=="point")

temp.df<-simulated_point%>%select(-"geom")


#we want both examples to have almost no mortality at no MHW
unique_simsp<-unique(sim_pred_bound$species)


models <- by(simulated_point, simulated_point$species, function(sub_df) lm(log10(time) ~ temp, data = sub_df))




ggplot()+
  geom_line(data=simulated_line,aes(x=temp,y=time,color=species))+
  scale_y_log10()+
  geom_point(data=simulated_point,aes(x=temp,y=time,color=species))+
  theme_classic()+
  ylab("Time to death (minutes)")+xlab("Assay Temperature (°C)")+
  scale_color_manual(name="Species Strategy",
                     labels=c("Acute Tolerator","Chronic Strategy","Mixed Strategy"),
                     values=c("royalblue","tomato3","forestgreen"))+  theme(legend.text.align = 0)+
  scale_x_continuous(breaks=seq(26,46,by=4))+
  
  #chronic-acute
  geom_point(aes(x=equ_res1b$tempatdeath,y=equ_res1b$timetodeath))+
  geom_segment(aes(x=25,xend=equ_res1b$tempatdeath,y=equ_res1b$timetodeath,yend=equ_res1b$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  geom_segment(aes(x=equ_res1b$tempatdeath,xend=equ_res1b$tempatdeath,y=0,yend=equ_res1b$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  annotate("text",x=27,y=8000,label=paste0(round((equ_res1b$timetodeath/24),digits=1)," Hours, \n ",equ_res1b$tempatdeath,"°C "),
           size=4)+
  #chronic-mixed
  geom_point(aes(x=equ_res2$tempatdeath,y=equ_res2$timetodeath))+
  geom_segment(aes(x=25,xend=equ_res2$tempatdeath,y=equ_res2$timetodeath,yend=equ_res2$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  geom_segment(aes(x=equ_res2$tempatdeath,xend=equ_res2$tempatdeath,y=0,yend=equ_res2$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  annotate("text",x=41,y=0.5,label=paste0(round((equ_res2$timetodeath/24),digits=1)," Hours,  ",equ_res2$tempatdeath,"°C "),
           size=4)+
  #acute-mixed
  geom_point(aes(x=equ_res3$tempatdeath,y=equ_res3$timetodeath))+
  geom_segment(aes(x=25,xend=equ_res3$tempatdeath,y=equ_res3$timetodeath,yend=equ_res3$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  geom_segment(aes(x=equ_res3$tempatdeath,xend=equ_res3$tempatdeath,y=0,yend=equ_res3$timetodeath),color="black",linetype="dashed",inherit.aes=F)+
  annotate("text",x=38,y=9000,label=paste0(round((equ_res3$timetodeath/24),digits=1)," Hours,  ",round(equ_res3$tempatdeath,digits=1),"°C "),
           size=4)


merged_temp_tri<-merged_temp_tri_alphas[[2]]

rm(merged_temp_tri_alphas)
#   user  system elapsed 
#  71.71   60.69 16739  


unique_simsp<-unique(temp.df$species)
final_surv_crosssim_species<-vector("list",length(unique_simsp))
r<-vector("list",length(merged_temp_tri))

numCores <- detectCores()
cl <- makeCluster(numCores-4)
registerDoParallel(cl)

sys.tm <- system.time({
  for (yy in 1:length(unique_simsp)) {
    temp.df_sp <- temp.df %>% dplyr::filter(species == unique_simsp[yy])
    
    results <- foreach::foreach(i = 1:length(merged_temp_tri), .packages = c('tidyverse', 'ggplot2'), .combine = 'rbind') %dopar% {
      tmpi <- rezende_min(merged_temp_tri[[i]], temp.df_sp,tc=26)
      temp_plat_rep <- purrr::map_dfr(seq_len(2), ~merged_temp_tri[[i]])
      r_subset <- cbind(tmpi, temp_plat_rep)
      
      r_subset <- r_subset %>% dplyr::filter(type == "cum.mort")
      
      final_surv_crosssim_species_n <- data.frame(
        'cum.mort' = min(na.omit(r_subset$mort)),
        "magnitude" = max(r_subset$magnitude),
        "duration" = max(r_subset$duration),
        "rising_slope" = max(r_subset$rising_slope),
        "falling_slope" = max(r_subset$falling_slope),
        "area" = max(r_subset$area),
        "rep" = max(r_subset$rep),
        "mean_i" = max(r_subset$mean_i),
        "hdd" = max(r_subset$hdd),
        "duration_hob" = max(r_subset$duration_hob),
        "peak_date" = max(r_subset$peak_date),
        "category" = max(r_subset$category),
        "imax_hob" = max(r_subset$imax_hob),
        "p_moderate" = max(r_subset$p_moderate),
        "p_strong" = max(r_subset$p_strong),
        "p_severe" = max(r_subset$p_severe),
        "p_extreme" = max(r_subset$p_extreme),
        "species" = unique_simsp[yy],
        "alpha" = 12.75
      )
      
      return(final_surv_crosssim_species_n)
    }
    
    final_surv_crosssim_species[[yy]] <- results
  } # loop for each species
})

parallel::stopCluster(cl)
sys.tm





saveRDS(final_surv_crosssim_species,"mhwsim/data/final_surv_crosssim_species_try.RDs")
#25772



