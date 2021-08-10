library(SOWPtrend)
library(tidyverse
        )

calc_BCI<- function(x, interval=0.90){
  c(map_estimate(x), ci(x, method="HDI", ci=interval)$CI_low,ci(x, method="HDI", ci=interval)$CI_high )}

data_name <- "SOWP.csv"
path_name <- system.file("extdata", data_name, package = "SOWPtrend", mustWork = T)
SOWP.df <- read.csv(path_name, header = T, stringsAsFactors=F)


SOWP_full.df <- create_df(SOWP.df %>% 
                            select(-Strata) %>% 
                            dplyr::rename(Strata=BCR)) %>%
  mutate(logA= log(5*5), 
         Strata=as.factor(Strata)) #clarifying that all plots are 5x5
  
#Setting the "Strata to BCR"
#Setting the "Strata to BCR"


#Now we need to add in the 2019- 2021 data (or really lack there of)

SOWP_2021 <- SOWP_full.df %>%
  complete(Year = 1971:2021, nesting(Strata,logA,  SiteID, SpeciesID), fill=list(flown=F, TIP=NA))


MALL.df <- select_sp(x=SOWP_2021, sp="MALL", year=1971:2021)


MALL.jags.data <- BWSPowerAnalysis::dataForSimulation(x=MALL.df) #USing bayesian power packages functions to prepare the data 
MALL.jags.data <- MALL.jags.data[c(1:11, 16:17)] #Because we aren't going to do the simulation, just project forward, we can remove some of the created parameters
MALL.jags.data$mask_TIP <- ifelse(is.na(MALL.jags.data$TIP), 0, 1)

str(MALL.jags.data)

my.params <-   c("mu.int","sigma.int.plot","epsilon.int.plot", "alpha", "TIP.plot",
                 "mu.beta","sigma.beta.plot","epsilon.beta.plot","beta",
                 "epsilon.future",
                 "lambda.plot", "lambda.strata", 
                 "fit_seen","fit_seen_new", "TIP.plot_new"
)
my.inits= function(w){list(
  TIP = ifelse(is.na(w$TIP), round(mean(w$TIP, na.rm=T)), NA),
  
  mu.int = rnorm(w$n.strata, 0, 2),
  sigma.int.plot =  runif(1,0.01,0.8), 
  epsilon.int.plot = rnorm(w$n.plot,0,0.1),
  # epsilon.int.plot[w$n.plot] <- NA
  
  mu.beta  =  runif(w$n.strata,-0.05,0.05),
  sigma.beta.plot = runif(1,0.01,0.4), 
  # epsilon.beta.plot[w$n.plot] <-  NA
  epsilon.beta.plot = rnorm(w$n.plot,0,0.1),
  #epsilon.beta.plot.star = rnorm(w$n.plot,0,0.1),
  #theta.TIP = runif(1,0.10,1),
  rho.TIP = runif(w$n.obs,0.90,1.1)
  
  # TIP_new = ifelse(is.na(w$TIP), 1, w$TIP)
)}

MALL.model.results <- jagsUI::jags(data= MALL.jags.data, 
                                parameters.to.save =my.params, inits = function(x){my.inits(MALL.jags.data)},
                                model.file = file.path("jags", "historical_nb_projections.txt"),
                                #n.chains = 3, n.adapt = 0, n.iter=10, n.burnin = 0, n.thin = 1, parallel=T
                                #n.chains = 3, n.adapt = 0, n.iter=100, n.burnin = 0, n.thin = 1, parallel=T
                                n.chains=6, n.adapt=1000, n.iter= 13000, n.burnin =7000, n.thin=10, parallel = T
)
#constrained the prior on the mean slopes

mean(MALL.model.results$sims.list$fit_seen >MALL.model.results$sims.list$fit_seen_new)
plot(MALL.model.results$sims.list$fit_seen ~MALL.model.results$sims.list$fit_seen_new)
abline(a=0,b=1,col="red")

save(MALL.model.results, file= file.path("output", "MALL_output_Aug 6 2021.txt"))


# #DO the same for WODU
# 
# WODU.df <- select_sp(x=SOWP_2021, sp="WODU", year=1971:2021)
# 
# 
# WODU.jags.data <- BWSPowerAnalysis::dataForSimulation(x=WODU.df) #USing bayesian power packages functions to prepare the data 
# WODU.jags.data <- WODU.jags.data[c(1:11, 16:17)] #Because we aren't going to do the simulation, just project forward, we can remove some of the created parameters
# str(WODU.jags.data)
# 
# my.params <-   c("mu.int","sigma.int.plot","epsilon.int.plot", "alpha", "TIP.plot",
#                  "mu.beta","sigma.beta.plot","epsilon.beta.plot","beta",
#                  "epsilon.future",
#                  "lambda.plot", "lambda.strata", 
#                  "fit_seen","fit_seen_new"
# )
# my.inits= function(){BWSPowerAnalysis::create_inits(WODU.jags.data)}
# 
# WODU.model.results <- jagsUI::jags(data= WODU.jags.data, 
#                                    parameters.to.save =my.params, inits = my.inits,
#                                    model.file = "C:/Users/coxa/Documents/Requested Analyses/SOWPS 2019-2021 Projections/historical_nb_projections.txt",
#                                    #n.chains = 1, n.adapt = 0, n.iter=10, n.burnin = 0, n.thin = 1, parallel=T
#                                    #n.chains = 3, n.adapt = 0, n.iter=100, n.burnin = 0, n.thin = 1, parallel=T
#                                    n.chains=3, n.adapt=1000, n.iter= 20000, n.burnin =10000, n.thin=15, parallel = T
# )
# save(WODU.model.results, file= "C:/Users/coxa/Documents/Requested Analyses/SOWPS 2019-2021 Projections/WODU_output_July 28 2021.txt")
# 
# 
# #Do the same for CAGO
# CAGO.df <- select_sp(x=SOWP_2021, sp="CAGO", year=1971:2021)
# 
# 
# CAGO.jags.data <- BWSPowerAnalysis::dataForSimulation(x=CAGO.df) #USing bayesian power packages functions to prepare the data 
# CAGO.jags.data <- CAGO.jags.data[c(1:11, 16:17)] #Because we aren't going to do the simulation, just project forward, we can remove some of the created parameters
# str(CAGO.jags.data)
# 
# my.params <-   c("mu.int","sigma.int.plot","epsilon.int.plot", "alpha", "TIP.plot",
#                  "mu.beta","sigma.beta.plot","epsilon.beta.plot","beta",
#                  "epsilon.future",
#                  "lambda.plot", "lambda.strata", 
#                  "fit_seen","fit_seen_new"
# )
# my.inits= function(){BWSPowerAnalysis::create_inits(CAGO.jags.data)}
# 
# CAGO.model.results <- jagsUI::jags(data= CAGO.jags.data, 
#                                    parameters.to.save =my.params, inits = my.inits,
#                                    model.file = "C:/Users/coxa/Documents/Requested Analyses/SOWPS 2019-2021 Projections/historical_nb_projections.txt",
#                                    #n.chains = 1, n.adapt = 0, n.iter=10, n.burnin = 0, n.thin = 1, parallel=T
#                                    #n.chains = 3, n.adapt = 0, n.iter=100, n.burnin = 0, n.thin = 1, parallel=T
#                                    n.chains=3, n.adapt=1000, n.iter= 20000, n.burnin =10000, n.thin=15, parallel = T
# )
# save(CAGO.model.results, file= "C:/Users/coxa/Documents/Requested Analyses/SOWPS 2019-2021 Projections/CAGO_output_July 28 2021.txt")





output <- MALL.model.results$sims.list

pred.df <- data.frame(Year= rep(min(SOWP.df$Year):2021, times=2), 
                      Strata = rep(c(12, 13), each= MALL.jags.data$n.year), 
                      Density = NA, 
                      LCI= NA, 
                      UCI=NA, 
                      lambda= NA, 
                      lambda_LCL= NA, 
                      lambda_UCL=NA, 
                      Density_new=NA, 
                      Density_new_LCL=NA, 
                      Density_new_UCL=NA)

#BCR 12
#Sum across all plots
output$strataDensity12 <- apply(output$TIP.plot[, which(MALL.jags.data$p.strata==1), ],  c(1,3), 
      function(x){
        sum(x)/(sum(MALL.jags.data$p.strata==1))})
pred.df[which(pred.df$Strata==12), 3:5]<- t(apply(output$strataDensity12, 2, startR::params_CI))
pred.df[which(pred.df$Strata==12), 6:8]<- t(apply(output$lambda.strata[,1,], 2, startR::params_CI))


output$strataDensity12_new <- apply(output$TIP.plot_new[, which(MALL.jags.data$p.strata==1), ],  c(1,3), 
                                function(x){
                                  sum(x)/(sum(MALL.jags.data$p.strata==1))})
pred.df[which(pred.df$Strata==12), 9:11]<- t(apply(output$strataDensity12_new, 2, startR::params_CI))


#BCR 13 

output$strataDensity13 <- apply(output$TIP.plot[, which(MALL.jags.data$p.strata==2), ],  c(1,3), 
      function(x){
        sum(x)/(sum(MALL.jags.data$p.strata==2))})

pred.df[which(pred.df$Strata==13), 3:5]<-t(apply(output$strataDensity13, 2, startR::params_CI))
pred.df[which(pred.df$Strata==13), 6:8]<- t(apply(output$lambda.strata[,2,], 2, startR::params_CI))


output$strataDensity13_new <- apply(output$TIP.plot_new[, which(MALL.jags.data$p.strata==2), ],  c(1,3), 
                                    function(x){
                                      sum(x)/(sum(MALL.jags.data$p.strata==1))})
pred.df[which(pred.df$Strata==13), 9:11]<- t(apply(output$strataDensity13_new, 2, startR::params_CI))

pred.df <- pred.df %>% full_join(SOWP_full.df  %>% group_by(Year) %>% dplyr::summarise(Flown= sum(flown)>0))

ggplot(  pred.df %>% filter(Flown==T | Year %in% c(2019:2021)), 
  aes(x=Year, fill=factor(Strata)))+
   geom_line(data= pred.df %>% filter(Flown==T | Year %in% c(2019:2021)), 
             aes(y=Density, color=factor(Strata)))+
   geom_ribbon(data=pred.df %>% filter(Flown==T | Year %in% c(2019:2021)), aes(ymin=LCI, ymax=UCI), alpha=0.3)+
  geom_line(data=pred.df , aes(y=lambda, color=factor(Strata)))+
  geom_ribbon(data=pred.df , aes(ymin=lambda_LCL, ymax=lambda_UCL), alpha=0.3)+
  
  facet_grid(~factor(Strata))+
  labs(y="TIP per plot", color="BCR", fill="BCR")+
  geom_vline(xintercept=2019, linetype="dashed")+
  ggthemes::theme_few()
ggsave(file.path("Requested Analyses","SOWPS 2019-2021 Projections", "figures", "TIP SOWP for flown years.jpeg"),
       width=8, height=4, units="in")



ggplot()+
  #geom_count(data= SOWP_2021, aes(x=Year, y=TIP, color=Strata))+ 
  geom_line(data=pred.df , aes(x=Year, y=lambda, color=factor(Strata)))+ 
  geom_ribbon(data=pred.df , aes(x=Year, ymin=lambda_LCL, ymax=lambda_UCL, color=factor(Strata)), alpha=0.3)+
  facet_grid(~Strata)+
  geom_smooth(data= SOWP_2021, aes(x=Year, y=TIP, color=Strata), method= "glm", method.args=list(family= "quasipoisson"), linetype="dashed")
  


ggplot(pred.df, aes(x=Year, fill=Strata, color=Strata))+
  geom_line(aes(y=Density))
  

ggplot()+
  geom_line(data=pred.df, aes(x=Year, y=Density, color=factor(Strata)))+
  geom_line(data=pred.df, aes(x=Year, y=lambda, color=factor(Strata)))+
  geom_ribbon(data=pred.df, aes(x=Year, ymin=lambda_LCL, ymax=lambda_UCL, fill=factor(Strata)), alpha=0.3)+
  
  
  #geom_count(data= SOWP_full.df, aes(x=Year, y=TIP))+
  geom_ribbon(data=pred.df, aes(x=Year, ymin=LCI, ymax=UCI, fill=factor(Strata)), alpha=0.3)+
  facet_grid(~factor(Strata))+
  labs(y=" TIP per plot", color="BCR", fill="BCR")+
  geom_vline(xintercept=2019, linetype="dashed")+
  ggthemes::theme_few()
ggsave(file.path("Requested Analyses","SOWPS 2019-2021 Projections", "figures", "TIP SOWP for all years.jpeg"),
       width=8, height=4, units="in")

SOWP_full.df  %>% group_by(Year) %>% dplyr::summarise(Flown= sum(flown)>0)

