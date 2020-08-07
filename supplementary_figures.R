library(tidyverse)
library(lubridate)
library(fitdistrplus)
library(viridis)

#import data
case_data <- read_csv(file = "data/case_data.csv")
age_transmission_pairs <- read_csv(file = "data/age_pairs.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")
bar_data <- read_csv("data/bar_data.csv")

#recode case_data dates
case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

#Supplementary Figure 1 excluding N=106 outlier edited in post
case_data %>%
  filter(cluster.id != 0) %>%
  group_by(cluster.id) %>%
  mutate(cluster.size = n()) %>%
  arrange(epi.date) %>%
  slice(1L) %>%
  ggplot() +
  geom_jitter(aes(x = epi.date, y = cluster.size), alpha = 0.8, size = 2, shape = 21, width = 0, height = 0.8) +
  geom_smooth(method = lm, aes(x = epi.date, y = cluster.size), se = T, size = 0.5, colour = "black", alpha = 0.2) + 
  scale_x_date("Cluster index date by onset", date_breaks = "1 week", date_labels = "%d %b", minor_breaks = NULL) +
  scale_y_continuous("Final cluster size", expand = c(0,0), limits = c(0,25)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = "NULL", axis.text.x = element_text(angle = 45, hjust = 1))


##Stats for Supplementary Figure 1
case_data %>%
  filter(cluster.id != 0) %>%
  group_by(cluster.id) %>%
  mutate(cluster.size = n()) %>%
  arrange(epi.date) %>%
  slice(1L) %>%
  lm(as.numeric(epi.date) ~ cluster.size, data = .) %>%
  summary()

#Supplementary Table 1 

#Descriptive Cluster analysis
case_data %>%
  group_by(cluster.category) %>%
  summarise(n = n(), pct = n/1038)

#find overseas cluster values
case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster of imported cases") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()


#find local cluster values
case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster initiated by local case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()

#find imported cluster values
case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  n_groups()

case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  median()

case_data %>%
  filter(cluster.category == "Cluster initiated by imported case") %>%
  group_by(cluster.id) %>%
  group_size() %>%
  range()



#asymptomatic contigency table (Supplementary Table 2)
case_data %>%
  group_by(cluster.category, symptomatic) %>%
  summarise(n = n()) %>%
  spread(key = symptomatic, value = n)

#Supplementary Figure 2
ggplot(data = case_data) +
  geom_histogram(aes(x = epi.date, fill = fct_rev(symptomatic)), colour = 'black', binwidth = 1, alpha = 0.9) +
  labs(x="Illness Onset Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_x_date(date_breaks = "1 week", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Daily SARS-CoV-2 Infections in Hong Kong (N)", expand = c(0,0), limits = c(0,50)) +
  scale_fill_grey(start = 0.9, end = 0.6)

#####Supplementary Figure 3
#Figure 3A
ggplot(data=age_transmission_pairs) +
  geom_bar(aes(x = agegroup, fill = transmission), color = "black", position = position_dodge2(preserve = "single")) +
  theme_classic() +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Frequency", expand = c(0,0)) +
  labs(fill = "Pair") +
  theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_fill_grey()

#Figure 3B (missing age for 1 case responsible for 11 secondary cases)
ggplot(data = transmission_pairs) +
  geom_count(aes(x = agegroup.infector, y = agegroup.infectee), colour = "black", alpha = 0.4) +
  geom_smooth(method = lm, aes(x=agegroup.infector, y=agegroup.infectee), se = T, size = 0.5, colour = "black", alpha = 0.2) +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Infectee Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  theme_classic() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90), legend.position = "none")


#Significance tests for Figure 3A and 3B
age_transmission_pairs %>%
  mutate(transmission = as_factor(transmission)) %>%
  t.test(age ~ transmission, data = .)

summary(lm(agegroup.infector ~ agegroup.infectee, data = transmission_pairs))


###Supplementary figure 4 and supplementary table 4 (wave 1 wave 2)
## WAVE ONE
#count number of offspring per individual infector for wave one before march first
offspring_w1 <- transmission_pairs %>%
  dplyr::select(infector.case, infector.epi.date) %>%
  filter(infector.epi.date < "2020-03-01") %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee_w1 <- transmission_pairs %>%
  filter(infectee.epi.date < "2020-03-01") %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector_w1 <- transmission_pairs %>%
  filter(infector.epi.date < "2020-03-01") %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

duplicate_w1 <- infector_w1 %>%
  left_join(., infectee_w1, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct()

nterminal_infectees_w1 <- infectee_w1 %>% 
  dplyr::select(value) %>%
  filter(!value %in% duplicate_w1$value) %>%
  transmute(case.no = as.numeric(value)) %>%
  nrow() + 11

#create vector of complete offspring distribution for wave one with terminal cases having zero secondary cases
complete_offspringd_w1 <- enframe(c(offspring_w1$n, rep(0,nterminal_infectees_w1)))

#fit negative binomial distribution to wave one offspring distribution
nbfit_w1 <- complete_offspringd_w1 %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

#bootstrap analysis for wave one (Supplementary Table 4)
nbfit_boot_w1 <- summary(bootdist(nbfit_w1)) 

#calculation of proportion of cases who do not spread to anyone from nbfit and nbfit_boot of wave 1
dnbinom(0, size = 0.3981823, mu = 0.6064477)
dnbinom(0, size = 0.2071336, mu = 0.3838060)
dnbinom(0, size = 0.9461882, mu = 0.8586690)

#Plot Supplementray Figure 4A
ggplot() +
  geom_histogram(aes(x=complete_offspringd_w1$value, y = ..density..), fill = "#dedede", colour = "Black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w1$estimate[[1]], mu = nbfit_w1$estimate[[2]])), size = 1.5) +
  stat_smooth(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w1$estimate[[1]], mu = nbfit_w1$estimate[[2]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, size = 0.5, colour = 'black') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Wave 1 Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Proportion",limits = c(0,0.8),expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)


### WAVE TWO
#count number of offspring per individual infector for wave two after march first. 
offspring_w2 <-transmission_pairs %>%
  dplyr::select(infector.case, infector.epi.date) %>%
  filter(infector.epi.date >= "2020-03-01") %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee_w2 <- transmission_pairs %>%
  filter(infectee.epi.date >= "2020-03-01") %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector_w2 <- transmission_pairs %>%
  filter(infector.epi.date >= "2020-03-01") %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

duplicate_w2 <- infector_w2 %>%
  left_join(., infectee_w2, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct()

nterminal_infectees_w2 <- infectee_w2 %>% 
  dplyr::select(value) %>%
  filter(!value %in% duplicate_w2$value) %>%
  transmute(case.no = as.numeric(value)) %>%
  nrow() + 35

#create vector of complete offspring distribution for wave one with terminal cases having zero secondary cases
complete_offspringd_w2 <- enframe(c(offspring_w2$n, rep(0,nterminal_infectees_w2)))

#fit negative binomial distribution to wave one offspring distribution
nbfit_w2 <- complete_offspringd_w2 %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

#bootstrap analysis for wave two
nbfit_boot_w2 <- summary(bootdist(nbfit_w2))

#calculation of proportion of cases who do not spread to anyone from nbfit and nbfit_boot of wave 2
dnbinom(0, size = 0.4425313, mu = 0.5705630)
dnbinom(0, size = 0.2731069, mu = 0.4240127)
dnbinom(0, size = 0.8230628, mu = 0.7329678)

#Plot Supplementary Figure 4B
ggplot() +
  geom_histogram(aes(x=complete_offspringd_w2$value, y = ..density..), fill = "#dedede", colour = "Black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w2$estimate[[1]], mu = nbfit_w1$estimate[[2]])), size = 1.5) +
  stat_smooth(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w2$estimate[[1]], mu = nbfit_w1$estimate[[2]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, size = 0.5, colour = 'black') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Wave 2 Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Proportion",limits = c(0,0.8),expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)

#SUPPLEMENTARY FIGURE 5
offspring <- transmission_pairs %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

transmission_pairs %>%
  dplyr::select(infector.case, infector.epi.date) %>%
  distinct() %>%
  left_join(., offspring, by = "infector.case") %>%
  arrange(desc(n)) %>%
  mutate(infector.epi.date = dmy(infector.epi.date)) %>%
  ggplot() +
  geom_bar(aes(x = infector.epi.date, y = n), stat = 'identity', fill = "#dedede", color = "black", size = 0.3) +
  scale_x_date("Index Illness Onset Date",
               date_breaks = "1 week", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Secondary Cases / Day", expand = c(0,0), limits = c(0,20)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, axis.text.x = element_text(angle = 45, hjust = 1))


####SUPPLEMENTARY FIGURE 6 & SUPPLEMENTARY TABLE 3 & 5
###Supplementary Figure 6A and Table 3
#fit lognormal distribution by maximum likeihood
lgfit <- transmission_pairs %>%
  filter(onset.diff != 'NA', 
         onset.diff >=1 ) %>%
  pull(onset.diff) %>%
  fitdist(data = ., distr = 'lnorm')

#fit gamma distribution by maximum likeihood
gfit <- transmission_pairs %>%
  filter(onset.diff != 'NA', 
         onset.diff >=1 ) %>%
  pull(onset.diff) %>%
  fitdist(data = ., distr = 'gamma')

#fit weibull distribution by maximum likeihood
wfit <- transmission_pairs %>%
  filter(onset.diff != 'NA', 
         onset.diff >=1 ) %>%
  pull(onset.diff) %>%
  fitdist(data = ., distr = 'weibull')

summary(lgfit)
summary(gfit)
summary(wfit)

#bootstrapped analysis
gfit_boot <- summary(bootdist(gfit))
wfit_boot <- summary(bootdist(wfit))
lgfit_boot <- summary(bootdist(lgfit))

#Supplementary Figure 6A
transmission_pairs %>%
  filter(!is.na(onset.diff)) %>%
  ggplot() +
  geom_histogram(aes(x = onset.diff, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1) +
  stat_function(fun = dgamma, args = list(shape = gfit$estimate[[1]], rate = gfit$estimate[[2]]), size = 0.8, linetype = 1) +
  stat_function(fun = dweibull, args = list(shape = wfit$estimate[[1]], scale = wfit$estimate[[2]]), size = 0.8, linetype = 2) +
  stat_function(fun = dlnorm, args = list(meanlog = lgfit$estimate[[1]], sdlog = lgfit$estimate[[2]]), size = 0.8, linetype = 3) +
  scale_x_continuous("Serial Interval (Days)") +
  scale_y_continuous("Proportion", expand = c(0,0), limits = c(0,0.20)) +
  theme_classic() +
  theme(aspect.ratio = 1)

###Supplementary Figure 6B and Table 5
#count number of offspring per individual infector
offspring <- transmission_pairs %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee <- transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  transmission_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

duplicate <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct()

nterminal_infectees <- infectee %>% 
  dplyr::select(value) %>%
  filter(!value %in% duplicate$value) %>%
  transmute(case.no = as.numeric(value)) %>%
  nrow() + 46 #46 Sporadic Local cases without links to additional transmission

#create vector of complete offspring distribution with terminal cases having zero secondary cases
complete_offspringd <- enframe(c(offspring$n, rep(0,nterminal_infectees)))

#fit  distriubtions by maximum likelihood

nbfit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

gefit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'geom')

pfit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'pois')

summary(nbfit)
summary(gefit)
summary(pfit)

#bootsprapped analysis
nbfit_boot <- summary(bootdist(nbfit))
gefit_boot <- summary(bootdist(gefit))
pfit_boot <- summary(bootdist(pfit))

#Figure 6B
ggplot() +
  geom_histogram(aes(x=complete_offspringd$value, y = ..density..), fill = "#dedede", colour = "black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dgeom(x = 0:11, prob = gefit$estimate[[1]])), color = 'black', size = 2, shape = 'triangle') +
  stat_smooth(aes(x = 0:11, y = dgeom(x = 0:11, prob = gefit$estimate[[1]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, colour = 'black', size =0.8) +
  geom_point(aes(x = 0:11, y = dpois(x = 0:11, lambda = pfit$estimate[[1]])), color = 'black', size = 2, shape = 'square') +
  stat_smooth(aes(x = 0:11, y = dpois(x = 0:11, lambda = pfit$estimate[[1]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, colour = 'black', size = 0.8, linetype = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Density", limits = c(0,0.7), expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)

####Supplementary Figure 7 (dependant on negative binomial analysis above)

#Plot joint estiamte of R and k (Figure 7A)
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.3)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,1.0), breaks = seq(0,1.0, by = 0.2))

#Plot historgram of R estiamtes (Figure 7B)
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,150))

#Plot historgram of k estiamtes (Figure 7C)
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,150))


#Supplementary Figure 8 (dependant on analysis from Supplementary figure 4 above wave 1 and wave 2)
#Plot joint estiamte of R and k wave 1 (Figure 8A)
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4),breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2))

#Plot historgram of R estiamtes (Figure 8B)
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,300))

#Plot historgram of k estiamtes (Figure 8C)
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4), breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,400))

#Plot joint estiamte of R and k wave 2 (Figure 8D)
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4),breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2))

#Plot historgram of R estiamtes (Figure 8E)
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,300)) 

#Plot historgram of k estiamtes (Figure 8F)
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4), breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,400))


###SUPPLEMENTARY FIGURE 9


##Figure 9A
bar_data %>%
  mutate(cluster.generation = as.ordered(cluster.generation)) %>%
  ggplot() +
  geom_bar(aes(x = epi.date, fill = cluster.generation), width = 0.9) +
  scale_x_date(name = "Onset Date",
               date_breaks = "2 days", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Case Count", expand = c(0,0), breaks = seq(0,16, by = 2), limits = c(0,16)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1 )) +
  scale_fill_viridis_d()

#Figure 9B
bar_data %>%
  ggplot() +
  geom_bar(aes(x = age.group, fill = cluster.generation)) +
  scale_x_continuous(name = "Age Group", expand = c(0,0), breaks = 1:18,
                     labels = c("0-4",
                                "5-9",
                                "10-14",
                                "15-19",
                                "20-24",
                                "25-29",
                                "30-34",
                                "35-39",
                                "40-44",
                                "45-49",
                                "50-54",
                                "55-59",
                                "60-64",
                                "65-69",
                                "70-74",
                                "75-79",
                                "80-84",
                                "85+")) +
  scale_y_continuous("Case Count", expand = c(0,0), breaks = seq(0,25, by = 5), limits = c(0,25)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, 
        legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1 )) +
  scale_fill_viridis_d()

#Figure 9C (additional editing for production done in post)
bar_data %>%
  mutate(cluster.generation = as.ordered(cluster.generation)) %>%
  dplyr::select(-epi.date, -age.group) %>%
  filter(!is.na(onset.date)) %>%
  gather(date.type, date, -c(case.no, cluster.generation)) %>%
  mutate(case.no = factor(case.no, levels=rev(unique(case.no)), ordered=TRUE)) %>%
  mutate(date = as.Date(date, origin = "1970-01-01")) %>%
  ggplot() +
  geom_line(aes(x=rev(case.no), y=date, colour = cluster.generation), size = 3) +
  labs(x=NULL, y="Onset-to-Confirmation") +
  scale_y_date(date_breaks = "3 day", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_x_discrete(breaks = NULL) +
  theme_classic() +
  theme(legend.position = 'none', panel.grid.major.y = element_line(colour = 'grey')) +
  scale_color_viridis_d()

##Results for Supplementary Tables 6 & 7 in respective files hypothetical analysis and quarantine.

