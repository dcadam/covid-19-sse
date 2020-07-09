library(tidyverse)
library(fitdistrplus)

#######SERIAL INTERVAL AND OFFSPRING DISTRIBUTION ANALYSIS

transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

###SERIAL INTERVAL ANALYSIS

#fit normal distribution to serial intervals
nfit <- transmission_pairs %>%
  filter(onset.diff != 'NA') %>%
  pull(onset.diff) %>%
  fitdist(data = ., distr = 'norm')

#Plot serial interval with lognormal distributions
ggplot(data = transmission_pairs) +
  geom_histogram(aes(x = onset.diff, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1) + stat_function(fun = dnorm, args = list(mean = nfit$estimate[[1]], sd = nfit$estimate[[2]]), size = 0.8, linetype = 2) +
  scale_x_continuous("Serial Interval (Days)", limits = c(-10,30), breaks = seq(-10, 30, by =5), expand = c(0,0)) +
  scale_y_continuous("Proportion", expand = c(0,0), limits = c(0,0.20)) +
  theme_classic() +
  theme(aspect.ratio = 1)


transmission_pairs %>%
  dplyr::select(onset.diff) %>%
  arrange(onset.diff) %>%
  filter(!is.na(onset.diff))

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

summary(nfit)
summary(lgfit)
summary(gfit)
summary(wfit)

#bootstrapped analysis
nfit_boot <- summary(bootdist(nfit))
gfit_boot <- summary(bootdist(gfit))
wfit_boot <- summary(bootdist(wfit))
lgfit_boot <- summary(bootdist(lgfit))

#plot serial interval with normal, gamma and weibull intervals
ggplot(data = transmission_pairs) +
  geom_histogram(aes(x = onset.diff, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1) +
  stat_function(fun = dgamma, args = list(shape = gfit$estimate[[1]], rate = gfit$estimate[[2]]), size = 0.8, linetype = 1) +
  stat_function(fun = dweibull, args = list(shape = wfit$estimate[[1]], scale = wfit$estimate[[2]]), size = 0.8, linetype = 2) +
  stat_function(fun = dlnorm, args = list(meanlog = lgfit$estimate[[1]], sdlog = lgfit$estimate[[2]]), size = 0.8, linetype = 3) +
  scale_x_continuous("Serial Interval (Days)", limits = c(-5,30), breaks = seq(-5, 30, by =5), expand = c(0,0)) +
  scale_y_continuous("Proportion", expand = c(0,0), limits = c(0,0.20)) +
  theme_classic() +
  theme(aspect.ratio = 1)


####OFFSPRING DISTRIBUTION ANALYSIS
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
  nrow() + 46 #46 Sporadic Local cases without links additional transmission

#create vector of complete offspring distribution with terminal cases having zero secondary cases
complete_offspringd <- enframe(c(offspring$n, rep(0,nterminal_infectees)))

#fit negative binomial distribution to the final offspring distribution
nbfit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

summary(nbfit)
#bootstrap analysis
nbfit_boot <- summary(bootdist(nbfit))

#plot offspring distribution with negative binomial parameters
ggplot() +
  geom_histogram(aes(x=complete_offspringd$value, y = ..density..), fill = "#dedede", colour = "Black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit$estimate[[1]], mu = nbfit$estimate[[2]])), size = 1.5) +
  stat_smooth(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit$estimate[[1]], mu = nbfit$estimate[[2]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, size = 0.5, colour = 'black') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Proportion", limits = c(0,0.7), expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)

#fit geometric distriubtion by maximum likelihood
gefit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'geom')

pfit <- complete_offspringd %>%
  pull(value) %>%
  fitdist(., distr = 'pois')

summary(gefit)
summary(pfit)

#bootsprapped analysis
gefit_boot <- summary(bootdist(gefit))
pfit_boot <- summary(bootdist(pfit))

#plot offspring distribution with geometric and Poisson distribution
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

#Negative Binomial bootstrap, plot parameters
#Plot joint estiamte of R and k
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.3)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,0.8), breaks = seq(0,0.8, by = 0.2)) ->PA

#Plot historgram of k estiamtes
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,150)) ->PC

#Plot historgram of R estiamtes
nbfit_boot$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0)) ->PB




###PAIR STATISTICS and Supplementary Figures

age_transmission_pairs <- read_csv(file = "data/age_pairs.csv")

age_transmission_pairs %>%
  mutate(transmission = as_factor(transmission)) %>%
  t.test(age ~ transmission, data = .)

###Plot age figures

ggplot(data=age_transmission_pairs) +
  geom_bar(aes(x = agegroup, fill = transmission), color = "black", position = position_dodge2(preserve = "single")) +
  theme_classic() +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Frequency", expand = c(0,0)) +
  labs(fill = "Pair") +
  theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_fill_grey()


#linear regression of by age of infector/infectee
summary(lm(agegroup.infector ~ agegroup.infectee, data = transmission_pairs))

#plot age-age matrix
ggplot(data = transmission_pairs) +
  geom_count(aes(x = agegroup.infector, y = agegroup.infectee), colour = "black", alpha = 0.4) +
  geom_smooth(method = lm, aes(x=agegroup.infector, y=agegroup.infectee), se = T, size = 0.5, colour = "black", alpha = 0.2) +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Infectee Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  theme_classic() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90), legend.position = "none")


### WAVE ONE
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

#bootstrap analysis for wave one
nbfit_boot_w1 <- summary(bootdist(nbfit_w1))

#plot offspirng disirtibution for wave one with negative binomonal parameters
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

#plot offspirng disirtibution for wave two with negative binomonal parameters
ggplot() +
  geom_histogram(aes(x=complete_offspringd_w2$value, y = ..density..), fill = "#dedede", colour = "Black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w2$estimate[[1]], mu = nbfit_w1$estimate[[2]])), size = 1.5) +
  stat_smooth(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit_w2$estimate[[1]], mu = nbfit_w1$estimate[[2]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, size = 0.5, colour = 'black') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Wave 2 Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Proportion",limits = c(0,0.8),expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)


#Plot joint estiamte of R and k wave 1
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4),breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2))

#Plot historgram of k estiamtes
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4), breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,400))

#Plot historgram of R estiamtes
nbfit_boot_w1$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,300))

#Plot joint estiamte of R and k wave 2
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +  
  geom_point(aes(x = size, y = mu), shape = 21, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4),breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2))

#Plot historgram of k estiamtes
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = size), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(expression(paste(italic("k"))), expand = c(0,0),limits = c(0,2.4), breaks = seq(0,2.4, by = 0.4)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,400))

#Plot historgram of R estiamtes
nbfit_boot_w2$estim %>%
  as_tibble() %>%
  ggplot() +
  geom_histogram(aes(x = mu), bins = 25, fill = "#dedede", colour = "black") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_x_continuous("R", expand = c(0,0),limits = c(0,1.2), breaks = seq(0,1.2, by = 0.2)) +
  scale_y_continuous("Sampling distribution", expand = c(0,0), limits = c(0,300)) 

#Supplementary Figure 5

transmission_pairs %>%
  dplyr::select(infector.case, infector.epi.date) %>%
  distinct() %>%
  left_join(., offspring, by = "infector.case") %>%
  arrange(desc(n)) %>%
  mutate(infector.epi.date = dmy(infector.epi.date)) %>%
  ggplot() +
  geom_bar(aes(x = infector.epi.date, y = n), stat = 'identity', fill = "#dedede", color = "black", size = 0.3) +
  scale_x_date("Index Illness Onset Date",
               date_breaks = "3 days", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Secondary Cases / Day", expand = c(0,0), limits = c(0,20)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, axis.text.x = element_text(angle = 45, hjust = 1))

###NEGATIVE BINOMIAL MODEL OF TRANSMISSION PAIRS EXPOSURES
library(MASS)

#Compare family and social setting exposures. 
transmission_pairs %>%
  filter(cluster.risk != "travel",
         cluster.risk != "work") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  glm.nb(.$n ~ .$cluster.risk, data = .) %>%
  summary()


#Compare work and social setting exposures. 
transmission_pairs %>%
  filter(cluster.risk != "travel",
         cluster.risk != "family") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  glm.nb(.$n ~ .$cluster.risk, data = .) %>%
  summary()

#compare work and family
transmission_pairs %>%
  filter(cluster.risk != "travel",
         cluster.risk != "social") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  glm.nb(.$n ~ .$cluster.risk, data = .) %>%
  summary()

#calculate means for resolved pairs
transmission_pairs %>%
  filter(cluster.risk != "travel") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  ungroup() %>%
  group_by(cluster.risk) %>%
  summarise(n = mean(n))


transmission_pairs %>%
  group_by(infector.case) %>%
  count() %>%
  filter(n >= 2) %>%
  arrange(desc(n))


