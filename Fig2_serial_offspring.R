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

#fit results
summary(nfit)

#bootstrapped analysis
nfit_boot <- summary(bootdist(nfit))

#Plot serial interval with normal distribution (Figure 2A)
ggplot(data = transmission_pairs) +
  geom_histogram(aes(x = onset.diff, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1) + 
  stat_function(fun = dnorm, args = list(mean = nfit$estimate[[1]], sd = nfit$estimate[[2]]), size = 0.8, linetype = 2) +
  scale_x_continuous("Serial Interval (Days)", limits = c(-10,30), breaks = seq(-10, 30, by =5), expand = c(0,0)) +
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

#calculation of proportion of cases who do not spread to anyone from nbfit and nbfit_boot
dnbinom(0, size = 0.4258355, mu = 0.5828315)
dnbinom(0, size = 0.2869252, mu = 0.6694716)
dnbinom(0, size = 0.4520273, mu = 0.7176502)

#plot offspring distribution with negative binomial parameters (Figure 2B)
ggplot() +
  geom_histogram(aes(x=complete_offspringd$value, y = ..density..), fill = "#dedede", colour = "Black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit$estimate[[1]], mu = nbfit$estimate[[2]])), size = 1.5) +
  stat_smooth(aes(x = 0:11, y = dnbinom(x = 0:11, size = nbfit$estimate[[1]], mu = nbfit$estimate[[2]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, size = 0.5, colour = 'black') +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Secondary Cases / Index", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Proportion", limits = c(0,0.7), expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)


###NEGATIVE BINOMIAL REGRESSION ANALYSISI OF TRANSMISSION PAIRS BY EXPOSURE SETTING
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




