library(tidyverse)
library(lubridate)
library(fitdistrplus)
library(viridis)

#import data
case_data <- read_csv(file = "data/case_data.csv")
age_transmission_pairs <- read_csv(file = "data/age_pairs.csv")
sex_transmission_pairs <- read_csv(file = "data/sex_transmission_pairs.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")
bar_data <- read_csv("data/bar_data.csv")

#recode case_data dates
case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

#Extended Data Figure 1 excluding N=106 outlier edited in post
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


##Stats for Extended Data Figure 1
case_data %>%
  filter(cluster.id != 0) %>%
  group_by(cluster.id) %>%
  mutate(cluster.size = n()) %>%
  arrange(epi.date) %>%
  slice(1L) %>%
  lm(as.numeric(epi.date) ~ cluster.size, data = .) %>%
  summary()


#Extended Data Figure 2
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

###Extended Data Figure 3
##Figure 3A
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

#Figure 3B
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

#Figure 3C (additional editing for production done in post)
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

####EXTENDED DATA FIGURE 4 
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

#Extended Data Figure 4A
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

###Extended Data Figure 4B
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

#Figure 4B
ggplot() +
  geom_histogram(aes(x=complete_offspringd$value, y = ..density..), fill = "#dedede", colour = "black", binwidth = 1) +
  geom_point(aes(x = 0:11, y = dgeom(x = 0:11, prob = gefit$estimate[[1]])), color = 'black', size = 2, shape = 'triangle') +
  stat_smooth(aes(x = 0:11, y = dgeom(x = 0:11, prob = gefit$estimate[[1]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, colour = 'black', size =0.8) +
  geom_point(aes(x = 0:11, y = dpois(x = 0:11, lambda = pfit$estimate[[1]])), color = 'black', size = 2, shape = 'square') +
  stat_smooth(aes(x = 0:11, y = dpois(x = 0:11, lambda = pfit$estimate[[1]])), method = 'lm', formula = y ~ poly(x, 9), se = FALSE, colour = 'black', size = 0.8, linetype = 3) +
  expand_limits(x = 0, y = 0) +
  scale_x_continuous("Secondary Cases / Infector", expand = c(0, 0), breaks = 0:11)  +
  scale_y_continuous("Density", limits = c(0,0.7), expand = c(0, 0)) +
  theme_classic() +
  theme(aspect.ratio = 1)

#####Extended Data Figure 5
#Figure 5A
ggplot(data=age_transmission_pairs) +
  geom_bar(aes(x = agegroup, fill = transmission), color = "black", position = position_dodge2(preserve = "single")) +
  theme_classic() +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Frequency", expand = c(0,0)) +
  labs(fill = "Pair") +
  theme(aspect.ratio = 1, legend.position = c(0.95, 0.95), legend.title = element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_fill_grey()

#Figure 5B (missing age for 1 case responsible for 11 secondary cases)
ggplot(data = transmission_pairs) +
  geom_count(aes(x = agegroup.infector, y = agegroup.infectee), colour = "black", alpha = 0.4) +
  geom_smooth(method = lm, aes(x=agegroup.infector, y=agegroup.infectee), se = T, size = 0.5, colour = "black", alpha = 0.2) +
  scale_x_continuous(name = "Infector Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  scale_y_continuous(name = "Infectee Age", breaks = 1:18,
                     labels = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")) +
  theme_classic() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90), legend.position = "none")


#Significance tests for Figure 5A and 5B
age_transmission_pairs %>%
  mutate(transmission = as_factor(transmission)) %>%
  t.test(age ~ transmission, data = .)

summary(lm(agegroup.infector ~ agegroup.infectee, data = transmission_pairs))

chisq_test(x = sex_transmission_pairs$male, y = sex_transmission_pairs$transmission)