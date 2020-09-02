library(tidyverse)
library(lubridate)
library(fitdistrplus)
library(igraph)


#import data
case_data <- read_csv(file = "data/case_data.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

#recode case_data dates
case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))


##### SUPPLEMENTARY TABLE 1
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

#SUPPLEMENTARY TABLE 2
case_data %>%
  group_by(cluster.category, symptomatic) %>%
  summarise(n = n()) %>%
  spread(key = symptomatic, value = n)


#SUPPLEMENTARY TABLE 3
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

#SUPPLEMENTARY TABLE 4
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




#SUPPLEMENTARY TABLE 5
#count number of offspring per individual infector for wave one before pre march
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

#bootstrap analysis for pre march
nbfit_boot_w1 <- summary(bootdist(nbfit_w1)) 

### Post March
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

##Coefficient of Variation and Porportion 80%
#proportion resposible formula
propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}

#Empirical Estiamtes P80% Total
propresponsible(0.58,0.43,0.8)
propresponsible(0.45,0.29,0.8)
propresponsible(0.72,0.67,0.8)

#Branching Process p80% Total
propresponsible(0.74,0.33,0.8)
propresponsible(0.58,0.14,0.8)
propresponsible(0.97,0.98,0.8)

#Hypothetical Scenario One P80% Total
propresponsible(0.62,0.35,0.8)
propresponsible(0.49,0.25,0.8)
propresponsible(0.70,0.56,0.8)

#Hypothetical Scenario Two P80% Total
propresponsible(0.72,0.19,0.8)
propresponsible(0.53,0.13,0.8)
propresponsible(0.94,0.26,0.8)

#Empirical Estiamtes P80% Wave 1
propresponsible(0.61,0.40,0.8)
propresponsible(0.38,0.21,0.8)
propresponsible(0.88,0.86,0.8)

#Empirical Estiamtes P80% Wave 2
propresponsible(0.57,0.44,0.8)
propresponsible(0.42,0.27,0.8)
propresponsible(0.73,0.82,0.8)

#Branching Process p80% Wave 1
propresponsible(0.66,2.31,0.8)
propresponsible(0.44,0.22,0.8)
propresponsible(0.99,1*10^100,0.8)

#Branching Process p80% Wave2
propresponsible(0.77,0.20,0.8)
propresponsible(0.54,0.08,0.8)
propresponsible(1.13,0.63,0.8)


#Coeffcient of variation formula
nbinomcoefficient=function(r,k){
  var=r*(1+(r/k))
  sd=sqrt(var)
  sd/r
}

#Empirical Estiamtes CV
nbinomcoefficient(0.58,0.43)
nbinomcoefficient(0.45,0.29)
nbinomcoefficient(0.72,0.67)

#Hypothetical Scenario One CV
nbinomcoefficient(0.62,0.35)
nbinomcoefficient(0.49,0.25)
nbinomcoefficient(0.70,0.56)

#Hypothetical Scenario Two CV
nbinomcoefficient(0.72,0.19)
nbinomcoefficient(0.53,0.13)
nbinomcoefficient(0.94,0.26)

#Empirical Estiamtes CV Wave 1
nbinomcoefficient(0.61,0.40)
nbinomcoefficient(0.38,0.21)
nbinomcoefficient(0.88,0.86)

#Empirical Estiamtes CV Wave 2
nbinomcoefficient(0.57,0.44)
nbinomcoefficient(0.42,0.27)
nbinomcoefficient(0.73,0.82)

#Branching Process CV
nbinomcoefficient(0.74,0.33)
nbinomcoefficient(0.58,0.14)
nbinomcoefficient(0.97,0.98)

#Branching Process CV Wave 1
nbinomcoefficient(0.66,2.31)
nbinomcoefficient(0.44,0.22)
nbinomcoefficient(0.99,Inf)

#Branching Process CV Wave2
nbinomcoefficient(0.77,0.21)
nbinomcoefficient(0.54,0.08)
nbinomcoefficient(1.13,0.63)


####HYPOTHETICAL SECOANRIOS OFFSPRING DISTRIBUTION ANALYSIS
hypothetical_pairs <- read_csv(file = "data/hypothetical_pairs.csv")

hypothetical_pairs %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

####SCENARIO ONE
#count number of offspring per individual infector filtering out unresolved bar and band cases
offspring_1 <- hypothetical_pairs %>%
  filter(cluster.id != 77.2) %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee <- hypothetical_pairs %>%
  filter(cluster.id != 77.2) %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  hypothetical_pairs %>%
  filter(cluster.id != 77.2) %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infector.case')

duplicate <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct()

nterminal_infectees_1 <- infectee %>% 
  dplyr::select(value) %>%
  filter(!value %in% duplicate$value) %>%
  transmute(case.no = as.numeric(value)) %>%
  nrow() + 46 #46 Sporadic Local cases without links additional transmission

#create vector of complete offspring distribution with terminal cases having zero secondary cases
hypothetical_offspringd_1 <- enframe(c(offspring_1$n, rep(0,nterminal_infectees_1)))

#fit negative binomial distribution to the final offspring distribution
nbfit_1 <- hypothetical_offspringd_1 %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

summary(nbfit_1)
#bootstrap analysis
nbfit_boot_1 <- summary(bootdist(nbfit_1))

#count number of offspring per individual infector
offspring_2 <- hypothetical_pairs %>%
  dplyr::select(infector.case) %>%
  group_by(infector.case) %>%
  count() %>%
  arrange(desc(n))

#count number of terminal infectees including sporadic local cases
infectee <- hypothetical_pairs %>%
  dplyr::select(infector.case, infectee.case) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  hypothetical_pairs %>%
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
hypothetical_offspringd_2 <- enframe(c(offspring_2$n, rep(0,nterminal_infectees)))

#fit negative binomial distribution to the final offspring distribution
nbfit_2 <- hypothetical_offspringd_2 %>%
  pull(value) %>%
  fitdist(., distr = 'nbinom')

summary(nbfit_2)
#bootstrap analysis
nbfit_boot_2 <- summary(bootdist(nbfit_2))


#calculation of proportion of cases who do not spread to anyone from fit 1 and fit 2
dnbinom(0, size = 0.3511851, mu = 0.6170447)
dnbinom(0, size = 0.2494476, mu = 0.4720252)
dnbinom(0, size = 0.5296645, mu = 0.7753710)

dnbinom(0, size = 0.1852538, mu = 0.7214719)
dnbinom(0, size = 0.1342745, mu = 0.5286102)
dnbinom(0, size = 0.2570040, mu = 0.9408714)


#SUPPLEMENTARY TABLE 6
quarantine_list <-  case_data %>%
  dplyr::select(case.no, quarantine) %>%
  mutate(id = as.character(case.no)) %>%
  dplyr::select(id, quarantine)

temple_cluster <- read_csv(file = "data/temple_cluster.csv")
wedding_cluster <- read_csv(file = "data/wedding_cluster.csv")
bar_cluster <- read_csv(file = "data/bar_cluster.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

complete_pairs <- transmission_pairs %>%
  filter(cluster.id != 46,
         cluster.id != 77,
         cluster.id != 80) %>%
  mutate(from = infector.case, to = infectee.case) %>%
  dplyr::select(from, to, cluster.risk, cluster.generation) %>%
  full_join(., temple_cluster) %>%
  full_join(., wedding_cluster) %>%
  full_join(., bar_cluster)


###ENCODE CASE POSITIONS WITHIN THE ENTIRE NETWORK AS SOURCE, INTERMEDIATE, OR TERMINAL
infectee <- complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  gather() %>%
  filter(key == 'infectee.case')

infector <-  complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  gather() %>%
  filter(key == 'infector.case')

middle_cases <- infector %>%
  left_join(., infectee, by = 'value') %>%
  filter(key.y != 'NA') %>%
  dplyr::select(value) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "middle")

terminal_cases <- infectee %>% 
  dplyr::select(value) %>%
  filter(!value %in% middle_cases$position) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "terminal")

sources <- infector %>% 
  dplyr::select(value) %>%
  filter(!value %in% middle_cases$case.no) %>%
  distinct() %>%
  transmute(case.no = value, 
            position = "source")

#find those whole are terminal but are linked to the bar directly and exclude them due to their ambigious position in the chain
bars_direct <- complete_pairs %>%
  transmute(infector.case = from,
            infectee.case = to) %>%
  filter(infector.case == "BARS") %>%
  dplyr::select(infectee.case)

terminal_cases_xbar <- terminal_cases %>%
  filter(!case.no %in% bars_direct$infectee.case) %>%
  distinct() %>%
  transmute(case.no = case.no, 
            position = "terminal")

#combine source, middle, and terminal adding a position coulmn, then terminal Y/N, then quarantine column 
position_list <- sources %>%
  full_join(middle_cases) %>%
  full_join(terminal_cases_xbar)

#excluding cluster sources cases, cross table and count
position_list %>%
  filter(position != "source") %>%
  mutate(terminal = case_when(position == "terminal" ~ "Y",
                              TRUE ~ "N")) %>%
  mutate(id = case.no) %>%
  dplyr::select(-case.no) %>%
  left_join(quarantine_list, by = "id") %>%
  mutate(quarantine = case_when(is.na(quarantine) ~ "N", 
                                TRUE ~ quarantine)) %>%
  group_by(terminal, quarantine) %>%
  count()

#odds ratio
(45/144)/(1/46)

