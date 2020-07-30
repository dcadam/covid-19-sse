library(tidyverse)
library(fitdistrplus)


####HYPOTHETICAL OFFSPRING DISTRIBUTION ANALYSIS
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


###NEGATIVE BINOMIAL MODEL OF HYPOTHETICAL TRANSMISSION PAIRS by EXPOSURES
#calculate means number of secondary cases for additional scenario one pairs
hypothetical_pairs %>%
  filter(cluster.risk != "travel") %>%
  filter(cluster.id != 77.2) %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  ungroup() %>%
  group_by(cluster.risk) %>%
  summarise(n = mean(n))

##negative binomial analsysi for scenario one social vs family
hypothetical_pairs %>%
  filter(cluster.risk != "travel",
         cluster.risk != "work") %>%
  filter(cluster.id != 77.2) %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  glm.nb(.$n ~ .$cluster.risk, data = .) %>%
  summary()

#calculate means number of secondary cases for additional scenario two pairs
hypothetical_pairs %>%
  filter(cluster.risk != "travel") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  ungroup() %>%
  group_by(cluster.risk) %>%
  summarise(n = mean(n))

##negative binomial analsysi for scenario one social vs family
hypothetical_pairs %>%
  filter(cluster.risk != "travel",
         cluster.risk != "work") %>%
  group_by(infector.case, cluster.risk) %>%
  count() %>%
  glm.nb(.$n ~ .$cluster.risk, data = .) %>%
  summary()




