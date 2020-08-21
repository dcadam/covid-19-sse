#import libraries
library(tidyverse)
library(sandwich)
library(rstatix)

#import deidentified data with demographics and transmission settng (excludes travel)
setting_pairs <- read_csv("data/settings_pairs.csv") ##includes all cases by transmission setting i.e. infectors and infectees 
secondary_cases <- read_csv("data/secondary_cases.csv") #only infectors, previously used to caclaute n secondary.cases

#adjust variable types
setting_pairs <- setting_pairs %>%
  mutate(cluster.risk = as.factor(cluster.risk),
         male = as.factor(male))

secondary_cases <- secondary_cases %>%
  mutate(cluster.risk = as.factor(cluster.risk),
         male = as.factor(male))

#compare age and sex of all cases by setting (excluding travel)
kruskal_test(age ~ cluster.risk, data = setting_pairs)
wilcox_test(age ~ cluster.risk, data = setting_pairs)
chisq_test(x = setting_pairs$male, y = setting_pairs$cluster.risk)

setting_pairs %>%
  filter(!is.na(age)) %>%
  group_by(cluster.risk) %>%
  summarise(median = median(age))

#compare age and sex of infectors only
kruskal_test(age ~ cluster.risk, data = secondary_cases)
wilcox_test(age ~ cluster.risk, data = secondary_cases)
chisq_test(x = secondary_cases$male, y = secondary_cases$cluster.risk)

#three models for the number of secondary cases by transmission setting (cluster.risk) including age and sex
m1 <- glm.nb(secondary.cases ~ cluster.risk + age + male, data = secondary_cases)
m2 <- glm.nb(secondary.cases ~ cluster.risk + age, data = secondary_cases)
m3 <- glm.nb(secondary.cases ~ cluster.risk, data = secondary_cases)

summary(m1)
summary(m2)
summary(m3)

#Obtain robust standard errors for m2 and m3 (with and without age) becasue social variance > mean but not work or family. 
cov.m2 <- vcovHC(m2, type="HC0")
r.est.m2 <- cbind(Estimate= coef(m2),
               LL = coef(m2) - 1.96 * std.err,
               UL = coef(m2) + 1.96 * std.err)

cov.m3 <- vcovHC(m3, type="HC0")
r.est.m3 <- cbind(Estimate= coef(m3), 
                  LL = coef(m3) - 1.96 * std.err,
                  UL = coef(m3) + 1.96 * std.err)

## exponentiate estimates and bounds
rexp.est.s2 <- exp(r.est.m2)
rexp.est.s3 <- exp(r.est.m3)

#Final results with P value summary
summary(m2)
rexp.est.s2

summary(m3)
rexp.est.s3








