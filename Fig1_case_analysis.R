library(tidyverse)
library(lubridate)
library(ggsci)
library(forcats)

#import data
case_data <- read_csv(file = "data/case_data.csv")

case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

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

#Plot Main Figure 1
ggplot(data = case_data) +
  geom_bar(aes(x = epi.date, fill = cluster.category), colour = 'black', size = 0.05, binwidth = 1, alpha = 0.9) +
  labs(x="Illness Onset Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'right') +
  scale_x_date(date_breaks = "2 days", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Daily SARS-CoV-2 Cases in Hong Kong (N)", expand = c(0,0), limits = c(0,50)) +
  scale_fill_brewer(type = "div", palette = 4)

#Supplementary Figure 1 excluding N=106 outlier
case_data %>%
  filter(cluster.id != 0) %>%
  group_by(cluster.id) %>%
  mutate(cluster.size = n()) %>%
  arrange(epi.date) %>%
  slice(1L) %>%
  ggplot() +
  geom_jitter(aes(x = epi.date, y = cluster.size), alpha = 0.8, size = 2, shape = 21) +
  geom_smooth(method = lm, aes(x = epi.date, y = cluster.size), se = T, size = 0.5, colour = "black", alpha = 0.2) + 
  scale_x_date("Cluster index date by onset", date_breaks = "2 days", date_labels = "%d %b", minor_breaks = NULL) +
  scale_y_continuous("Final cluster size", expand = c(0,0), limits = c(0,25)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = "NULL", axis.text.x = element_text(angle = 45, hjust = 1))


#Supplementary Figure 2
ggplot(data = case_data) +
  geom_histogram(aes(x = epi.date, fill = fct_rev(symptomatic)), colour = 'black', binwidth = 1, alpha = 0.9) +
  labs(x="Illness Onset Date") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_x_date(date_breaks = "2 days", 
               date_labels = "%d %b", 
               minor_breaks = NULL) +
  scale_y_continuous("Daily SARS-CoV-2 Infections in Hong Kong (N)", expand = c(0,0), limits = c(0,50)) +
  scale_fill_grey(start = 0.9, end = 0.6)

#asymptomatic contigency table
case_data %>%
  group_by(cluster.category, symptomatic) %>%
  summarise(n = n()) %>%
  spread(key = symptomatic, value = n)




