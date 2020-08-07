library(tidyverse)
library(lubridate)
library(ggsci)

###delay analysis and figures

case_data <- read_csv("data/case_data.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

#Figure 4A - figures combined and colours customised in post
case_data %>%
  filter(cluster.id != 0) %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  filter(!is.na(delay)) %>%
  ggplot() +
  geom_histogram(aes(x = delay, y = ..density..),  fill = '#dedede', colour = "black", binwidth = 1) +
  scale_y_continuous("Frequency", expand = c(0,0), limits = c(0,0.16)) + 
  scale_x_continuous("Delay from onset-to-confirmation (days)", expand = c(0,0), breaks = seq(0,30, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = 'none')


case_data %>%
  filter(cluster.id != 0) %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  group_by(cluster.id) %>%
  mutate(n = n()) %>% 
  dplyr::select(n, delay) %>%
  ungroup() %>%
  mutate(n = as_factor(n)) %>% 
  ggplot() +
  geom_boxplot(aes(x=n, y=delay), fill = "#dedede") +
  theme_classic() +
  theme(aspect.ratio=1) +
  scale_y_continuous("Delay from onset-to-confirmation (days)", limits = c(0, 30), expand = c(0,0), breaks = seq(0,30, by = 3)) +
  xlab("Cluster size (n)") +
  coord_flip()

#Figure 4B
transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(!is.na(delay)) %>%
  ggplot() +
  geom_histogram(aes(x = delay, y = ..density..),  fill = '#dedede', colour = "black", binwidth = 1) +
  scale_y_continuous("Frequency", expand = c(0,0), limits = c(0,0.16)) + 
  scale_x_continuous("Delay from onset-to-isolation of index (days)", 
                     expand = c(0,0), limits = c(0,27), breaks = seq(0,27, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = 'none')

transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(delay != "NA") %>%
  arrange(desc(n)) %>%
  ggplot() +
  geom_smooth(method = lm, aes(x=delay, y = n), color = "black", alpha = 0.1, size = 0.7) +
  geom_jitter(aes(x = delay, y = n, colour = cluster.risk), height = 0.3, width = 0.3) +
  scale_y_continuous("Secondary Cases / Index", breaks = 1:11) +
  scale_x_continuous("Delay from onset-to-confirmation of index (days)", 
                     expand = c(0,0),                     
                     limits = c(0,27), breaks = seq(0,27, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = 'none') +
  scale_color_lancet()

#Calculatate delay stats
transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(delay != "NA") %>% 
  arrange(desc(n)) %>%
  lm(n ~ delay, data = .) %>%
  summary()

transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(!is.na(delay)) %>%
  pull(delay) %>%
  median()

case_data %>%
  filter(cluster.id != 0) %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  filter(!is.na(delay)) %>%
  pull(delay) %>%
  median()



