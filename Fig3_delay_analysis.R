library(tidyverse)
library(lubridate)

###delay analysis and figures

case_data <- read_csv("data/case_data.csv")
transmission_pairs <- read_csv(file = "data/transmission_pairs.csv")

case_data <- case_data %>%
  mutate(onset.date = dmy(onset.date),
         confirm.date = dmy(confirm.date),
         epi.date = dmy(epi.date))

#Figure 4A - figures combined and colours customised in post
case_data %>%
  filter(cluster.id != 0,
         cluster.category != "Cluster of imported cases") %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  filter(!is.na(delay)) %>%
  ggplot() +
  geom_histogram(aes(x = delay, y = ..density..),  fill = '#dedede', colour = "black", binwidth = 1) +
  scale_y_continuous("Frequency", expand = c(0,0), limits = c(0,0.20)) + 
  scale_x_continuous("Delay from onset-to-confirmation (days)", expand = c(0,0), breaks = seq(0,27, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = 'none')

case_data %>%
  filter(cluster.id != 0, 
         cluster.category != "Cluster of imported cases") %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  group_by(cluster.id) %>%
  mutate(n = n()) %>% 
  dplyr::select(n, delay) %>%
  ungroup() %>%
  mutate(n = as_factor(n)) %>%
  group_by(n) %>%
  filter(!is.na(delay)) %>%
  summarise(ymin = min(delay),
            ymax = max(delay),
            middle = median(delay),
            lower = quantile(delay,0.25),
            upper = quantile(delay,0.75)) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = n, ymin = ymin,ymax = ymax,middle = middle,upper = upper,lower= lower),
               stat = 'identity', 
               fill = "#dedede",
               width = 0.8) +
  theme_classic() +
  theme(aspect.ratio=1) +
  scale_y_continuous("Delay from onset-to-confirmation (days)", limits = c(0, 27), expand = c(0,0), breaks = seq(0,27, by = 3)) +
  xlab("Cluster size (n)") +
  coord_flip()

#Figure 4B
transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(!is.na(delay)) %>%
  ggplot() +
  geom_histogram(aes(x = delay, y = ..density..),  fill = '#dedede', colour = "black", binwidth = 1) +
  scale_y_continuous("Frequency", expand = c(0,0), limits = c(0,0.20)) + 
  scale_x_continuous("Delay from onset-to-isolation of infector (days)", 
                     expand = c(0,0), limits = c(0,27), breaks = seq(0,27, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 0.3, legend.position = 'none')

transmission_pairs %>%  
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(!is.na(delay)) %>%
  arrange(desc(n)) %>%
  ggplot() +
  geom_smooth(method = lm, aes(x=delay, y = n), color = "black", alpha = 0.1, size = 0.7) +
  geom_jitter(aes(x = delay, y = n, colour = cluster.risk), height = 0.3, width = 0.3) +
  scale_y_continuous("Secondary Cases / Infector", breaks = 1:11) +
  scale_x_continuous("Delay from onset-to-confirmation of infector (days)", 
                     expand = c(0,0),                     
                     limits = c(0,27), breaks = seq(0,27, by = 3)) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = c(0.85, 0.85), legend.title = element_blank())  #colours are modified custom in post 

#####Calculate delay stats
##median delay across all local clusters
case_data %>%
  filter(cluster.id != 0,
         cluster.category != "Cluster of imported cases") %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  filter(!is.na(delay)) %>%
  pull(delay) %>%
  median()

#Relationship between median delay within clusters and cluster size (excluding two largest)
case_data %>%
  filter(cluster.id != 0, 
         cluster.category != "Cluster of imported cases") %>%
  mutate(delay = as.numeric(confirm.date - onset.date, na.rm = T)) %>%
  group_by(cluster.id) %>%
  mutate(n = n()) %>% 
  dplyr::select(n, delay) %>%
  ungroup() %>%
  filter(!is.na(delay)) %>%
  group_by(n) %>%
  summarise(m.delay = median(delay)) %>%
  filter(n < 22) %>%
  lm(n ~ m.delay, data = .) %>%
  summary()

#median delay among infectors
transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = mean(delay.infector)) %>%
  filter(!is.na(delay)) %>%
  pull(delay) %>%
  median()

#relationship between delay in confirmation of infectors and number of secondary cases
transmission_pairs %>%
  group_by(infector.case, cluster.risk) %>%
  summarise(n = n(), delay = median(delay.infector)) %>%
  filter(delay != "NA") %>% 
  arrange(desc(n)) %>%
  lm(n ~ delay, data = .) %>%
  summary()






