library(tidyverse)
library(ggplot2)
library(ggthemes)
library(dplyr)




# mean_ddg_score_by_grouping_new_residues


# Reading the csv file 
entireProtein <- read.csv("entireProtein.csv")

# Removing the columns with sequences in the df
proteinshort <- entireProtein[, -c(7)]
proteinshort <- proteinshort[, (-c(11))]

# Creating the new df by grouping the old amino acids (the amino acids that are about to be mutated)
data <- proteinshort %>% 
  filter(!is.na(ddg_score), !is.na(new_aa_1l)) %>% 
  group_by(previous_aa) %>% 
  summarize(mean_ddg = mean(ddg_score), sd_ddg = sd(ddg_score)) %>% 
  arrange(desc(mean_ddg))

# Calculate mean of mean_ddg column
mean_ddg_mean <- mean(data$mean_ddg)

# Sort the data by mean_ddg in ascending order
data <- data[order(desc(data$mean_ddg)),]

# Create the bar chart
p <- ggplot(data, aes(x = factor(previous_aa, levels = previous_aa), y = mean_ddg)) +
  geom_bar(stat = "identity", fill = "#F28F3B") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 800, by = 100), limits = c(0, 800)) +
  theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_blank()) +
  ggtitle("Mean ddg scores grouped by the amino acid that was added") +
  labs(x = "Previous amino acid", y = "Mean ddg") +
  geom_hline(yintercept = mean_ddg_mean, color = "red", size = 1.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1, 1, 1, 0.5), "cm"))

p 

ggsave("mean_ddg_score_by_prev_aa.png", p, width = 8, height = 6, dpi = 300)


