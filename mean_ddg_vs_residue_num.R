library(tidyverse)
library(ggplot2)
library(ggthemes)
library(dplyr)

entireProtein <- read.csv("entireProtein.csv")

# Removing the columns with sequences in the df
proteinshort <- entireProtein[, -c(7)]
proteinshort <- proteinshort[, (-c(11))]


# Creating the new df by grouping the old amino acids (the amino acids that are about to be mutated)
data <- proteinshort %>% 
  filter(type == "mutant") %>% 
  group_by(residue_number) %>% 
  summarize(mean_ddg = mean(ddg_score), sd_ddg = sd(ddg_score))


# Calculate mean of mean_ddg column
mean_ddg_mean <- mean(data$mean_ddg)

# Sort the data by mean_ddg in ascending order
Ordereddata <- data[order(desc(data$mean_ddg)),]

# Create the bar chart
p_order <- ggplot(Ordereddata, aes(x = factor(residue_number, levels = residue_number), y = mean_ddg)) +
  geom_bar(stat = "identity", fill = "#F28F3B") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1700, by = 100), limits = c(0, 1700)) +
  theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_blank()) +
  ggtitle("Mean ddg scores for all of the 19 mutations at residue number x") +
  labs(x = "Residue Number", y = "Mean ddg") +
  geom_hline(yintercept = mean_ddg_mean, color = "red", size = 1.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1, 1, 1, 0.5), "cm"))

p_order


# image is simply too large to be produced
options(ragg.max_dim = 200000)
ggsave(filename = "charts/mean_ddg_score_per_residue_num_from_max_to_min.png", p_order, width = 100, height = 50, dpi = 300, limitsize = FALSE)


p_unorder <- ggplot(data, aes(x = factor(residue_number, levels = residue_number), y = mean_ddg)) +
  geom_bar(stat = "identity", fill = "#F28F3B") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1700, by = 100), limits = c(0, 1700)) +
  theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_blank()) +
  ggtitle("Mean ddg scores for all of the 19 mutations at residue number x") +
  labs(x = "Residue Number", y = "Mean ddg") +
  geom_hline(yintercept = mean_ddg_mean, color = "red", size = 1.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1, 1, 1, 0.5), "cm"))

p_unorder

ggsave(filename = "charts/mean_ddg_score_per_residue_num_from_1.png", p_unorder, width = 100, height = 50, dpi = 300, limitsize = FALSE)



###################
###################
# Producing the mean_ddg vs continous batches of residues



# Group data into batches of 35
batch_size <- 35
data$batch <- cut(as.numeric(as.character(data$residue_number)), 
                  breaks = seq(0, max(as.numeric(as.character(data$residue_number))), batch_size),
                  labels = FALSE)

# Calculate mean ddg for each batch
data_batched <- data %>%
  group_by(batch) %>%
  summarize(mean_ddg = mean(mean_ddg))

# Create the bar chart with batched x-axis labels
p_batch <- ggplot(data_batched, aes(x = factor(batch), y = mean_ddg)) +
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
  labs(x = "Residue Number (batches of 35)", y = "Mean ddg") +
  geom_hline(yintercept = mean(data$mean_ddg), color = "red", size = 1.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1, 1, 1, 0.5), "cm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_batch

ggsave(filename = file.path("charts", "mean_ddg_score_per_batch_of_residue_num.png") , plot = p_batch, width = 8, height = 6, dpi = 300)






