library(tidyverse)
library(ggplot2)
library(ggthemes)
library(dplyr)

# removing the seq headers for higher readability 
entireProtein <- read.csv("entireProtein.csv")
proteinshort <- entireProtein[, -c(7)]
proteinshort <- proteinshort[, (-c(11))]


# Create the data frame with the mean ddg scores, filtering out NA values
data <- proteinshort %>%
  filter(!is.na(ddg_score), !is.na(new_aa_1l)) %>% # remove NA values
  group_by(new_aa_1l) %>%
  summarize(mean_ddg = mean(ddg_score)) %>%
  arrange(desc(mean_ddg))

# Convert new_aa_1l to a factor, excluding NA
data$new_aa_1l <- factor(aa_and_ddg$new_aa_1l, exclude = NA)

# Plot the bar plot with The Economist theme
p <- ggplot(data, aes(x = reorder(new_aa_1l, -mean_ddg), y = mean_ddg)) +
  geom_bar(stat = "identity", fill = "#4575b4") + # blue color for the bars
  geom_hline(yintercept = mean(data$mean_ddg), color = "red", size = 1.5, linetype = "dashed") + # add mean line
  labs(title = "Mean ddg scores grouped by the amino acid that was mutated", 
       x = "Residue type", 
       y = "Mean ddg score",
       subtitle = "The new (mutated) amino acid chains were group together by which amino acid was added.\nAfterwards, the mean ddg value was calculated for the new residue.\nThis shows that on average a single point mutation with P as the new residue increases the ddg value the most.\nConversely, if the new residue is G the mean ddg value is unlikely to change.") +
  scale_y_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 250)) + # set the y-axis limits and intervals
  theme_economist() + # use The Economist theme
  coord_flip() +  # flip the x and y axes
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(hjust = 0, margin = margin(b = 5)),
        axis.text = element_text(face = "bold"))

p 
# Save the plot as a high-resolution PNG file
ggsave(filename = file.path("charts", "mean_ddg_score_vs_new_aa.png"), plot = p, width = 8, height = 6, dpi = 300)







# Sort the data by mean_ddg in ascending order
# data <- data[order(desc(data$mean_ddg)),]
# 
# # Create a variable to store the mean of the mean_ddg column
# mean_ddg_mean <- mean(data$mean_ddg)
# 
# # Create the bar chart
# p <- ggplot(data, aes(x = factor(new_aa_1l, levels = new_aa_1l), y = mean_ddg)) +
#   geom_bar(stat = "identity", fill = "#F28F3B") +
#   scale_y_continuous(expand = c(0, 0), breaks = seq(0, 800, by = 100), limits = c(0, 800)) +
#   theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"),
#         axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 12),
#         panel.background = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.x = element_line(color = "gray90", size = 0.5),
#         panel.grid.minor.x = element_blank()) +
#   ggtitle("Mean ddg by previous amino acid") +
#   labs(x = "Previous amino acid", y = "Mean ddg") +
#   geom_hline(yintercept = mean_ddg_mean, color = "red", size = 1.5, linetype = "dashed") +
#   theme(plot.margin = unit(c(1, 1, 1, 0.5), "cm"))
# p

