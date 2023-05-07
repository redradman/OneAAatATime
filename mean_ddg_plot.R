library(tidyverse)
library(ggplot2)
library(ggthemes)
library(dplyr)

entireProtein <- read.csv("entireProtein.csv")
proteinshort <- entireProtein[, -c(7)]
proteinshort <- proteinshort[, (-c(11))]
# Create the data frame with the mean ddg scores, filtering out NA values
aa_and_ddg <- proteinshort %>%
  filter(!is.na(ddg_score), !is.na(new_aa_1l)) %>% # remove NA values
  group_by(new_aa_1l) %>%
  summarize(mean_ddg = mean(ddg_score)) %>%
  arrange(desc(mean_ddg)) %>%
  print(n = 20)

# Convert new_aa_1l to a factor, excluding NA
aa_and_ddg$new_aa_1l <- factor(aa_and_ddg$new_aa_1l, exclude = NA)

# Plot the bar plot with The Economist theme
p <- ggplot(aa_and_ddg, aes(x = reorder(new_aa_1l, -mean_ddg), y = mean_ddg)) +
  geom_bar(stat = "identity", fill = "#4575b4") + # blue color for the bars
  labs(title = "Mean ddg scores by residue type", x = "Residue type", y = "Mean ddg score") +
  scale_y_continuous(limits = c(0, 2000)) + # set the y-axis limits
  theme_economist() + # use The Economist theme
  coord_flip() # flip the x and y axes


# Save the plot as a high-resolution PNG file
ggsave("plot.png", p, width = 8, height = 6, dpi = 300)
