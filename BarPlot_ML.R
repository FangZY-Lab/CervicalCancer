# Load required library
library(ggplot2)

# Create the data frame
data <- data.frame(
  Dataset = rep(c("GSE7803", "GSE9750", "GSE168244"), each = 3),
  Method = rep(c("KNN", "SVM", "RF"), times = 3),
  Accuracy = c(0.88, 0.90, 0.88, 0.85, 0.88, 0.88, 0.83, 0.86, 0.86)
)

# Optional: Define dummy error values for visualization (update with real values if available)
data$Error = c(0.02, 0.02, 0.02, 0.03, 0.02, 0.03, 0.03, 0.02, 0.02)
bar_colour=c("KNN" = "#e53528", "RF" = "#193e8f","SVM" = "grey")
# Create the grouped bar plot with error bars and grey background shades
ggplot(data, aes(x = Dataset, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5)  +
  scale_fill_manual(values = bar_colour) +
  labs(
    title = "",
    x = "",
    y = "Prediction Accuracy",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 1, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85"),
    panel.grid.minor.y = element_line(color = "gray90"),
    panel.background = element_rect(fill = "gray95", color = NA)
  )
