#Volcanoplot

Volcano_plot_data <- as.data.frame(results)

ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant")), size = 2) +
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggtitle("Volcano Plot") +
  geom_point(aes(color = ifelse(abs(log2FoldChange) > 2, "Significant", "Not Significant")), size = 2) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "red")
