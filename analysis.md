# The aanlysis in R
After Jesse went through the full data table and edited all the metadata to be more descriptive
```R
# Copy data from EXCEL
# data <- read.table(pipe("pbpaste"), sep = "\t", header = TRUE)

# Using only host-associated samples
data = read.csv("/Users/rfitak/Dropbox/DUKE/MTB_MGRAST/data.csv", header = T)
   # 15811 rows, 39 columns

# Principle components analysis
data.pca <- prcomp(data[,18:38], center = TRUE, scale = TRUE)
scores = as.data.frame(data.pca$x)
data = cbind(data, scores)

# Plot PCs
library(ggplot2)
ggplot(data, aes(x = PC1, y = PC2, color = Scientific.Name..if.given.)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point()

```
