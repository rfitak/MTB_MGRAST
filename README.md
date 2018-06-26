# MTB_MGRAST
Identifying magnetotactic bacteria (MTB) in microbiome projects from the MGRAST database

## Part 1:  Downloading metadata for all MGRast projects
First, we need to download all the emtadata for the projects in the MGRast database.  At the time of analysis, June 26, 2018, there are 54376 samples in 2175 projects.  The API tools in the [R package 'matR'](https://github.com/MG-RAST/matR) were used.
```R
# Install the matR package
remove.packages('matR')        # if necessary
install.packages('devtools')
library(devtools)
install_github(repo='MG-RAST/matR') 
library(matR)
dependencies()

# Load libraries
library(matR)
library(reshape2)

# Build file of Project and Sample IDs, Some projects have 0 samples, and are ommitted.
projs = rownames(dir.MGRAST())
ids = metadata(projs)
ids = melt(ids[which(lapply(ids, length) > 0)])
ids = apply(ids, 2, as.character)
colnames(ids) = c("Sample", "Project")
save(ids, file = "ids.Rdata")

# Build empty table to store metadata
names = colnames(metadata(ids[1,1], detail = "metadata"))
out = as.data.frame(matrix(nrow = nrow(ids), ncol = length(names)))
colnames(out) = names

# Loop through each sample, grabbing the metadata information then making a final data table
for (i in 1:nrow(out)){
   z = metadata(ids[i,1], detail = "metadata")
   z.match = which(colnames(z) %in% names)
   names.match = which(names %in% colnames(z))
   out[i, names.match] <- z[1, z.match]
   print(paste0("Finished sample ", i))
   if (i %in% seq(from = 0, to = nrow(out), by = 1000)) save(out, file = "MGRast-metadata.Rdata")
}

# Combine the metadata with the Project and Sample IDs table and save
metadata.ids = cbind(ids, out)
save(metadata.ids, file = "MGRast-metadata.final.Rdata")
```
