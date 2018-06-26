# MTB_MGRAST
Identifying magnetotactic bacteria (MTB) in microbiome projects from the MGRAST database
Special thanks to Victoria Hsiung for the original implementations of this code.

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

## Part 2:  Download normalized sequence counts for each sample
In this section, the raw sequence counts assigned to each genus are retrieved and normalized to a total of 10000 sequence reads. If the sequences are from a 16S amplicon metagenome, then the [RDP database](https://rdp.cme.msu.edu/) is used, otherwise the [RefSeq database](https://www.ncbi.nlm.nih.gov/refseq/) is used.  
The code, written and R, is implemented via the command line by simply giving the project number (`$PROJ`) as so:
```
Rscript mgrast.R $PROJ
```

In this way, it is easy to automate the above script on a cluster to analyze all 2176 projects.  The code in the `mgrast.R` script is:

```R
# Load library
library(matR)

# Custom funtion to get the mode
# This is used later the pick the most common sample type for each project
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Load metadata table made earlier
load("MGRast-metadata.final.Rdata")

# Get project ID (p)
args = commandArgs(trailingOnly = TRUE)
p = args[1]

# Get subset of samples within the project
s = as.vector(metadata.ids[which(metadata.ids$Project == p),]$Sample)

# Get the most common type of experiment for the samples in each project
# One of "WGS", "Amplicon", "Metabarcode", or "Other"
seqtype = getmode(metadata.ids[which(metadata.ids$Project == p),]$sequence_type)

# Get data for the samples in the project
if (seqtype == "WGS") {
   generaTable <- as.matrix(biomRequest(s, request = 'organism', hit_type = 'all', source = 'RefSeq', group_level = 'genus', evalue = 5, wait = TRUE))
} else {
   generaTable <- as.matrix(biomRequest(s, request = 'organism', hit_type = 'all', source = 'RDP', group_level = 'genus', evalue = 5, wait = TRUE))
}

# Normalize table to 10000 reads
out = t(apply(generaTable, 2, function(x) y = x / sum(x) * 10000))

# Save output to table
write.table(out, file = paste0(p, ".mgrast.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
```

