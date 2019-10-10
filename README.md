# Mining Magnetotactic Bacteria Abundance in Metagenomes Uploaded to MGRAST
Identifying magnetotactic bacteria (MTB) in microbiome projects from the MGRAST database
Special thanks to Victoria Hsiung for the original implementations of this code.

## Part 1:  Downloading metadata for all MGRast projects
First, we need to download all the metadata for the projects in the MGRast database.  At the time of analysis, June 26, 2018, there are 54376 samples in 2175 projects.  The API tools in the [R package 'matR'](https://github.com/MG-RAST/matR) were used.
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

The two R-data files can be downloaded from this GitHub repository:
- [ids.Rdata](./ids.Rdata)
- [MGRast-metadata.final.Rdata](./MGRast-metadata.final.Rdata)


## Part 2:  Download abundance counts for each MGRast sample using the python tools.

```bash
# Download and install the python tools
git clone http://github.com/MG-RAST/MG-RAST-Tools
cd MG-RAST-Tools
python setup.py build
sudo python setup.py install

# Loop through a file of sample IDs
c=1
while read i
   do
   mg-display-statistics.py --id "$i" --url=http://api.metagenomics.anl.go --stat species > $i.tsv
   echo "Finished sample number $c"
   c=$(( $c + 1 ))
done < file
```

Now process them in R  
```R
# Read in all tsv file names
tsv = list.files(".", pattern = ".tsv", full.names = F, recursive = T)

# Read in initial dataset
data = t(read.table(tsv[1], header = F, sep  ="\t", row.names = 1))
data = 10000 * data / sum(data)
data = data.frame(sample = sub(".tsv", "", tsv[1]), data)
rownames(data) = NULL

count = seq(1000, 55000, by = 1000)

# Merge with additional datasets
for (i in 2:length(tsv)){
   tmp = t(read.table(tsv[i], header = F, sep  ="\t", row.names = 1))
   tmp = 10000 * tmp / sum(tmp)
   tmp = data.frame(sample = sub(".tsv", "", tsv[i]), tmp)
   rownames(tmp) = NULL
   data = merge(data, tmp, all = T)
   if (i %in% count) save(data, file = "final.tmp.Rdata")
   r = nrow(data)
   c = ncol(data)
   print(paste0("Finished TSV file ", i, " ... ", r, " samples and ", c, " genera..."))
}

# Save final results
rownames(data) = data$sample
data$sample = NULL
save(data, file = "Full.species_table.Rdata")
write.table(data, file = "Full.species_table.tsv", quote = F, sep = "\t", header = T, row.names = F)

# Get subtable of only magnetotactic bacteria species
# load vector of species names
mtb = scan("MTB2", what = "character")
data.mtb = data[, which(colnames(data) %in% mtb)]
data.mtb = data.matrix(data.mtb)
data.mtb[is.na(data.mtb)] <- 0
data.mtb.clean = data.mtb[which(rowSums(data.mtb) > 0),]   # remove samples without any MTB

# Merge with metadata
out = as.data.frame(data.mtb.clean)
out$Total = rowSums(data.mtb.clean)
load("MGRast-metadata.final.Rdata")
tmp = metadata.ids[which(metadata.ids$Sample %in% rownames(out)), c(1,2,7,12,14,18,21,34,39,40,44,54,55,72)]
tmp = tmp[order(tmp$Sample),]
out = cbind(out[order(rownames(out)),], tmp)
save(out, file = "MTB.table.Rdata")
```

## Part 3:  Processing and plotting the normalized species counts
```R
load("MTB.table.Rdata")
data = out

# Principle components analysis
data.pca <- prcomp(data[,1:21], center = TRUE, scale = TRUE)
scores = as.data.frame(data.pca$x)
data = cbind(out, scores)

# Plot PCs
library(ggplot2)
ggplot(data, aes(x = PC1, y = PC2, color = metadata.env_package.data.env_package)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point()


# Build dendrogram and plot (select only "host-associated" samples)
library(dendextend)
data2 = subset(data, metadata.env_package.data.env_package == "host-associated")
d <- dist(data2[,1:21], method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
fit.d = as.dendrogram(fit)
fit.d %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 0.5) %>%
  set("leaves_col", as.numeric(as.factor(data2$metadata.sample.data.material))) %>%
  set("labels", "") %>%
  set("branches_k_color", k = 8) %>%
  plot
fit.d %>% rect.dendrogram(k = 8, border = "red", lty = 5, lwd = 2)
```

## __OLD CODE__
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
You can also download the file here:
- [mgrast.R](./mgrast.R)


Next, merge all the individual project's tsv files into a single, large data table for downstream use.
```R
# Read in all tsv file names
tsv = list.files(".", pattern = ".tsv", full.names = T, recursive = T)

# Read in initial dataset
data = read.table(tsv[1], header = T, sep  ="\t", row.names = 1)
data$samples = rownames(data); rownames(data) = NULL

# Merge with additional datasets
for (i in 2:length(tsv)){
   tmp = read.table(tsv[i], header=T, sep="\t", row.names = 1)
   tmp$samples = rownames(tmp); rownames(tmp) = NULL
   data = merge(data, tmp, all = T)
   save(data, file = "final.tmp.Rdata")
   r = nrow(data)
   c = ncol(data)
   print(paste0("Finished TSV file ", i, " ... ", r, " samples and ", c, " genera..."))
}

# Save final results
save(data, file = "final.table.Rdata")
write.table(data, file = "final.table.tsv", quote = F, sep = "\t", header = T, row.names = F)
```
