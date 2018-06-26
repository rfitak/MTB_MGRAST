library(matR)

# Code to search a particular set of projects in MGRast for magnetotactic bacteria
# arguments = 1) MGRast project, 2) count for output files

# Funtion to get the mode
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Load metadata table
load("/work/frr6/MGRAST/MGRast-metadata.final.Rdata")

# Get project ID and sequence type ("WGS", "Amplicon", "Metabarcode", "Other")
args = commandArgs(trailingOnly = TRUE)
p = args[1]

# Get subset of samples within project
s = as.vector(metadata.ids[which(metadata.ids$Project == p),]$Sample)

# Get the most common sequence type for a project
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

