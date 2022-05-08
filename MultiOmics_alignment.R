library(gplots)
library(ggplot2)
library(reshape2)
library(visNetwork)
source('MONA_function.R')


##########Read Data ##########
omicsdatas <- list()
## omicsdatas = list of matrix, with each matrix represents each omicsdata

for (i in 1:5) {
  omicsdatas[[i]] <- read.delim(sprintf('ExampleData/timepoint%d.csv', i))
}

omicstypes <- c('rna', 'rna', 'rna', 'rna', 'rna')
## omicstypes = vector of omicsdata type
## available option are c('rna', 'pro', 'cnv', 'mirna')


##### Determine how many pairwise alignments are to be performed #####
alignments <- determine_alignment(omicstypes)



########## Run pairwise alignments ##########
align_outputs <- list()
for (i in 1:nrow(alignments)) {
  omic.x <- alignments$omic.x[i]
  omic.y <- alignments$omic.y[i]
  omic.x.type <- alignments$omic.x.type[i]
  omic.y.type <- alignments$omic.y.type[i]
  align.type <- alignments$align.type[i]
  align_outputs[[i]] <- align_omics(omicsdatas[[omic.x]], omicsdatas[[omic.y]], omic.x.type, omic.y.type, align.type)
}

matchings   <- list()
samplecors <- list()
genecors   <- list()
for (i in 1:length(align_outputs)) {
  genecors[[i]] <- align_outputs[[i]]$geneCorPlot_2
  samplecors[[i]] <- align_outputs[[i]]$sampleCor_2
  matchings[[i]]  <- align_outputs[[i]]$matching_2
}


########## Network alignment by looking at all pairwise alignment ##########
comm_id <- network_alignment(omicstypes, alignments[1:7, ], matchings, samplecors)

colnames(comm_id) <- paste0(omicstypes, '_', 1:length(omicstypes))
comm_id <- as.data.frame(comm_id)

test <- apply(comm_id, 2, function(x) x == 1:nrow(comm_id))
comm_id$mismatch <- apply(test, 1, function(x) FALSE %in% x)

for (i in 1:ncol(comm_id)) { print(which(1:nrow(comm_id) != comm_id[, i]))}



########## Visualizing aligned network ##########
## Edge represents the most likely pairing with different omicsdata
## Red edge represents the omicsdata are from different patients labels

nodes <- data.frame()
for (i in 1:length(omicstypes)) {
  mnode <- data.frame(id = sprintf('%s(%d)|%s(%d)', omicstypes[i], i, colnames(omicsdatas[[i]]), 1:ncol(omicsdatas[[i]])), stringsAsFactors = F)
  mnode$type <- omicstypes[i]
  mnode$omic_id <- i
  nodes <- rbind(nodes, mnode)
}


edges <- data.frame()
for (a in 1:nrow(alignments)) {
  omic.x      <- alignments$omic.x[a]
  omic.y      <- alignments$omic.y[a]
  omic.x.type <- alignments$omic.x.type[a]
  omic.y.type <- alignments$omic.y.type[a]
  omic.x.lbl  <- rownames(samplecors[[a]])
  omic.y.lbl  <- colnames(samplecors[[a]])
  
  aligncomm <- abs(comm_id[, c(omic.x, omic.y)])
  for (i in 1:length(omic.x.lbl)) {
    ox.id <- which(aligncomm[, 1] == i)
    oy.id <- which(aligncomm[, 2] == i)
    edgepairs <- merge(data.frame(ox = ox.id), data.frame(oy=oy.id))
    
    if (nrow(edgepairs) == 0)    next
    edgepairs$from <- sprintf('%s(%d)|%s(%d)', omic.x.type, omic.x, omic.x.lbl[edgepairs$ox], edgepairs$ox)
    edgepairs$to   <- sprintf('%s(%d)|%s(%d)', omic.y.type, omic.y, omic.y.lbl[edgepairs$oy], edgepairs$oy)
    edgepairs$matching.type <- apply(edgepairs, 1, function(x) if (x[1] == x[2])  0  else  1)
    edgepairs$corr <- apply(edgepairs, 1, function(x) samplecors[[a]][as.numeric(x[1]), as.numeric(x[2])])
    
    edges <- rbind(edges, edgepairs[, 3:6])
  }
}

visualize_network(nodes, edges)
