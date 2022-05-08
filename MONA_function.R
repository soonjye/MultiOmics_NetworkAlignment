#### Check input data
determine_alignment <- function(omicstypes) {
  combinations <- combn(length(omicstypes), 2)
  alignments <- data.frame()
  for (i in 1:ncol(combinations)) {
    omic.x <- combinations[1, i]
    omic.y <- combinations[2, i]
    omic.x.type <- omicstypes[omic.x]
    omic.y.type <- omicstypes[omic.y]
    
    if (omic.x.type == 'mirna' & !(omic.y.type %in% c('rna')))  next
    if (omic.y.type == 'mirna' & !(omic.x.type %in% c('rna')))  next
    
    if (omic.y.type == omic.x.type) {
      align.type <- 1  ## direct correlation
    } else if (omic.y.type == 'mirna' | omic.x.type == 'mirna') {
      align.type <- 2  ## cis-regulation
    } else {
      align.type <- 1  ## direct correlation
    }
    
    alignments <- rbind(alignments, data.frame(omic.x, omic.y, omic.x.type, omic.y.type, align.type, stringsAsFactors = F))
  }
  return(alignments)
}

check_omics_data <- function(omic_data) {
  if (sum(is.na(omic_data)) > 0) {
    stop('Missing value in the data matrix')
  }
  
  var_feature <- which(apply(omic_data, 1, sd) != 0)
  omic_data <- omic_data[var_feature, ]
  return(omic_data)
}

check_patient_sample <- function(pt_samples, omicstype){
  omics_num <- 0
  colname_check <- paste0('omics', omics_num+1) %in% colnames(pt_samples)
  
  while (colname_check) {
    omics_num <- omics_num + 1
    colname_check <- paste0('omics', omics_num+1) %in% colnames(pt_samples)
  }
  
  if (omics_num != length(omicstype)) {
    stop(sprintf('Found %d omics data but %d omics type are given', omics_num, length(omicstype)))
  }
  
  return(omics_num)
}

check_sample_label <- function(omic_data, sample_label) {
  data_label <- colnames(omic_data)
  
  if (length(data_label) != sum(!is.na(sample_label))) {
    stop(sprintf('Found %d columns in data but there are %d sample labels', length(data_label), sum(!is.na(sample_label))))
  }
  
  if (length(setdiff(data_label, sample_label)) != 0) {
    stop(sprintf('Sample in the data matrix is not found in the sample labels'))
  }
  
  if (length(setdiff(sample_label, data_label)) != 0) {
    stop(sprintf('Sample in the sample labels is not found in the data matrix'))
  }
  
}




#### Alignment: clean data input
computeCorrelatedGene <- function(omics.x, omics.y, cutoff = 0.5){
  feature_pair <- list()
  omics.x <- as.matrix(omics.x)
  omics.y <- as.matrix(omics.y)
  
  intergene <- intersect(colnames(omics.x), colnames(omics.y))
  cat(sprintf('Compute Pearson Correlation of %d genes (present in both dataset)... ', length(intergene)))
  gene_cor  <- rep(0, length(intergene))
  names(gene_cor) <- intergene
  for (gene in intergene) {
    gene_cor[gene] <- cor(omics.x[, gene], omics.y[, gene], method='pearson')
  }
  
  gene_cor <- gene_cor[order(gene_cor, decreasing = T)]
  gene_cor <- data.frame(Var1=names(gene_cor), Var2=names(gene_cor), corr=gene_cor, stringsAsFactors = F)
  colnames(gene_cor) <- c('Var1', 'Var2', 'corr')
  cat('Completed \n')
  
  feature_pair$geneCor <- gene_cor
  feature_pair$geneCorPlot <- ggplot(gene_cor, aes(x=corr)) +
    geom_histogram(color='black', fill='white') +
    geom_vline(xintercept = cutoff, linetype="dashed") +
    labs(x='Pearson Correlation', title='Distribution of Gene-wise Pearson Correlation',
         subtitle=sprintf('Mean = %.3f; SD = %.3f', mean(gene_cor$corr, na.rm=T), sd(gene_cor$corr, na.rm=T)))
  
  corgene <- gene_cor[gene_cor$corr > cutoff & !(is.na(gene_cor$corr)), ]
  feature_pair$correlated_pair <- corgene
  cat(sprintf('Selected %d genes (correlation > %.2f) \n', nrow(corgene), cutoff))
  
  return(feature_pair)
}


computeCisRegulation <- function(omics.x, omics.y, type.x, type.y) {
  type.x <- tolower(type.x)
  type.y <- tolower(type.y)
  
  omics.x <- as.matrix(omics.x)
  omics.y <- as.matrix(omics.y)
  
  if (type.x == 'mirna') {
    omics.mirna <- omics.x
    omics.non.mirna <- omics.y
  } else if (type.y == 'mirna') {
    omics.mirna <- omics.y
    omics.non.mirna <- omics.x
  }
  
  feature_pair <- list()
  nearbypair <- read.table('cohost_mirna.txt', sep='\t', header=T, stringsAsFactors = F)
  if (grepl('^ENSG', colnames(omics.non.mirna)[1]))    gene.id.type <- 'ensembl_gene_id'    else    gene.id.type <- 'hgnc_symbol'
  nearbypair <- nearbypair[which(nearbypair[, gene.id.type] %in% colnames(omics.non.mirna) & nearbypair$mirna %in% colnames(omics.mirna)), ]
  nearbypair <- unique(nearbypair[, c(gene.id.type, 'mirna')])
  nearbypair$Rcor <- NA
  nearbypair$Scor <- NA
  nearbypair$pvalue <- NA
  for (i in 1:nrow(nearbypair)) {
    nearbypair$Rcor[i] <- cor(as.numeric(omics.non.mirna[, nearbypair[i, 1]]), as.numeric(omics.mirna[, nearbypair[i, 2]]))
    nearbypair$Scor[i] <- cor(as.numeric(omics.non.mirna[, nearbypair[i, 1]]), as.numeric(omics.mirna[, nearbypair[i, 2]]), method = 'spearman')
    nearbypair$pvalue[i] <- cor.test(as.numeric(omics.non.mirna[, nearbypair[i, 1]]), as.numeric(omics.mirna[, nearbypair[i, 2]]))$p.value
  }
  nearbypair$qvalue <- p.adjust(nearbypair$pvalue, method="BH", n=nrow(nearbypair))
  
  feature_pair$geneCor <- nearbypair
  feature_pair$geneCorPlot <- ggplot(nearbypair, aes(x=Rcor)) +
    geom_histogram(color='black', fill='white') +
    geom_vline(xintercept = 0.2, linetype="dashed") +
    labs(x='Pearson Correlation', title='Distribution of miRNA vs host Gene Pearson Correlation',
         subtitle=sprintf('Mean = %.3f; SD = %.3f', mean(nearbypair$Rcor, na.rm=T), sd(nearbypair$Rcor, na.rm=T))) 
  
  corpairs <- nearbypair[which(nearbypair$Rcor > 0.2), ]
  corpairs <- corpairs[order(-corpairs$Rcor), ]
  
  if (type.x == 'mirna') {
    corpairs <- corpairs[, c(2, 1, 3:6)]
  }
  colnames(corpairs)[1:2] <- c('Var1', 'Var2')
  
  feature_pair$correlated_pair <- corpairs
  return(feature_pair)
}


computeRankdist <- function(corsample){
  rnadistpro <- corsample
  prodistrna <- corsample
  for (i in 1:nrow(corsample)){
    rnadistpro[i,] <- exp(scale(rnadistpro[i,])) / sum(exp(scale(rnadistpro[i,])))
    prodistrna[,i] <- exp(scale(prodistrna[,i])) / sum(exp(scale(prodistrna[,i])))
  }
  rankdist <- sqrt(rnadistpro * prodistrna)
  
  return(rankdist)
}


stableMarriage <- function(cormatrix){
  probmatrix <- computeRankdist(cormatrix)
  sampleN <- nrow(probmatrix)
  
  rankrow <- t(apply(probmatrix, 1, function(x) rank(-x, ties.method='first')))
  rankcol <- apply(probmatrix, 2, function(x) rank(-x, ties.method='first'))
  
  matches <- data.frame(omics1=1:nrow(probmatrix), omics1_label=NA, omics2=NA, omics2_label=NA)
  matches$self_cor   <- diag(cormatrix)
  matches$self_score <- diag(rankcol + rankrow)
  matches$match_cor   <- NA
  matches$match_score <- NA
  singlemales <- 1:nrow(probmatrix)
  femalerings <- rep(NA, ncol(probmatrix))
  
  while (length(singlemales) != 0){
    for (ppsmale in singlemales){
      propose <- 1
      single <- TRUE
      while (single == TRUE) {
        ppsfmle <- which(rankrow[ppsmale,] == propose)
        engaged <- femalerings[ppsfmle]
        if (is.na(engaged) || rankcol[engaged, ppsfmle] > rankcol[ppsmale, ppsfmle]) {
          matches$omics2[ppsmale] <- ppsfmle
          femalerings[ppsfmle]    <- ppsmale
          matches$match_score[ppsmale] <- rankrow[ppsmale, ppsfmle] + rankcol[ppsmale, ppsfmle]
          matches$match_cor[ppsmale]   <- cormatrix[ppsmale, ppsfmle]
          singlemales             <- setdiff(singlemales, ppsmale)
          if (!(is.na(engaged))) singlemales <- c(singlemales, engaged) 
          single <- FALSE
        } else {
          propose <- propose + 1
        }
      }
    }
  }
  
  matches$omics1_label <- rownames(probmatrix)[matches$omics1]
  matches$omics2_label <- colnames(probmatrix)[matches$omics2]
  
  nonmatch <- which(matches$omics1 != matches$omics2)
  nonmatch <- c(nonmatch, which(matches$omics1 == matches$omics2 & matches$match_score > max(2, sampleN*0.1)))
  matches$mismatch_status <- 0
  matches$mismatch_status[nonmatch] <- 1
  return(matches)
}


align_omics <- function(omicdata.x, omicdata.y, omictype.x, omictype.y, align.type, cutoff = 0.5) {
  align_output <- list()            # initialize output with a list of intermediate results
  omicdata.x <- scale(t(omicdata.x))
  omicdata.x <- t(na.omit(t(omicdata.x)))
  omicdata.y <- scale(t(omicdata.y))
  omicdata.y <- t(na.omit(t(omicdata.y)))
  
  ## if mirna is not being compared, then only extract overlapped features
  if (align.type == 1) {
    feature_label <- intersect(colnames(omicdata.x), colnames(omicdata.y))
    cat('Number of variables present in both omics data =', length(feature_label), '\n\n')
    
    omicdata.x <- omicdata.x[, feature_label]
    omicdata.y <- omicdata.y[, feature_label]
  }
  
  ## First iteration of matching
  cat('First iteration: Compute Gene Correlation ... \n')
  if (align.type == 2) {
    feature_cor <- computeCisRegulation(omicdata.x, omicdata.y, omictype.x, omictype.y)
  } else {
    feature_cor <- computeCorrelatedGene(omicdata.x, omicdata.y)
  }
  
  align_output <- c(align_output, feature_cor)
  names(align_output) <- paste0(names(align_output), '_1')
  
  corsample <- cor(t(omicdata.x[, align_output$correlated_pair_1$Var1]), t(omicdata.y[, align_output$correlated_pair_1$Var2]))
  align_output$sampleCor_1 <- corsample
  cat('First iteration: Perform data matching ... \n')
  
  matcher <- stableMarriage(corsample)
  align_output$matching_1 <- matcher
  
  nonmatch <- which(matcher$mismatch_status == 1)
  cat('First iteration: Found', length(nonmatch), 'pair mismatch\n\n')
  
  
  ## Second iteration of matching
  if (length(nonmatch) > 0) {
    if (align.type == 2) {
      feature_cor <- computeCisRegulation(omicdata.x[-nonmatch, ], omicdata.y[-nonmatch, ], omictype.x, omictype.y)
    } else {
      feature_cor <- computeCorrelatedGene(omicdata.x[-nonmatch, ], omicdata.y[-nonmatch, ])
    }
    
    names(feature_cor) <- paste0(names(feature_cor), '_2')
    align_output <- c(align_output, feature_cor)
    
    corsample <- cor(t(omicdata.x[, align_output$correlated_pair_2$Var1]), t(omicdata.y[, align_output$correlated_pair_2$Var2]))
    align_output$sampleCor_2 <- corsample
    
    cat('Second iteration: Perform data matching ... \n')
    matcher <- stableMarriage(corsample)
    align_output$matching_2 <- matcher
    
    nonmatch <- which(matcher$mismatch_status == 1)
    cat('Second iteration: Found', length(nonmatch), 'pair mismatch\n\n')
    
  } else {
    names(feature_cor) <- paste0(names(feature_cor), '_2')
    align_output <- c(align_output, feature_cor)
    align_output$sampleCor_2 <- align_output$sampleCor_1
    align_output$matching_2 <- align_output$matching_1
  }
  
  return(align_output)
}



network_alignment <- function(omicstypes, alignments, matchings, samplecors, pctthreshold = 0.9) {
  sampleN     <- nrow(samplecors[[1]])
  num_omics   <- length(omicstypes)
  
  comm_id     <- matrix(1:sampleN, nrow=sampleN, ncol=num_omics)
  comm_confd  <- matrix(0, nrow=sampleN, ncol=num_omics)
  comm_assign <- matrix(0, nrow=sampleN, ncol=num_omics)
  comm_prio   <- matrix(0, nrow=sampleN, ncol=num_omics)
  
  for (i in 1:ncol(comm_confd)) {
    alignids <- which(alignments$omic.x == i | alignments$omic.y == i)
    num_compare <- length(alignids)
    
    for (k in 1:num_compare) {
      alignid <- alignids[k]
      comm_confd[, i] <- comm_confd[, i] + matchings[[alignid]]$mismatch_status
      comm_prio[, i]  <- comm_prio[, i] + diag(samplecors[[alignid]])
    }
    
    comm_confd[, i] <- comm_confd[, i]/k
  }
  comm_confd <- 1 - comm_confd
  comm_assign[which(comm_confd == 1)] <- 1
  
  while (sum(comm_assign == 0) > 0) {
    tobeassigned <- as.data.frame(which(comm_assign != 1, arr.ind = TRUE))
    tobeassigned$val <- apply(tobeassigned, 1, function(x) comm_confd[x[1], x[2]])
    tobeassigned$priority <- apply(tobeassigned, 1, function(x) comm_prio[x[1], x[2]])
    tobeassigned <- tobeassigned[order(tobeassigned$priority, tobeassigned$val), ]
    i <- tobeassigned$col[1]
    j <- tobeassigned$row[1]
    
    aligncor <- matrix(0, nrow=sampleN, ncol=num_omics)
    for (k in 1:num_omics) {
      alignidx <- which(alignments$omic.x == i & alignments$omic.y == k)
      if (length(alignidx) > 0) {
        aligncor[, k] <- samplecors[[alignidx]][j, comm_id[,k]]
      }
      
      alignidy <- which(alignments$omic.x == k & alignments$omic.y == i)
      if (length(alignidy) > 0) {
        aligncor[, k] <- samplecors[[alignidy]][comm_id[,k], j]
      }
    }
    #aligncor <- t(apply(aligncor, 1, function(x) x*(comm_confd[j, ]+0.001)))
    corsums <- rowSums(aligncor)
    
    if (max(corsums) == 0)    pct <- 0    else   pct <- corsums[j] / max(corsums)
    potentiallbl <- which.max(corsums)
    singlenodes <- setdiff(1:sampleN, comm_id[, i] * comm_assign[, i])
    
    if (potentiallbl == j | pct >= pctthreshold | (!(potentiallbl %in% singlenodes) & pct > pctthreshold/2)) {
      comm_id[j, i] <- j
      cat(sprintf("Sample %d from omics %d: Retained from %d (self-align pct = %.2f)\n", j, i, potentiallbl, pct))
      comm_confd[j, i] <- sum(apply(aligncor, 2, function(x) x[j] == max(x)))/num_omics
      comm_assign[j, i] <- 1
      snatched <- setdiff(which(comm_id[, i] == j), j)
      if (length(snatched) > 0) {
        comm_assign[snatched, i] <- 0
        comm_id[snatched, i] <- snatched
        comm_confd[snatched, i] <- comm_confd[snatched, i] - 1/num_omics
      }
      pl <- j
    } else {
      comm_id[j, i] <- potentiallbl
      cat(sprintf("Sample %d from omics %d: label corrected to Sample %d (self-align pct = %.2f)\n", j, i, potentiallbl, pct))
      comm_confd[j, i] <- sum(apply(aligncor, 2, function(x) x[potentiallbl] == max(x)))/num_omics
      comm_assign[j, i] <- 1
      pl <- potentiallbl
    }
    
    
    # check if this increase confidence of potentiallbl
    iprime <- setdiff(which(comm_assign[pl, ] == 0), i)
    
    for (ip in iprime) {
      singlenodes <- tobeassigned$row[tobeassigned$col == ip]
      
      alignids <- which(alignments$omic.x == ip | alignments$omic.y == ip)
      num_compare <- length(alignids)
      
      selfalign <- 0
      for (k in 1:num_compare) {
        alignid <- alignids[k]
        if (alignments$omic.x[alignid] == ip) {
          relabel <- which(comm_id[, alignments$omic.y[alignid]] == pl)
          tocompare <- which.max(samplecors[[alignid]][pl, ])
          if (tocompare %in% relabel)    selfalign <- selfalign + 1
          
        } else if (alignments$omic.y[alignid] == ip) {
          relabel <- which(comm_id[, alignments$omic.x[alignid]] == pl)
          tocompare <- which.max(samplecors[[alignid]][, pl])
          if (tocompare %in% relabel)    selfalign <- selfalign + 1
        }
      }
      
      comm_confd[pl, ip] <- selfalign / k
      if (comm_confd[pl, ip] == 1) {
        comm_assign[pl, ip] <- 1
      }
    }
    cat('\n')
    
  }
  
  return(comm_id)
}



visualize_network <- function(vis.nodes, vis.links) {
  vis.nodes$title  <- vis.nodes$id # Text on click
  vis.nodes$label  <- vis.nodes$id # Node label
  vis.nodes$color.background <- substr(rainbow(max(vis.nodes$omic_id))[vis.nodes$omic_id], 1, 7)
  
  vis.links$color.color <- c(NA, "red")[vis.links$matching.type+1]
  vis.links$width <- c(NA, 5)[vis.links$matching.type+1]
  vis.links$title <- round(vis.links$corr, 4)
  visNetwork(vis.nodes, vis.links)
}
