require(ggplot2)
# calculate p-value by using t_test
rma.DEF <- function(sample1, sample2, exp_data, samples=NULL, p_cutoff=0.05, FC_cutoff=2, isPaired=FALSE) {
  #Parameters:
  # sample1: the sample name of 1st condition
  # sample2: the sample name of 2nd condition
  # exp_data: expression data of all samples
  # samples: all sample names
  # p_cutoff: Threshold of p value
  # FC_cutoff: Threshold of Fold change
  # isPaired: boolen. Does sample1 and sample2 paired
  #Returns:
  # A list containing significant gene set (sigGene) and its coressponding expr data (sigData), the down and up 
  # regulate gene set (DownRegulatedGene, UpRegulatedGene), also a heatmap (Heatmap) and a volcano image (Volcano)
  
  #sample_new_index = c(rep(0,length(sample1)),rep(1,length(sample2)))
  sample_new = c(sample1,sample2)
  exp_data_new = exp_data[,sample_new]
  
  n <- nrow(exp_data_new)
  p_values <- rep(0,n)
  FC_values <- rep(1,n)
  probset = rownames(exp_data_new)
  names(p_values) = probset
  names(FC_values) = probset
  
  nSample1 <- length(sample1); nSample2 <- length(sample2)
  for (i in 1:n) {
    if((sum(duplicated(exp_data_new[i,sample1])) == nSample1-1) | (sum(duplicated(exp_data_new[i,sample2])) == nSample2-1)) {
      bool <- FALSE
    } else {
      bool <- (shapiro.test(exp_data_new[i,sample1])$p.value < 0.05) | (shapiro.test(exp_data_new[i,sample2])$p.value < 0.05)
    }
    
    if(bool) {
      p_values[i] <- wilcox.test(exp_data_new[i,sample1],exp_data_new[i,sample2],paired = isPaired)$p.value
    } else {
      if(var.test(exp_data_new[i,sample1], exp_data_new[i,sample2])$p.value < 0.05) {
        p_values[i] <- t.test(exp_data_new[i,sample1],exp_data_new[i,sample2], paired = isPaired, var.equal = FALSE)$p.value
      } else {
        p_values[i] <- t.test(exp_data_new[i,sample1],exp_data_new[i,sample2], paired = isPaired, var.equal = TRUE)$p.value
      }
    }
    # pay attention this
    FC_values[i] <- (mean(exp_data_new[i,sample1]) / mean(exp_data_new[i,sample2]))
  }
  
  par(mfrow = c(1,2))
  hist(-log10(p_values),breaks = 100,main = 'P value')
  hist(log2(FC_values),breaks = 100,main = 'FC')

  threshold <- as.factor((FC_values > FC_cutoff | FC_values < 1/FC_cutoff) & p_values < p_cutoff)
  gplt <- ggplot(data=NULL, aes(x=log2(FC_values), y=-log10(p_values), colour=threshold))+geom_point()+labs(title="Ranksum test Volcanoplot",x="log2(fold change)", y="-log10(p value)")
  
  pos <- intersect(which(p_values < p_cutoff),union(which(FC_values > FC_cutoff),which(FC_values < 1/FC_cutoff)))
  sig_data <- as.data.frame(exp_data_new[pos,])
  sig_data$AveExpr1 <- rowMeans(sig_data[, sample1])
  sig_data$AveExpr2 <- rowMeans(sig_data[, sample2])
  
  sig_Gene <- cbind(data.frame(sig_data[, c('AveExpr1', 'AveExpr2')], FC_values = FC_values[pos], p_values = p_values[pos]))
  
  down_regulated_gene <- probset[intersect(which(p_values < p_cutoff),which(FC_values < 1/FC_cutoff))]
  up_regulated_gene <- probset[intersect(which(p_values < p_cutoff),which(FC_values > FC_cutoff))]
  
  hm <- heatmap(
    x = as.matrix(sig_data[,sample_new]), scale = "row", Rowv = NULL, Colv = NULL, margins = c(15,15),
    col = colorRampPalette(c("green", "black", "red"))(100),
    labRow = FALSE#, labCol = as.list(samples[which(samples %in% colnames(sig_data))])
  )
  
  return(list(sigGene=sig_Gene, sigData=sig_data[,1:6], DownRegulatedGene=down_regulated_gene, 
              UpRegulatedGene=up_regulated_gene, Heatmap=hm, Volcano=gplt))
}

# one-way anova
perform.aov <- function(product, Conditon) {
  #Parameters:
  # product: the data used to perform anova of one obj
  # Condition: a factor containing status of all samples
  allCond <- levels(Conditon)
  dup_status <- FALSE
  for (Cond in allCond) {
    LEN <- sum(Conditon == Cond) - 1
    dup_status <- dup_status | (sum(duplicated(product[Conditon == Cond])) == LEN)
  }
  
  if(dup_status) {
    bool <- FALSE
  } else {
    bool <- FALSE
    for (Cond in allCond) {
      bool <- bool | (shapiro.test(product[Conditon == Cond])$p < 0.05)
    }
  }
  
  if(bool) {
    hsd <- NA
    F_value <- NA
    p_value <- kruskal.test(x = product, g = Conditon)$p.value
  } else {
    if(bartlett.test(product~Conditon)$p.value < 0.05) {
      fit <- oneway.test(product~Conditon)
      hsd <- NA
      p_value <- fit$p.value
      F_value <- fit$statistic[[1]]
    } else {
      fit <- aov(product~Conditon)
      hsd <- TukeyHSD(fit)
      p_value <- summary(fit)[[1]][1,5]
      F_value <- summary(fit)[[1]][1,4]
    }
  }
  return(list(P.value=p_value, F.value=F_value, HSD=hsd))
}

# calculate p-value by using anova
rma.aovDEF <- function(exp_data, Conditon, cutoff = 0.05) {
  #Parameters:
  # exp_data: expression data of all samples
  # Condition: status of all samples
  # cutoff: Threshold of p value
  #Returns:
  # A list containing significant gene set (sigGene) and its coressponding expr data (sigData), also HSD (sigHSD)
  n <- nrow(exp_data)
  p_values <- rep(0,n)
  F_values <- rep(1,n)
  HSDs <- rep(1,n)
  probset = rownames(exp_data)
  names(p_values) = probset
  names(F_values) = probset
  
  for (i in 1:n) {
    value <- perform.aov(exp_data[i,], Conditon)
    p_values[i] <- value[[1]]
    F_values[i] <- value[[2]]
    HSDs[i] <- value[[3]]
    # print(HSDs[i])
  }
  
  par(mfrow = c(1,2))
  hist(-log10(p_values),breaks = 100,main = 'P value')
  hist(F_values,breaks = 100,main = 'F_value')
  
  pos <- which(p_values < cutoff)
  sig_data <- as.data.frame(exp_data[pos,])
  sig_Gene <- cbind(data.frame(F_values = F_values[pos], p_values = p_values[pos]))
  sig_HSD <- HSDs[pos]
  
  return(list(sigGene=sig_Gene, sigData=sig_data, sigHSD=sig_HSD))
}
