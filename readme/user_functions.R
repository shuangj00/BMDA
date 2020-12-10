BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
  }
  return(PPI_sorted[k])
}

# =============================================================================================
# The function of estimating size factors from count data
# ---------------------------------------------------------------------------------------------
# Input:  Y, a n-by-p count matrix, where n is the number of samples and p is the number of 
#         features
# Input:  method, a categorical variable chosen from
#         TC:   total counts
#         CSS:  cumulative-sum scaling by Paulson et al., 2013
#         TMM:  trimmed mean by M-values by Robinson and Oshlack, 2010
#         RLE:  relative log expression by Anders and Huber, 2010
#         Q75:  the 75th percentile by Bullard et al., 2010
#         none: all size factors are set to 1
# Input:  constraint, a boolean variable, with TRUE rescaling size factors to their sum equal 1
# ---------------------------------------------------------------------------------------------
# Output: s, a n-dimensional vector, where each element is the estimated size factor for the 
#         corresponding sample
# =============================================================================================
size_factor_estimator = function(Y, method = c("TSS", "CSS", "TMM", "RLE", "Q75")) {
  n <- dim(Y)[1];
  s <- rep(NA, n);
  if (method == "TSS") {
    s <- rowSums(Y);
  } else if (method == "CSS") {
    p <- cumNormStatFast_copy(t(Y));
    for (i in 1:n) {
      temp <- Y[i, which(Y[i,] != 0)];
      s[i] <- sum(Y[i, which(Y[i,] <= quantile(temp, p))]);
    }
  } else if (method == "TMM") {
    s <- calcNormFactors_copy(t(Y), method = method, refColumn = NULL, logratioTrim = 0.3, 
                              sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10);
  } else if (method == "RLE") {
    s <- calcNormFactors_copy(t(Y), method = method);
    # }
  } else if (method == "Q75") {
    s <- apply(Y, 1, quantile, 0.75);
    s <- pmax(s, 1);
    # s <- calcNormFactors_copy(t(Y), method = "upperquartile");
  }
  # if (constraint) {
  #   s <- s/sum(s);
  # }
  names(s) <- rownames(Y);
  return(s);
}
# =============================================================================================
# The functions attached to size_factor_estimator(...)
# =============================================================================================
rescale <- function(x, constraint = c("sum", "product")) {
  if (constraint == "sum") {
    return(x/sum(x));
  } else if (constraint == "product") {
    return(x/exp(mean(log(x))))
  }
} 
size_normalization = function(x){
  log_si = log(x) - mean(log(x))
  exp(log_si)
}
# Download from https://github.com/HCBravoLab/metagenomeSeq/blob/master/R/cumNormStatFast.R on 
# April 12, 2018
cumNormStatFast_copy <- function(mat, pFlag = FALSE, rel = 0.1, ...){
  smat = lapply(1:ncol(mat), function(i) {
    sort(mat[which(mat[, i] > 0), i], decreasing = TRUE)
  })
  leng = max(sapply(smat, length))
  if(any(sapply(smat, length) == 1)) {
    stop("Warning sample with one or zero features")
  }
  smat2 = array(NA, dim = c(leng, ncol(mat)))
  for(i in 1:ncol(mat)){
    smat2[leng:(leng - length(smat[[i]]) + 1), i] = smat[[i]]
  }
  rmat2 = sapply(1:ncol(smat2), function(i){
    quantile(smat2[, i], p = seq(0, 1, length.out = nrow(smat2)), na.rm = TRUE)
  })
  smat2[is.na(smat2)] = 0
  ref1 = rowMeans(smat2)
  ncols = ncol(rmat2)
  diffr = sapply(1:ncols, function(i) {
    ref1 - rmat2[, i]
  })
  diffr1 = matrixStats::rowMedians(abs(diffr))
  if(pFlag == TRUE){
    plot(abs(diff(diffr1))/diffr1[-1], type = "h", ...)
    abline(h = rel)
    axis(1, at = seq(0, length(diffr1), length.out = 5), labels = seq(0, 1, length.out = 5))
  }
  x= which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
  if(x <= 0.50){
    message("Default value being used")
    x = 0.50
  }
  return(x)
}

# Download from https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R on April 12, 2018
calcNormFactors_copy <- function(object, lib.size = NULL, method = c("TMM", "RLE", "upperquartile",
                                                                     "none"), 
                                 refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, 
                                 doWeighting = TRUE, Acutoff = -1e10, p = 0.75, ...) {
  x <- as.matrix(object)
  if (any(is.na(x))) {
    stop("NA counts not permitted")
  }
  if (is.null(lib.size)) {
    lib.size <- colSums(x)
  }
  if (any(is.na(lib.size))) {
    stop("NA lib.sizes not permitted")
  }
  method <- match.arg(method)
  allzero <- .rowSums(x>0, nrow(x), ncol(x)) == 0
  if(any(allzero)) {
    x <- x[!allzero, , drop = FALSE]
  }
  if(nrow(x) == 0 || ncol(x) == 1) {
    method = "none"
  }
  f <- switch(method,
              TMM = {
                f75 <- .calcFactorQuantile_copy(data = x, lib.size = lib.size, p = 0.75)
                if(is.null(refColumn)) {
                  refColumn <- which.min(abs(f75 - mean(f75)))
                }
                if(length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x)) {
                  refColumn <- 1
                }
                f <- rep(NA, ncol(x))
                for(i in 1:ncol(x)) {
                  f[i] <- .calcFactorWeighted_copy(obs = x[, i], ref = x[, refColumn], 
                                                   libsize.obs = lib.size[i], 
                                                   libsize.ref = lib.size[refColumn], 
                                                   logratioTrim = logratioTrim, sumTrim = sumTrim, 
                                                   doWeighting = doWeighting, Acutoff = Acutoff)
                }
                # f
                f*lib.size
              },
              RLE = .calcFactorRLE_copy(x),
              upperquartile = .calcFactorQuantile_copy(x, lib.size, p = p),
              none = rep(1, ncol(x))
  )
  f
}

.calcFactorRLE_copy <- function (data) {
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile_copy <- function (data, lib.size, p = 0.75)
{
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
}

.calcFactorWeighted_copy <- function(obs, ref, libsize.obs = NULL, libsize.ref = NULL, 
                                     logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
                                     Acutoff = -1e10) {
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs)) {
    nO <- sum(obs) 
  } else {
    nO <- libsize.obs
  }
  if (is.null(libsize.ref)) {
    nR <- sum(ref)
  } else {
    nR <- libsize.ref
  }
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  if(max(abs(logR)) < 1e-6) {
    return(1)
  }
  n <- length(logR)
  loL <- floor(n*logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n*sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)
  if (doWeighting) {
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep], na.rm = TRUE)
  } else {
    f <- mean(logR[keep], na.rm=TRUE)
  }
  if(is.na(f)) {
    f <- 0
  }
  2^f
}



# =============================================================================================
# The function of comparing ZINB-DPP model with the Kruskal-Wallis test
# ---------------------------------------------------------------------------------------------
# Input:  taxa.plot: name of the taxon that we want to compare the performance of 2 methods
#         HH: n-by-p matrix of the estimated zero-inflated index PPI
#         AA: n-by-p matrix of the estimated log alpha(ij)
# ---------------------------------------------------------------------------------------------
# Output: a figure compareing the abundance of the given taxon between 2 patient groups
#         
# =============================================================================================
log_dist_compare = function(taxa.plot, HH = Hmat, AA = Amat_est, z_v = z_v){
  check.nam.ind = which(taxa.output.name == taxa.plot)
  check.nam.ind_comp = which(colnames(Y_composition) == taxa.plot)
  
  data_temp = data.frame(Group = factor(z_v, levels = c(0, 1), labels = c("nonCRC", "CRC")), 
                         Abundance = log(AA[,check.nam.ind ]), 
                         HH = HH[,check.nam.ind])
  Bayes = ggplot(aes(x = Group, y = Abundance), data = subset(data_temp,HH<0.5 ) ) + 
    geom_violin(position=position_dodge(width=0.5)) + 
    stat_summary(fun.y = mean, geom = "point", size = 2, colour = "red")+
    xlab("") + ylab(expression(atop("Latent relative abundance", paste("in log scale log(", alpha["ij"], ")"))))
  
  # Relative Abundance #
  
  data_compo = data.frame(Group = factor(z_v, levels = c(0, 1), labels = c("nonCRC", "CRC")), Fraction = Y_composition[,check.nam.ind_comp]+ 1e-12)
  compostional = ggplot(aes(x = Group, y = Fraction), data = data_compo) + 
    geom_violin(width = 1,position=position_dodge(width=0.5)) + 
    stat_summary(fun.y = median, geom = "point", size = 2, colour = "blue") +
    xlab("") +
    ylab(expression(paste("Composition y"['ij'],"/y"['.j']))) +
    scale_y_continuous(trans='log10')
  
  title1 = textGrob(taxa.plot, gp=gpar(fontsize = 20, fontface="bold"))
  title2 = textGrob(paste0( "Kruskal-Wallis test: p-value = ", round(pval_nonpara[check.nam.ind_comp], 4), " \n ZINB-DPP: PPI = ", round(PPI_keep[check.nam.ind], 3)),
                    gp=gpar(fontsize = 20))
  ggres = grid.arrange(compostional, Bayes, nrow = 1, top = title1, 
                       bottom = title2)
  return(ggres)
}


# =============================================================================================
# The function to generate the adjacent matrix based on the taxonomic tree structure
# ---------------------------------------------------------------------------------------------
# Input:  Y: the sample-by-taxon count matrix for the bottom-level (e.g., species level) taxa
#         S: the taxonomic tree information matrix (each row is a higher-level taxon, and the
#            entry is the corresponding column index of lower level taxa)
# ---------------------------------------------------------------------------------------------
# Output: a large adjacent matrix representing the taxanomic tree
#         
# =============================================================================================

S2adj = function(Y, S) {
  taxa_names <- c(colnames(Y), rownames(S));
  levels <- c("k", "p", "c", "o", "f", "g", "s");
  p <- length(taxa_names);
  adj <- matrix(0, nrow = p, ncol = p);
  colnames(adj) <- taxa_names;
  rownames(adj) <- taxa_names;
  SS <- cbind(S, rep(0, dim(S)[1]));
  SS <- rbind(cbind(1:dim(Y)[2], matrix(0, nrow = dim(Y)[2], ncol = dim(SS)[2] - 1)), SS);
  rownames(SS) <- taxa_names;
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      level_i <- which(levels == substr(taxa_names[i], 1, 1));
      level_j <- which(levels == substr(taxa_names[j], 1, 1));
      if (abs(level_i - level_j) == 1) {
        taxa_i <- unique(SS[taxa_names[i],]);
        taxa_j <- unique(SS[taxa_names[j],]);
        if (length(intersect(taxa_i, taxa_j)) == length(taxa_i) || length(intersect(taxa_i, taxa_j)) == length(taxa_j)) {
          adj[i, j] <- 1;
          adj[j, i] <- 1;
        }
      }
    }
  }
  return(adj);
}

# =============================================================================================
# The function of generating the cladogram to visualize the selected taxa on the phylogenetic tree
# ---------------------------------------------------------------------------------------------
# Input:  controln: name of the control-enriched taxa
#         casen: name of the case-enriched taxa
#         hl_size: size of the highlighted taxon on the cladogram
#         line_wd: line width of the tree branch
# ---------------------------------------------------------------------------------------------
# Output: cladogram plot
#         
# =============================================================================================

tree.check = function(controln, casen, hl_size = 4, 
                      line_wd = 0.5){
  library(cowplot)
  library(ggpubr)
  library(ggtree)
  library(data.tree)
  library(microbiomeViz)
  if (!requireNamespace("microbiomeViz", quietly = TRUE)) {
    stop("Package \"microbiomeViz\" needed for this function to work. Please install it (devtools::install_github(\"zhanxw/microbiomeViz\")).",
         call. = FALSE)
  }
  
  taxa = taxa.full
  if (class(taxa) == "taxonomyTable") {
    taxa <- data.frame(taxa@.Data)
  }
  
  physeq = phyloseq(tax_table(taxa.full)) %>% fix_duplicate_tax()
  tr = parsePhyloseq(physeq, use_abundance = FALSE, node.size.scale = 0.00001)
  p = tree.backbone(tr, size = line_wd) # size: branch width
  if (is.null(controln)) {
    print(p)
    return(p)
  }
  
  # Specify node attributes to highlight
  dd <- data.frame(taxa = controln)
  dd2 = data.frame(taxa = casen)
  dd = rbind(dd, dd2)
  
  # generate plot #
  p <- p %<+% dd
  dd.node1 <- p$data$node[p$data$label %in% controln]
  dd.node2 <- p$data$node[p$data$label %in% casen]
  # node1: control; node2: case #
  p <- p +
    geom_point2(aes(subset = (node %in% dd.node2), col = "royalblue2", size = hl_size, alpha = 1)) + 
    geom_point2(aes(subset = (node %in% dd.node1), col = "maroon3", size = hl_size, alpha = 1)) +
    scale_size_identity() + scale_alpha_identity() + scale_color_identity()
  plot(p)
}


# =============================================================================================
# The function of generating the marginal posterior probability of inclusion (PPI)
# ---------------------------------------------------------------------------------------------
# Input:  gamma.PPI: estimated PPIs given by the model
#         tax.name.full: all names of the taxa shown in the PPI plot
#         (must be in the same order with the taxa in the gamma.PPI)
#         alpha.sig: Bayesian false discovery rate
# ---------------------------------------------------------------------------------------------
# Output: PPI plot
#         
# =============================================================================================


gamma.PPI.draw.pick = function(gamma.PPI, tax.name.full, pick.name, alpha.sig = 0.05){
  library(ggplot2)
  p.size = length(tax.name.full)
  dim(gamma.PPI) = c( 1,length(gamma.PPI))
  if(is.null(colnames(gamma.PPI))){
    colnames(gamma.PPI) = tax.name.full
    taxa.nam.ori = tax.name.full
  }else{
    taxa.nam.ori = colnames(gamma.PPI)
  }
  # (1)order the taxa by taxonomic level
  rank.od = c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  ind.by = match(substring(taxa.nam.ori, 1, 3), rank.od)
  ord = order(ind.by)
  ord.nam = taxa.nam.ori[ord]
  PPI.order = gamma.PPI[ord]; dim(PPI.order) = c( 1,length(gamma.PPI)); colnames(PPI.order) = ord.nam
  
  # 
  rank.group.size = rep(NA, length(rank.od))
  for(i in 1:length(rank.od)){
    rank.group.size[i] = length(grep(rank.od[i], ord.nam, value=TRUE))
    x.str = c(0, cumsum(rank.group.size)+1)[1:length(rank.od)] 
    x.end = cumsum(rank.group.size) + 1
  }
  x.taxa1  = grep(pick.name[1], ord.nam);x.taxa2  = grep(pick.name[2], ord.nam);
  rect.df = data.frame(x1 = x.str, x2 = x.end, y1 = rep(0, length(rank.od)), y2 = rep(1, length(rank.od)),
                       Rank = c("kingdom", "phylum", "class", "order", 
                                "family", "genus", "species"))
  
  cutoff.value = BayFDR(gamma.PPI, alpha.sig)
  cut.below = gamma.PPI[gamma.PPI == cutoff.value]
  cutoff.draw = max(cut.below) + 0.001
  ggplot() +
    geom_rect(data=rect.df,
              aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Rank), 
              alpha=0.35, size = 0) +
    geom_segment( color='gray60',aes(x=seq_along(PPI.order),
                                     xend=seq_along(PPI.order), y=0, yend = as.vector(PPI.order))) +
    geom_point(aes(x=seq_along(PPI.order), y= as.vector(PPI.order))) +
    xlab("Taxon Index") + 
    ylab(expression(paste("Posterior probability of inclusion ",pi,"(",gamma[j], "=1|.)"))) +
    #ggtitle(expression(paste("Marginal Posterior Probability of Inclusion for ",gamma))) +
    geom_hline(aes(yintercept = cutoff.draw), colour="#990000", linetype="dashed") +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)), 
          axis.title.y = element_text(size = rel(1.3)), 
          legend.title=element_text(size=rel(1.5)),
          legend.text=element_text(size=rel(1.5))) +
    scale_fill_brewer( limits=c("kingdom","phylum","class","order", "family", "genus","species"), 
                       palette="Spectral") + 
    geom_segment(aes(x=x.taxa1, xend=x.taxa1, y=1.06, yend=1.02), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "royalblue2") + 
    geom_segment(aes(x=x.taxa2, xend=x.taxa2, y=1.06, yend=1.02), 
                 arrow = arrow(length = unit(0.2, "cm")), color = 'maroon3')
  
  
}

# =============================================================================================
# The function of generating the 95% cedible interval of log(alpha_ij) given by the ZINB-DPP 
# ---------------------------------------------------------------------------------------------
# Input:  PPI: estimated PPIs given by the model
#         log_r_mat: an output from the ZINB-DPP model, with each row being a result from a single
#                    iteration. Each element in the row represents the estimated group difference, 
#                    which is log(alpha_j1/alpha_j0)
#         tax.name.full: all names of the taxa shown in the PPI plot
#         (must be in the same order with the taxa in the gamma.PPI)
#         alpha.sig: Bayesian false discovery rate
#         other parameters are used in generating the final figure
# ---------------------------------------------------------------------------------------------
# Output: PPI plot
#         
# =============================================================================================
CIalpha.ggplot = function(PPI, tax.name.full, alpha.mat, alpha.sig = 0.05, 
                          inputname.1 = "Control", inputname.2 = "Case", shift_right = 5.5,
                          shift_ang = 30, size_y = 1, taxa.font = 1.5){
  cutoff = BayFDR(PPI, alpha.sig)
  taxon.sel = which(PPI > cutoff)
  taxon_name = tax.name.full[taxon.sel]
  alpha.sele = alpha.mat[,taxon.sel ]
  alpha.sele = as.matrix(alpha.sele)
  taxon.plot.df = data.frame(taxon_name = taxon_name,
                             taxon_mean = apply(alpha.sele, 2,mean),
                             taxon_ql = apply(alpha.sele, 2, quantile,probs = 0.025),
                             taxon_qu = apply(alpha.sele, 2, quantile,probs = 0.975));
  plot.df = taxon.plot.df[with(taxon.plot.df, order(-taxon.plot.df$taxon_mean)), ]
  plot.df$taxon_name <- factor(plot.df$taxon_name, levels = plot.df$taxon_name)
  pos.count = sum(plot.df$taxon_mean > 0)
  neg.count = sum(plot.df$taxon_mean < 0)
  plot.df$sign = c(rep("pos", pos.count),  rep("neg", neg.count))
  if(length(unique(plot.df$sign)) !=1 ){
    ggplot(data=plot.df, aes(y=taxon_name, x=taxon_mean, xmin=taxon_ql, xmax=taxon_qu, color = sign)) + 
      geom_errorbarh(height=0.2, size=1) +
      geom_point(data=plot.df,aes(y=taxon_name, x=taxon_mean), size=2, shape=21, fill = 'white') +
      scale_y_discrete(limits = rev(levels(plot.df$taxon_name))) +
      geom_vline(xintercept=0, linetype="dotted") +
      theme(
        plot.margin=unit(c(5.5, 5.5, 5.5, shift_right), "points"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = shift_ang,face = "italic", hjust = 1, size = rel(taxa.font)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(size_y)),
        legend.title=element_text(size=rel(1.5)),
        legend.text=element_text(size=rel(1.5))
      ) +
      guides(fill=FALSE) +
      labs(x = expression(paste("Posterior effect size log(", alpha[j*2], "/", alpha[j*1], "|.)")),
           y = " ", color = " ") +
      scale_color_manual(labels = c(paste0("Taxa enriched in \n ", inputname.1, " group"), 
                                    paste0("Taxa enriched in \n", inputname.2, " group")), 
                         values = c("maroon3", "royalblue2"))+
      coord_flip()
  }else{
    col.idx = ifelse(unique(plot.df$sign) == "neg","maroon3", "royalblue2" )
    plt.label = ifelse(unique(plot.df$sign) == "neg", 
                       paste0("Taxa enriched in \n ", inputname.1, " group"), 
                       paste0("Taxa enriched in \n", inputname.2, " group"))
    ggplot(data=plot.df, aes(y=taxon_name, x=taxon_mean, xmin=taxon_ql, xmax=taxon_qu, color = sign)) + 
      geom_errorbarh(height=0.2, size=1) +
      geom_point(data=plot.df,aes(y=taxon_name, x=taxon_mean), size=2, shape=21, fill = 'white') +
      scale_y_discrete(limits = rev(levels(plot.df$taxon_name))) +
      geom_vline(xintercept=0, linetype="dotted") +
      theme(
        plot.margin=unit(c(5.5, 5.5, 5.5, shift_right), "points"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = shift_ang,face = "italic", hjust = 1, size = rel(taxa.font)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(size_y)),
        legend.title=element_text(size=rel(1.5)),
        legend.text=element_text(size=rel(1.5))
      ) +
      guides(fill=FALSE) +
      labs(x = expression(paste("Posterior effect size log(", alpha[j*2], "/", alpha[j*1], "|.)")),
           y = " ", color = " ") +
      scale_color_manual(labels = plt.label, 
                         values = col.idx)+
      coord_flip()
  }
  
}

gen_zinb_raw = function(zt, p_diff, p_remain, direction_t, seed,
                        low = 1, high = 2, phi_dispersion = 10, zero_pct = 0.5){
  if(length(direction_t) != p_diff){
    stop("check input information")
  }
  n = length(zt);
  K = 2
  N = floor(runif(n, 5000, 10000));
  set.seed(seed);
  
  # size factor #
  s = runif(n, 0.5, 4);
  # dispersion parameter psi in the NB distribution #
  psi_main = rexp(p_remain, 1/phi_dispersion);
  psi_diff = rexp(p_diff, 1/phi_dispersion);
  
  # mean vector #
  Mt = matrix(NA, nrow = K, ncol = p_remain);
  Mt[, ] = rep(runif(p_remain, 0, 2), each = K);
  
  Mt_diff = matrix(NA, nrow = K, ncol = p_diff);
  dir_count = 1
  for(j in 1:p_diff){
    sigma_b = runif(1, low, high)
    temp = seq(1 - (K - 1)/2*sigma_b, 
               1 - (K - 1)/2*sigma_b + (K - 1)*sigma_b, 
               by = sigma_b);
    if(direction_t[dir_count] == 1){
      Mt_diff[1:K, j] = temp
    }else{
      Mt_diff[K:1, j] = temp
    }
    dir_count = dir_count + 1
  }
  enrich_info = ifelse(direction_t > 0, "grp1_enrich", "grp0_enrich")
  A_main <- matrix(NA, nrow = n, ncol = p_remain);
  for (i in 1:n) {
    A_main[i,] <- exp(rnorm(p_remain, Mt[zt[i] + 1,], runif(n = 1, low, high)/10));
  }
  
  A_diff <- matrix(NA, nrow = n, ncol = p_diff);
  
  for (i in 1:n) {
    A_diff[i,] <- exp(rnorm(p_diff, Mt_diff[zt[i] + 1,], runif(n = 1, low, high)/10));
  }
  
  Y_main <- matrix(NA, nrow = n, ncol = p_remain);
  rownames(Y_main ) <- paste0(1:n, "-grp", zt + 1);
  Y_diff <- matrix(NA, nrow = n, ncol = p_diff);
  rownames(Y_diff ) <- paste0(1:n, "-grp", zt + 1);
  
  for (j in 1:p_remain) {
    for (i in 1:n) {
      Y_main[i, j] <- rnbinom(1, mu = s[i]*A_main[i, j], size = psi_main[j]);
    }
  }
  for (j in 1:p_diff) {
    for (i in 1:n) {
      Y_diff[i, j] <- rnbinom(1, mu = s[i]*A_diff[i, j], size = psi_diff[j]);
    }
  }
  
  H = matrix(0L, nrow = n, ncol = p_remain);
  H[sample(1:(n*p_remain), floor(n * p_remain * zero_pct))] = 1;
  Y_main[which(H == 1)] = 0;
  
  H = matrix(0L, nrow = n, ncol = p_diff);
  H[sample(1:(n*p_diff), floor(n * p_diff * zero_pct))] = 1;
  Y_diff[which(H == 1)] = 0;
  
  return(list(Y_main = Y_main, Y_diff = Y_diff, zt = zt, 
              st = s, Mt_main = Mt, Mt_diff = Mt_diff,
              A_main = A_main, A_diff = A_diff, enrich_info = enrich_info ));
}

gen_zinb_info = function(taxa_name_diff, 
                         taxa_name_all,
                         info_list){
  Y_diff = info_list$Y_diff
  Y_remain = info_list$Y_main
  
  gamma_v = rep(0, ncol(Y_diff) + ncol(Y_remain))
  da_idx = unlist(sapply(taxa_name_diff, function(x){which(x == taxa_name_all)}))
  da_idx = as.numeric( da_idx)
  gamma_v[da_idx ] = 1
  
  taxa_name_remain = taxa_name_all[-da_idx]
  Y_mat = cbind(Y_diff, Y_remain)
  colnames(Y_mat) = c(taxa_name_diff, taxa_name_remain)
  rownames(Y_mat) = rownames(Y_diff)
  
  reorder_taxa_idx = match(colnames(Y_mat), taxa_name_all)
  Y_mat = Y_mat[, order(reorder_taxa_idx)]
  return(list(Y = Y_mat, gamma_true = gamma_v))
}

