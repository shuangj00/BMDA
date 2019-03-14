# ============================================================================================================
# code for Figure 5 in the main text 
# Load libraries
require(Rcpp);
require(ggplot2);
require(cowplot);
require(grid);
require(gridExtra);
# ============================================================================================================
# load the MCMC output by the ZINB-DPP model on the CRC study #
source("user_functions.R")
load("/result/zinb-dpp_vs_KW.Rdata.Rdata")
# H: H = {h_ij}, an n-by-p matrix with each entry representing the probability 
#    of the element being true 0

# Alpha_mat: Alpha_mat = {alpha_ij}, an n-by-p matrix with each entry being the estimated 
#            relative abundance alpha_ij

# Y_full: full observed count matrix (n-by-pp, pp = 3940) containing the taxonomic count
#         for all taxa sequenced in the study, from species to kingdom level

# flag: a binary vector for all pp = 3940 taxa in the full count matrix, with 1 representing the 
#       taxa that have less than 3 count at one or both patient groups. These taxa would be removed 
#       when fitting the model

# PPI_zind: an n-by-1 vector of the marginal posterior porbability of inclusion (PPI) given by 
#           ZINB-DPP model. Note that we only include the PPI for all taxa we considered in the model
# ============================================================================================================

# ZINB-DPP model #
PPI = PPI_zinb

# Kruskal-Wallis test #

# (a) convert to compositional data #
Y_full_red = Y_full[, flag == 0]
Y_composition = t(apply(Y_full_red, 1, function(x){x/sum(x)}))
# get group information # 
group_info = rownames(Y_full_red)
z_vector = substr(group_info, nchar(group_info), nchar(group_info))
z_vector = as.numeric(z_vector)

# (b) KW test with BH adjusted p-values
pval_nonpara = apply(Y_composition, 2, 
                     function(x){kruskal.test(x ~ z_vector)$p.value})
pval_nonpara = p.adjust(pval_nonpara, method = "BH")

# generate the plots for distribution comparison #
taxa.compare = colnames(Y_full_red)
taxa_for_plot = c("g__Campylobacter", "s__Peptostreptococcus_anaerobius",
                  "g__Enterobacteriaceae_noname", "g__Pseudoflavonifractor")
for(i in 1:4){
  tt = taxa_for_plot[i]
  dist_compare(taxa.plot = tt, taxa.output.name = taxa.compare,
               z_v = z_vector, A_matrix = Alpha_mat, H_matrix = H)
}




