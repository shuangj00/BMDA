# code for Figure 2 in the main text 

# ============================================================================================================
# ============================================================================================================
# Load library, function and data
# ============================================================================================================
# ============================================================================================================
require(ggplot2);
require(cowplot);
theme_set(theme_cowplot())
source('utility/user_functions.R');
load("data/Zeller_reporduce.Rdata")


# ========================================================================================
# ========================================================================================
# generate the plots 
# ========================================================================================
# ========================================================================================

H_est = Reduce("+", H_list)/length(H_list)
A_est = Reduce("+", A_list)/length(A_list)


# get the result for KW test #
Y_full = Y_full[, match(taxa.output.name, colnames(Y_full))]
Y_composition = t(apply(Y_full, 1, function(x){x/sum(x)}))
# p-values from the KW test #
pval_nonpara = apply(Y_composition, 2, 
                     function(x){kruskal.test(x~phenotype.z)$p.value})
pval_nonpara = p.adjust(pval_nonpara, method = "BH")


# generate the plot #
PPI_keep = apply(gammaPPI_mat,2, mean )
candi1 = log_dist_compare("s__Clostridium_symbiosum", HH = H_est, AA = A_est, z_v = phenotype.z)
candi2 = log_dist_compare("s__Eubacterium_hallii", HH = H_est, AA = A_est, z_v = phenotype.z)
candi3 = log_dist_compare("s__Lachnospiraceae_bacterium_5_1_63FAA", HH = H_est, AA = A_est, z_v = phenotype.z)
candi4 = log_dist_compare("s__Streptococcus_salivarius", HH = H_est, AA = A_est, z_v = phenotype.z)
candi5 = log_dist_compare("s__Eubacterium_ventriosum", HH = H_est, AA = A_est, z_v = phenotype.z)

