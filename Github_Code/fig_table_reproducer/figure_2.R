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
taxa_for_plot = c("f__Synergistaceae", "s__Peptostreptococcus_anaerobius",
                  "g__Enterobacteriaceae_noname", "s__Anaerococcus_vaginalis")
taxa_plot_l = list()
for(i in 1:4){
  tt = taxa_for_plot[i]
  taxa_plot_l[[i]] = log_dist_compare(taxa.plot = tt)
}