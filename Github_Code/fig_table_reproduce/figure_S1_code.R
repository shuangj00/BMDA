# code for Figure S2 in the main text 

# ============================================================================================================
# ============================================================================================================
# Load library and functions
# ============================================================================================================
# ============================================================================================================
require(ggplot2);
source('user_functions.R');

# ============================================================================================================
# ============================================================================================================
# load the selected taxa by different models 
# ============================================================================================================
# ============================================================================================================
load("/result/comparison_Zeller.Rdata")


# ============================================================================================================
# ============================================================================================================
# generate the plot (the code for plotting may take a while)
# ============================================================================================================
# ============================================================================================================
Zeller_other_ggtree = list()
# frequentist methods #
for(ii in 1:5){
  print(ii); 
  if(length(top_sig_other[[ii]]$nam)==0){
    hl= NULL
    gg_tmp = tree.compare(taxa_selected = hl)
  }else{
    hl = top_sig_other[[ii]]$nam
    gg_tmp = tree.compare(taxa_selected = hl)
  }
  Zeller_other_ggtree[[ii]]= gg_tmp
}
# for ZINB-DPP #
tree_zinbdpp = tree.compare(taxa_selected = taxa.zinbdpp)

# for DM #
PPI_DM = apply(PPI_dm_matrix, 2, mean)
c_gamma = BayFDR(PPI_DM, 0.001)
taxa.name.plot = taxa.name.Zeller[flag == 0]
tree_dm = tree.compare(taxa_selected = taxa.name.plot[PPI_DM > c_gamma])

###### tree plot #######
post_plot = plot_grid(Zeller_other_ggtree[[2]],Zeller_other_ggtree[[3]],
                      Zeller_other_ggtree[[4]],Zeller_other_ggtree[[5]],
                      tree_dm, tree_zinbdpp,
                      nrow = 2, labels = c("Kruskal-Wallis","DESeq2","edgeR","metagenomeSeq","DM", "ZINB-DPP"))
plot(post_plot)
