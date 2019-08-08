# code for Figure 1 in the main text 

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

# Bayes FDR #
sig_level = 0.01

# select the differentially abundant taxa #
PPI_keep = apply(gammaPPI_mat,2, mean )
c_gamma = BayFDR(PPI_keep, alpha = sig_level)
taxa.sel.names = taxa.output.name[PPI_keep > c_gamma]
taxa.sele.num = length(taxa.output.name[PPI_keep > c_gamma])
taxa.sele.num # we selected 33 differential taxa #

# 1.95% credible interval #
# get CI plot #
CI_res = CIalpha.ggplot(PPI = PPI_keep,  
                        tax.name.full = taxa.output.name, 
                        alpha.mat = lalpha_mat, 
                        alpha.sig = sig_level,
                        inputname.1 = "non-CRC", inputname.2="CRC", 
                        shift_right = 130, size_y = 1.5,
                        shift_ang = 30, taxa.font = 1)
plot(CI_res)

# 2.marginal posterior of inclusion #
PPI_res = gamma.PPI.draw.pick(gamma.PPI = PPI_keep, 
                    tax.name.full = taxa.output.name,
                    pick.name = c("s__Fusobacterium_nucleatum", "g__Enterobacteriaceae_noname"),
                    alpha.sig = sig_level)
plot(PPI_res)


# 3.cladogram #
taxa.all = taxa.output.name
taxa.sig = taxa.all[PPI_keep > c_gamma]

taxa.sig.control = taxa.sig[c( 2, 7, 18, 21 )]
taxa.sig.case = taxa.sig[-c(2, 7, 18, 21)]

tree.check(controln = taxa.sig.control, casen = taxa.sig.case, hl_size = 2)


