# code for Figure 1 in the main text 

# ============================================================================================================
# ============================================================================================================
# Load library and functions
# ============================================================================================================
# ============================================================================================================
require(ggplot2);
require(cowplot);
source('user_functions.R');

# ============================================================================================================
# ============================================================================================================
# load the simulated data and the true information #
# ============================================================================================================
# ============================================================================================================
# count matrix #
Y = read.table("simulation_single.txt")
Y = as.matrix(Y)
n = dim(Y)[1]; p = dim(Y)[2]
# simulated size factor #
s_true = read.table("true_size_factor.txt")
s_true = as.vector(s_true)
# group index and ture discriminating features
group_info = rownames(Y)
z_vector = substr(group_info, nchar(group_info), nchar(group_info))
z_vector = as.numeric(z_vector) - 1
gamma_true = rep(0, dim(Y)[2])
gamma_true[grep(".TP", colnames(Y))] = 1


# ============================================================================================================
# ============================================================================================================
# load the MCMC results #
# ============================================================================================================
# ============================================================================================================
load("/data/simulation_example_zinbK=2sig=2.Rdata")


# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
zinbdpp_PPI = zinbdpp_output$gamma_ppi
cgamma = BayFDR(zinbdpp_PPI, 0.05)
ppi_df = data.frame(x = seq(1, length(zinbdpp_PPI)),
                    y = zinbdpp_PPI)
ppi_point_df = data.frame(x = which(gamma_true == 1),
                          y = zinbdpp_PPI[which(gamma_true == 1)])
# (1) posterior probability of inclusion # 
ggplot(ppi_df, 
       aes(x = x,xend = x, y = 0, yend = y) ) + 
  geom_segment() +
  geom_point(data = ppi_point_df, aes(x, y), color = 'red')+
  geom_hline(yintercept = cgamma, linetype = "dashed", color = 'red') + 
  xlab("Feature index") + ylab("Marginal PPI")

# true positive #
sum(zinbdpp_PPI > cgamma & gamma_true == 1)
# true negatige #
sum(zinbdpp_PPI <= cgamma & gamma_true == 0)
# false positive #
sum(zinbdpp_PPI > cgamma & gamma_true == 0)
# false negative #
sum(zinbdpp_PPI <= cgamma & gamma_true == 1)

# (2) size factor estimation with Dirichlet process prior #
s_mcmc = zinbdpp_output$s
si_df_dpp = data.frame(st = s_true, 
                       dpp = apply(s_mcmc, 2, mean),
                       low = apply(s_mcmc, 2, quantile, 0.025),
                       up = apply(s_mcmc, 2, quantile, 0.975))
colnames(si_df_dpp)[1] = "st"

# (3) size factor estimations given by other methods #
sizefactor_mat = matrix(0, nrow = 5, ncol = n)
methods = c("TSS", "CSS", "TMM", "RLE", "Q75")
for(ii in 1:5){
  mt = methods[ii]
  sizefactor_mat[ii, ] = size_factor_estimator(Y = Y, method = mt )
}
rownames(sizefactor_mat) = methods
# normalize the size factor estimation #
size_trans = apply(sizefactor_mat, 1,size_normalization )

si_est_df = data.frame(s_t = s_true, 
                       tss = size_trans[, 1], 
                       css = size_trans[, 2],
                       tmm = size_trans[, 3],
                       q75 = size_trans[, 5])
colnames(si_est_df)[1] = "s_t"

sub_dpp = ggplot(si_df_dpp, aes(st, dpp)) +
  geom_errorbar(aes(ymin = low,ymax = up),width = 0.2)+
  geom_point(colour = "black", size = 2.5, shape=1) +
  ylim(0, 5)+xlim(0, 5)+xlab("True")+ylab("Estimated")+ggtitle("DPP") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
sub_tss = ggplot(si_est_df, aes(s_t, tss)) + geom_point(colour = "black", size = 2.5, shape=1) +
  ylim(0, 5)+xlim(0, 5)+xlab("True")+ylab("Estimated")+ggtitle("TSS") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
sub_q75 = ggplot(si_est_df, aes(s_t, q75)) + geom_point(colour = "black", size = 2.5, shape=1) +
  ylim(0, 5)+xlim(0, 5)+xlab("True")+ylab("Estimated")+ggtitle("Q75") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
sub_tmm = ggplot(si_est_df, aes(s_t, tmm)) + geom_point(colour = "black", size = 2.5, shape=1) +
  ylim(0, 5)+xlim(0, 5)+xlab("True")+ylab("Estimated")+ggtitle("TMM") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
sub_css = ggplot(si_est_df, aes(s_t, css)) +geom_point(colour = "black", size = 2.5, shape=1) +
  ylim(0, 5)+xlim(0, 5)+xlab("True")+ylab("Estimated")+ggtitle("CSS") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

top_row = plot_grid(sub_tss, sub_q75, sub_tmm,align = 'h', ncol = 3 )
bottom_row = plot_grid(sub_css,sub_dpp,align = 'h', ncol = 3 )
plot_grid(top_row, bottom_row, ncol = 1)
