# ***README***
# BMDA
# The following script is used to fit the zero-inflated negative binomial model with Dirichlet process prior
# proposed in the submitted manuscript titled "Bayesian Modeling of Microbiome Data for Differential Abundance
# Analysis".

# The following code contains data loading, model fitting and results visualizing
# ***END***

# ========================================================================================
# ========================================================================================
# Load libraries
# ========================================================================================
# ========================================================================================

require(Rcpp);
require(ggplot2);

# ========================================================================================
# ========================================================================================
# Load functions
# ========================================================================================
# ========================================================================================
source('utility/user_functions.R');
Rcpp::sourceCpp('utility/core_zinb_x5.cpp');

# ========================================================================================
# ========================================================================================
# Load data
# ========================================================================================
# ========================================================================================
# count matrix #
Y = read.table("data/simulation_single.txt")
Y = as.matrix(Y)
n = dim(Y)[1]; p = dim(Y)[2]

# group index and ture discriminating features
group_info = rownames(Y)
z_vector = substr(group_info, nchar(group_info), nchar(group_info))
z_vector = as.numeric(z_vector) - 1
gamma_true = rep(0, dim(Y)[2])
gamma_true[grep(".TP", colnames(Y))] = 1



# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
s_ini = rep(1, n);  # Initialize the size factor for the Dirichlet process prior
iter = 20000;   # Number of iterations
S = matrix(1, 1, 1); G = matrix(1, 1, 1)  # no phylogenetic tree structure



# ========================================================================================
# ========================================================================================
# Load hyperparameters
# ========================================================================================
# ========================================================================================
b = 1; h = 10;



# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm 
# ========================================================================================
# ========================================================================================

# random seed for MCMC #
set.seed(198642)

# run zinb-dpp model #
zinbdpp_output = zinb_model_estimator(Y = Y, z = z_vector, s = s_ini, 
                                      iter = iter, DPP = TRUE, 
                                      b = b, h = h,
                                      S = S, G = G,store = TRUE,
                                      aggregate = FALSE, MRF = F)


# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
zinbdpp_PPI = zinbdpp_output$gamma_ppi
c_gamma = BayFDR(PPI = zinbdpp_PPI, alpha = 0.05)
ppi_df = data.frame(x = seq(1, length(zinbdpp_PPI)),
                    y = zinbdpp_PPI)
ppi_point_df = data.frame(x = which(gamma_true == 1),
                          y = zinbdpp_PPI[which(gamma_true == 1)])
# (1) posterior probability of inclusion # 
ggplot(ppi_df, 
       aes(x = x,xend = x, y = 0, yend = y) ) + 
  geom_segment() +
  geom_point(data = ppi_point_df, aes(x, y), color = 'red')+
  geom_hline(yintercept = c_gamma, linetype = "dashed", color = 'red') + 
  xlab("Feature index") + ylab("Marginal PPI")

# true positive #
sum(zinbdpp_PPI > c_gamma & gamma_true == 1)
# true negatige #
sum(zinbdpp_PPI <= c_gamma & gamma_true == 0)
# false positive #
sum(zinbdpp_PPI > c_gamma & gamma_true == 0)
# false negative #
sum(zinbdpp_PPI <= c_gamma & gamma_true == 1)


