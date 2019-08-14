# code for reproducing the results for the CRC study

# ========================================================================================
# ========================================================================================
# Load library and functions
# ========================================================================================
# ========================================================================================
require(Rcpp);
source('utility/user_functions.R');

# ========================================================================================
# ========================================================================================
# load data 
# ========================================================================================
# ========================================================================================

# count matrix #
Y = read.table("Zeller_countmatrix.txt")
Y = as.matrix(Y)
n = dim(Y)[1]; p = dim(Y)[2]
# phylogenetic tree structure #
S = read.table("Zeller_phylogenetic_tree.txt", row.names = 1)
S = as.matrix(S)
# group index 
group_info = rownames(Y)
z_vector = substr(group_info, nchar(group_info), nchar(group_info))
z_vector = as.numeric(z_vector)

# ========================================================================================
# ========================================================================================
# load function
# ========================================================================================
# ========================================================================================

sourceCpp("core_zinb_x5.cpp")

# ========================================================================================
# ========================================================================================
# load algorithm settings 
# ========================================================================================
# ========================================================================================

s_ini = rep(1, n);  # Initialize the size factor for the Dirichlet process prior
iter = 20000;   # Number of iterations
b = 1; 
h = 50;
ff = 0.5;


# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm 
# ========================================================================================
# ========================================================================================
gamma_ppi_list = list()
flag_record = c()
for(cc in 1:4){
  zinbdpp_output = zinb_model_estimator(Y = Y, z = z_vector,
                                        s = s_ini, 
                                        iter = iter, DPP = TRUE,
                                        b = b, h = h,
                                        store = TRUE,
                                        aggregate = TRUE,S = S, 
                                        MRF = TRUE,G = S,
                                        count_th = 3, phi_low = 0.5)
  gamma_ppi_list[[cc]] = zinbdpp_output$gamma_ppi
  if(cc == 1){
    flag_record = zinbdpp_output$flag
  }
}

# get the averaged PPI over 4 chains #
gamma_ppi = Reduce("+", gamma_ppi_list)/length(gamma_ppi_list)
gamma_ppi = gamma_ppi[flag_record == 0]
c_gamma = BayFDR(gamma_ppi, 0.001)
# control the Bayesian FDR #
plot(gamma_ppi, type = 'h', ylim = c(0, 1), ylab = "PPI", xlab = "Taxon Index")
abline(h = c_gamma, lty = 2, col = 'red')
