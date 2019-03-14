# ============================================================================================================
# code for Figure 2 in the main text 
# Load libraries
require(ggplot2);
require(cowplot);


# ============================================================================================================
# load the averaged results for all methods in all settings #
load("/result/simulated_data_result.Rdata")

# ============================================================================================================
# (1) data simulated from the ZINB model #
MCC_plot = MCC_plot_zinb; AUC_plot = AUC_plot_zinb

MCC_plot$Case = as.numeric(MCC_plot$Case)
AUC_plot$Case = as.numeric(AUC_plot$Case)

# reorder the methdos in the plot #
ori_level = levels(MCC_plot$Method)
MCC_plot$Method = factor(MCC_plot$Method,
                         levels = ori_level[c(2,1,3,4,6,7,5)])
AUC_plot$Method = factor(AUC_plot$Method,
                         levels = ori_level[c(2,1,3,4,6,7,5)])

# generate plot #
mcc_zinb = ggplot(MCC_plot, 
                aes(x=Case, y=mcc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels= x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        legend.position="none", 
        axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+ggtitle("DM Simulation - MCC")+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab("Simulation setting")+ ylab("MCC")+
  ylim(0,1)

auc_zinb = ggplot(AUC_plot, 
                  aes(x=Case, y=auc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels=x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        legend.position="none", 
        axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab("Simulation setting")+ ylab("AUC")+
  ylim(0.5,1)+ggtitle("ZINB Simulation - AUC")

# ============================================================================================================
# (2) data simulated from the DM model #
MCC_plot = MCC_plot_dm; AUC_plot = AUC_plot_dm

MCC_plot$Case = as.numeric(MCC_plot$Case)
AUC_plot$Case = as.numeric(AUC_plot$Case)

# reorder the methdos in the plot #
ori_level = levels(MCC_plot$Method)
MCC_plot$Method = factor(MCC_plot$Method,
                         levels = ori_level[c(2,1,3,4,6,7,5)])
AUC_plot$Method = factor(AUC_plot$Method,
                         levels = ori_level[c(2,1,3,4,6,7,5)])


# generate plot #
mcc_dm = ggplot(MCC_plot, 
                aes(x=Case, y=mcc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels= x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        legend.position="none", 
        axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+ggtitle("DM Simulation - MCC")+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab("Simulation setting")+ ylab("MCC")+
  ylim(0,1)

auc_dm = ggplot(AUC_plot, 
                aes(x=Case, y=auc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels= x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        legend.position="none", 
        axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab("Simulation setting")+ ylab("AUC")+
  ylim(0.9, 1)+
  ggtitle("DM Simulation - AUC")

# generate the common legend #
legend.eg = ggplot(AUC_plot, 
                   aes(x=Case, y=auc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab("Simulation setting")+ ylab("AUC (average over 50 replicates)")+
  ylim(0.9,1)

top_row <- plot_grid(auc_dm, mcc_dm, auc_zinb,mcc_zinb,nrow = 2, ncol = 2,
                     labels = c('a', 'b','c','d'), label_size = 20, align = "v")
p.full <- plot_grid(top_row, nrow = 2, rel_heights = c(1, .1),
                     get_legend(legend.eg + 
                                  theme(legend.direction = "horizontal",
                                        legend.justification="center" ,
                                        legend.box.just = "bottom")))
plot(p.full)


