# ============================================================================================================
# code for Figure 3 in the main text 
# Load libraries
require(ggplot2);



# ============================================================================================================
# load the averaged results for all methods in all settings #
load("/result/synthetic_Skin_plot.Rdata")


# ============================================================================================================
MCC_plot$Case = as.numeric(MCC_plot$Case)
AUC_plot$Case = as.numeric(AUC_plot$Case)

# reorder the methdos in the plot #
ori_level = levels(MCC_plot$Method)
MCC_plot$Method = factor(MCC_plot$Method,
                         levels = ori_level[c(1,7,2,3,4,6,5)])
AUC_plot$Method = factor(AUC_plot$Method,
                         levels = ori_level[c(1,7,2,3,4,6,5)])

# AUC plot #
ggplot(AUC_plot, aes(x=Case, y=auc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels=x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab(" ")+ ylab("AUC")+
  ylim(0.5,1)+ggtitle("Synthetic Data - AUC")

# MCC plot #
ggplot(MCC_plot,  aes(x=Case, y=mcc_mean, color=Method, shape=Method)) + 
  geom_line(linetype = 1) + 
  geom_point(size = 4) + 
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels= x_tick)+
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.1))
  )+ggtitle("Synthetic Data - MCC")+
  scale_shape_manual(values=c(16, 24, 3,4,5,13,9))+
  xlab(" ")+ ylab("MCC")+
  ylim(0,1)
