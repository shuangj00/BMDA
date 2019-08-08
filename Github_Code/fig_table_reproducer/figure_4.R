# code for Figure 4 in the main text 

# ============================================================================================================
# ============================================================================================================
# Load library and functions
# ============================================================================================================
# ============================================================================================================
require(ggplot2);
require(cowplot);
theme_set(theme_cowplot())
require(ggbiplot);
load("data/Zeller_downstream.Rdata")

# ============================================================================================================
# ============================================================================================================
# PCA plot #
# ============================================================================================================
# ============================================================================================================
# PCA for full data #
Yfull = Y.zeller.clr
zeller.grp = ifelse(zeller.v == 1, "CRC", "nonCRC")
data.pca.full <- prcomp(Yfull, center = T,scale. = T)

g.pca.FULL <- ggbiplot::ggbiplot(data.pca.full, obs.scale = 1, var.scale = 1, 
                                 groups = zeller.grp, 
                                 ellipse = TRUE, ellipse.prob = 0.95, 
                                 circle = FALSE, var.axes = F)
g.pca.FULL <- g.pca.FULL  + scale_color_manual(values=c( "royalblue2", "maroon3"))+
  xlab("PC1 (5.5% explained variance)") + ylab("PC2 (4.3% explained variance)")
g.pca.FULL <- g.pca.FULL + theme(plot.title = element_text(hjust = 0.5),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("PCA using all taxa") 
plot(g.pca.FULL)


# PCA using species detected by ZINB-DPP #

sub.taxa.ind = colnames(Yfull) %in% zeller.sig.taxa
Ysub = Yfull[, sub.taxa.ind]
clr.sub.pca <- prcomp(Ysub, center = T,scale. = T)

g.pca.dpp <- ggbiplot::ggbiplot(clr.sub.pca, obs.scale = 1, var.scale = 1, 
                                groups = zeller.grp, 
                                ellipse = TRUE, ellipse.prob = 0.95, 
                                circle = FALSE, var.axes = F)
g.pca.dpp <- g.pca.dpp +   scale_color_manual(values=c( "royalblue2", "maroon3"))+
  xlab("PC1 (30.2% explained variance)") + ylab("PC2 (14.2% explained variance)")
g.pca.dpp <- g.pca.dpp + theme(plot.title = element_text(hjust = 0.5),legend.position="none",
                               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))  +  
  ggtitle("PCA using the 11 taxa \n identified by the ZINB-DPP model")

plot(g.pca.dpp)


# ============================================================================================================
# ============================================================================================================
# prediction result by random forest #
# ============================================================================================================
# ============================================================================================================

ggplot(prediction.acc.dataframe, aes(x=Method, y=PE*100, fill=Method)) +
  geom_boxplot() + theme(aspect.ratio = 1, 
                         axis.text.x = element_text(angle = 30, hjust = 1),
                         legend.position = 'none',
                         plot.title = element_text(hjust = 0.5),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") + ylab("Prediction Error (%)")  +
  guides(fill=guide_legend(title=" ")) + 
  ggtitle("Prediction performance using \nthe selected taxa by ZINB-DPP and other method")




