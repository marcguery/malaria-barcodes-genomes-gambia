#####################READ DATA#################################
barcode.edges <- databarcodes
barcode.edges$fract_sites_IBD[barcode.edges$N_comp_sites<minsites] <- -1
barcode.edges$fract_sites_IBD[barcode.edges$N_informative_sites<mininfsites] <- -1
availbarcodepairs <- nrow(barcode.edges[barcode.edges$fract_sites_IBD>-1,])
nrow(wgs.edges[wgs.edges$fract_sites_IBD>-1,])
######################################################

#####################MERGE#################################
colnames(wgs.edges)[which(colnames(wgs.edges)=="fract_sites_IBD")] <- "fract_sites_IBD.wgs"
colnames(barcode.edges)[which(colnames(barcode.edges)=="fract_sites_IBD")] <- "fract_sites_IBD.barcode"

edges.both <- merge(barcode.edges[,c("ID1", "ID2", "fract_sites_IBD.barcode")], pairwise.edges,
                    by = c("ID1", "ID2"))
edges.both <- merge(edges.both, wgs.edges[,c("ID1", "ID2", "fract_sites_IBD.wgs")],
               by = c("ID1", "ID2"), all.y = T)

edges.both$pairwise <- 1-edges.both$pairwise
edges.both$pairwise[edges.both$fract_sites_IBD.barcode==-1] <- -1
nrow(edges.both[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode>-1 & edges.both$fract_sites_IBD.wgs>-1,])
length(unique(c(edges.both$ID1, edges.both$ID2)))
nrow(edges.both[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode==-1 & edges.both$fract_sites_IBD.wgs>-1,])
nrow(edges.both[is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.wgs>-1,])
######################################################

#####################SENSITIVITY AND SPECIFICITY#####################
rsq <- function(x, y) summary(lm(y~x))$r.squared
lm05 <- lm(edges.both$fract_sites_IBD.barcode[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode > 0.5 & edges.both$fract_sites_IBD.wgs > 0.5] ~ 
             edges.both$fract_sites_IBD.wgs[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode > 0.5 & edges.both$fract_sites_IBD.wgs > 0.5])

summary(lm05)
rsq05 <- rsq(edges.both$fract_sites_IBD.wgs[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode > 0.5 & edges.both$fract_sites_IBD.wgs > 0.5],
             edges.both$fract_sites_IBD.barcode[!is.na(edges.both$fract_sites_IBD.barcode) & edges.both$fract_sites_IBD.barcode > 0.5 & edges.both$fract_sites_IBD.wgs > 0.5])

accuracy <- data.frame(cutoff = c(),
                       sensitivity = c(), specificity = c(),
                       precision = c(), npv = c(),
                       pos = c(), tp = c(), neg = c(), tn = c())


accuracy.IBD <- accuracy
for (ctoff in seq(0.01,1,0.01)){
  tp <- length(which(edges.both$fract_sites_IBD.barcode >= ctoff & edges.both$fract_sites_IBD.wgs >= 0.5))
  fp <- length(which(edges.both$fract_sites_IBD.barcode >= ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5))
  fn <- length(which(edges.both$fract_sites_IBD.barcode >= 0 & edges.both$fract_sites_IBD.barcode < ctoff & edges.both$fract_sites_IBD.wgs >= 0.5))
  tn <- length(which(edges.both$fract_sites_IBD.barcode >= 0 & edges.both$fract_sites_IBD.barcode < ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5))
  se <- tp/(tp+fn)
  sp <- tn/(fp+tn)
  precision <- tp/(tp+fp)
  npv <- tn/(fn+tn)
  pos <- tp + fp
  neg <- tn + fn
  propkept <- length(which(barcode.edges$fract_sites_IBD.barcode>ctoff))/nrow(barcode.edges[barcode.edges$fract_sites_IBD.barcode>0,])
  accuracy.IBD <- rbind(accuracy.IBD, data.frame(cutoff = ctoff,
                                                 sensitivity = se,
                                                 specificity = sp,
                                                 precision = precision,
                                                 npv = npv,
                                                 pos = pos,
                                                 tp = tp,
                                                 neg = neg,
                                                 tn = tn,
                                                 propkept = propkept))
}
df <- with(subset(accuracy.IBD, subset = cutoff == 0.5),
           matrix(c(tp, (pos - tp),
                    (neg - tn), tn),
                  byrow = T, nrow = 2, ncol = 2))
sum(df[,1])/sum(df)
sum(df[1,])/sum(df)
agreement.chi <- mcnemar.test(df)
pagreement <- (df[1,1] + df[2,2]) / sum(df)
prelated <- (sum(df[1,]) / sum(df)) * (sum(df[,1]) / sum(df))
punrelated <- (sum(df[2,]) / sum(df)) * (sum(df[,2]) / sum(df))
agreement.kappa <- (pagreement - (prelated + punrelated)) / (1 - (prelated + punrelated))


ctoff <- 0.5
edges.both$accuracy.IBD <- NA
edges.both$accuracy.IBD[which(edges.both$fract_sites_IBD.barcode >= ctoff & edges.both$fract_sites_IBD.wgs >= 0.5)] <- "tp"
edges.both$accuracy.IBD[which(edges.both$fract_sites_IBD.barcode >= ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5)] <- "fp"
edges.both$accuracy.IBD[which(edges.both$fract_sites_IBD.barcode >= 0 & edges.both$fract_sites_IBD.barcode < ctoff & edges.both$fract_sites_IBD.wgs >= 0.5)] <- "fn"
edges.both$accuracy.IBD[which(edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.barcode < ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5)] <- "tn"
edges.both$accuracy.IBD <- factor(edges.both$accuracy.IBD, levels = c("tp", "fp", "tn", "fn"),
                             labels = c("True Positive", "False Positive",
                                        "True Negative", "False Negative"))
rmse.ibd <- with(edges.both[edges.both$fract_sites_IBD.barcode > -1,],
             sqrt(sum((fract_sites_IBD.barcode - fract_sites_IBD.wgs)^2, na.rm = T)/length(fract_sites_IBD.barcode)))
rmse.ibd.sup05 <- with(edges.both[edges.both$fract_sites_IBD.barcode >= 0.5,],
                   sqrt(sum((fract_sites_IBD.barcode - fract_sites_IBD.wgs)^2, na.rm = T)/length(fract_sites_IBD.barcode)))
rmse.ibd.inf05 <- with(edges.both[edges.both$fract_sites_IBD.barcode >=0 & edges.both$fract_sites_IBD.barcode < 0.5,],
                   sqrt(sum((fract_sites_IBD.barcode - fract_sites_IBD.wgs)^2, na.rm = T)/length(fract_sites_IBD.barcode)))
agreement.df <- data.frame(minsites = minsites,
                      mininfsites = mininfsites,
                      bcodepairs = availbarcodepairs,
                      rsquared05 = round(rsq05, 4),
                      sens05 = round(accuracy.IBD$sensitivity[accuracy.IBD$cutoff == 0.5], 4),
                      spec05 = round(accuracy.IBD$specificity[accuracy.IBD$cutoff == 0.5], 4),
                      prec05 = round(accuracy.IBD$precision[accuracy.IBD$cutoff == 0.5], 4),
                      agreement.kappa05 = round(agreement.kappa, 4),
                      rmse = round(rmse.ibd, 4),
                      rmsesup05 = round(rmse.ibd.sup05, 4),
                      rmseinf05 = round(rmse.ibd.inf05, 4))

accuracy.pairwise <- accuracy
edges.both.sub <- edges.both[!is.na(edges.both$pairwise),]
for (ctoff in seq(0.01,1,0.01)){
  tp <- length(which(edges.both.sub$pairwise >= ctoff & edges.both.sub$fract_sites_IBD.wgs >= 0.5))
  fp <- length(which(edges.both.sub$pairwise >= ctoff & edges.both.sub$fract_sites_IBD.wgs >= 0 & edges.both.sub$fract_sites_IBD.wgs < 0.5))
  fn <- length(which(edges.both.sub$pairwise >= 0 & edges.both.sub$pairwise < ctoff & edges.both.sub$fract_sites_IBD.wgs >= 0.5))
  tn <- length(which(edges.both.sub$pairwise>=0 & edges.both.sub$pairwise < ctoff & edges.both.sub$fract_sites_IBD.wgs >= 0 & edges.both.sub$fract_sites_IBD.wgs < 0.5))
  se <- tp/(tp+fn)
  sp <- tn/(fp+tn)
  precision <- tp/(tp+fp)
  npv <- tn/(fn+tn)
  pos <- tp + fp
  neg <- tn + fn
  propkept <- length(which(pairwise.edges$pairwise>ctoff))/nrow(pairwise.edges)
  accuracy.pairwise <- rbind(accuracy.pairwise, data.frame(cutoff = ctoff,
                                                  sensitivity = se,
                                                  specificity = sp,
                                                  precision = precision,
                                                  npv = npv,
                                                  pos = pos,
                                                  tp = tp,
                                                  neg = neg,
                                                  tn = tn,
                                                  propkept = propkept))
}

ctoff <- 0.84
edges.both$accuracy.pairwise <- NA
edges.both$accuracy.pairwise[edges.both$pairwise >= ctoff & edges.both$fract_sites_IBD.wgs >= 0.5] <- "tp"
edges.both$accuracy.pairwise[edges.both$pairwise >= ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5] <- "fp"
edges.both$accuracy.pairwise[edges.both$pairwise >= 0 & edges.both$pairwise < ctoff & edges.both$fract_sites_IBD.wgs >= 0.5] <- "fn"
edges.both$accuracy.pairwise[edges.both$pairwise>=0 & edges.both$pairwise < ctoff & edges.both$fract_sites_IBD.wgs >= 0 & edges.both$fract_sites_IBD.wgs < 0.5] <- "tn"
edges.both$accuracy.pairwise <- factor(edges.both$accuracy.pairwise, levels = c("tp", "fp", "tn", "fn"),
                             labels = c("True Positive", "False Positive",
                                        "True Negative", "False Negative"))
##########################################

#####################DEVIATION#####################
edges.both$deviation.IBD <- edges.both$fract_sites_IBD.barcode-edges.both$fract_sites_IBD.wgs
edges.both$deviation.IBD[edges.both$fract_sites_IBD.barcode==-1 | edges.both$fract_sites_IBD.wgs==-1] <- NA
summary(abs(edges.both$deviation.IBD))
length(which(abs(edges.both$deviation.IBD)<=0.1))/length(which(abs(edges.both$deviation.IBD)>=0))
length(which(edges.both$fract_sites_IBD.barcode>=0.5 & edges.both$fract_sites_IBD.wgs<0.4))
length(which(edges.both$fract_sites_IBD.barcode<0.5 & edges.both$fract_sites_IBD.barcode>-1))
length(which(edges.both$fract_sites_IBD.barcode<0.5 & edges.both$fract_sites_IBD.barcode>-1 & edges.both$fract_sites_IBD.wgs > 0.5))
##########################################

#####################PLOT#################################

dir.create(paste0("out/IBDcorr/", minsites, "minsites"), showWarnings = F)
accuracy.IBD$type <- "IBD"
accuracy.pairwise$type <- "pairwise"
accuracy.both <- rbind(accuracy.IBD,accuracy.pairwise)
accuracy.both <- melt(accuracy.both[,c("type", "cutoff",
                                       "sensitivity", "specificity",
                                       "precision")], id.vars = c("type","cutoff"))
gg <- ggplot(data = accuracy.both)+
  geom_line(aes(x = cutoff, y = value, linetype = type, color = variable))+
  scale_color_brewer("Classification metric",
                     palette = "Set1",
                     guide = guide_legend(nrow = 3, order = 2))+
  scale_linetype_manual("Barcode comparison method",
                        values = c("IBD" = 1, "pairwise" = 2),
                        guide = guide_legend(nrow = 3, order = 1))+
  scale_x_continuous("Value of comparison method corresponding to 0.5 in genomic IBD",
                     n.breaks = 15, expand = c(0.01,0.01))+
  scale_y_continuous("Value of classification metric",
                     n.breaks = 15, expand = c(0.01,0.01))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-vs-pairwise-accuracy.png"),
       width = 11, height = 5.5, dpi = 300)

gg <- ggplot(edges.both[edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = fract_sites_IBD.barcode,
                 fill = accuracy.IBD),
             show.legend = T, pch = 21)+
  scale_fill_manual(name = "Classification",
                     values = c("True Positive" = "green4", "True Negative" = "red4",
                                "False Positive" = "red", "False Negative" = "green2"),
                     labels = c("True Positive" = "Truly related", 
                                "True Negative" = "Truly unrelated",
                                "False Positive" = "Falsely related", 
                                "False Negative" = "Falsely unrelated"),
                     guide = guide_legend(ncol = 2, override.aes = c(size = 5)))+
  geom_segment(data=data.frame(1),
               x = 0, y = 0,
               xend = Inf, yend = Inf, linetype = 2, alpha = 0.7, linewidth = 1, color = "black")+
  xlab("Genome-IBD")+
  ylab("Barcode-IBD")+
  scale_x_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top",
        legend.key = element_blank())
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-correlation.png"),
       width = 10, height = 7, dpi = 500)

gg <- ggplot(edges.both[!is.na(edges.both$accuracy.IBD) & edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = fract_sites_IBD.barcode,
                 fill = factor(round(abs(deviation.IBD)*5)/5)),
             show.legend = T, pch = 21)+
  scale_fill_brewer(name = "Difference in IBD value",
                     palette = "RdYlGn",
                     labels = c("< 0.1", "< 0.3",
                                "< 0.5", "< 0.7", "< 0.9"),
                     direction = -1,
                     guide = guide_legend(nrow = 2, override.aes = c(size = 5)))+
  geom_segment(data=data.frame(1),
               x = 0, y = 0,
               xend = Inf, yend = Inf, linetype = 2, alpha = 0.7, linewidth = 1, color = "black")+
  xlab("Genome-IBD")+
  ylab("Barcode-IBD")+
  scale_x_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top",
        legend.key = element_blank())
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-deviation.png"),
       width = 10, height = 7, dpi = 500)

gg <- ggplot(edges.both[!is.na(edges.both$accuracy.IBD) & edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = fract_sites_IBD.barcode), fill = "black", pch = 21)+
  geom_segment(data=data.frame(1),
               x = 0, y = 0,
               xend = Inf, yend = Inf, linetype = 2, alpha = 0.7, linewidth = 1, color = "black")+
  xlab("Genome-IBD")+
  ylab("Barcode-IBD")+
  scale_x_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        legend.key = element_blank())
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-bw.png"),
       width = 10, height = 7, dpi = 500)

gg <- ggplot(edges.both[!is.na(edges.both$accuracy.IBD) & edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_histogram(aes(x = abs(deviation.IBD)),
                 binwidth = 0.025)+
  scale_x_continuous("Deviation to WGS IBD", limits = c(-0.05, 1.05), n.breaks = 10, expand = c(0,0))+
  scale_y_continuous("Number of pairs", n.breaks = 10, expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-deviationprop.png"),
       width = 10, height = 7, dpi = 500)

gg <- ggplot(edges.both[!is.na(edges.both$accuracy.IBD) & edges.both$fract_sites_IBD.barcode>=0.5 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_histogram(aes(x = abs(deviation.IBD)),
                 binwidth = 0.05)+
  scale_x_continuous("Deviation to WGS IBD", limits = c(-0.05, 1.05), n.breaks = 10, expand = c(0,0))+
  scale_y_continuous("Number of pairs", n.breaks = 10, expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-deviationprop0.5.png"),
       width = 10, height = 7, dpi = 500)

gg <- ggplot(edges.both[!is.na(edges.both$accuracy.IBD) & edges.both$fract_sites_IBD.barcode>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = fract_sites_IBD.barcode,
                 fill = factor(round(abs(deviation.IBD)*5)/5)),
             show.legend = T, pch = 21)+
  scale_fill_brewer(name = "Deviation",
                    palette = "RdYlGn",
                    labels = function(x){
                      paste0("<= ", as.character(as.numeric(x)+0.5/5))
                    },
                    direction = -1,
                    guide = guide_legend(ncol = 2))+
  xlab("Genome-IBD")+
  ylab("Barcode-IBD")+
  scale_x_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1.01), n.breaks = 10, expand = c(0.01,0))+
  theme(panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        text = element_text(size = 17),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
gg
ggsave(paste0("out/IBDcorr/", minsites, "minsites/ibd-deviation.png"),
       width = 10, height = 7, dpi = 500)


gg <- ggplot(edges.both[edges.both$pairwise>=0 & edges.both$fract_sites_IBD.wgs >= 0,])+
  geom_point(aes(x = fract_sites_IBD.wgs, y = pairwise,
                 color = accuracy.pairwise),
             show.legend = T)+
  scale_color_manual(name = "Classification",
                     values = c("True Positive" = "green4", "True Negative" = "red4",
                                "False Positive" = "green2", "False Negative" = "red"),
                     labels = c("Truly related", "Truly unrelated",
                                "Falsely related", "Falsely unrelated"))+
  geom_segment(data=data.frame(1),
               x = 0, y = 0,
               xend = Inf, yend = Inf, linetype = 2, alpha = 0.8, linewidth = 0.75)+
  geom_hline(yintercept = 0.84, color = "blue", linetype = 3)+
  xlab("Genome-IBD")+
  ylab("Barcode-pairwise")+
  scale_x_continuous(limits = c(0,1), n.breaks = 10, expand = c(0.01,0))+
  scale_y_continuous(limits = c(0,1), n.breaks = 10, expand = c(0.01,0))
gg

ggsave(paste0("out/IBDcorr/", minsites, "minsites/pairwise-correlation.png"),
       width = 10, height = 7, dpi = 500)
############################

