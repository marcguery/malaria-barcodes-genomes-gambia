###########################GLOBAL VARIABLES###########################
colalpha2hex <- function(colour, transparency){
  colour.rgb <- grDevices::col2rgb(colour)
  mixedcol <- colorspace::coords(colorspace::mixcolor(transparency,
                                                      colorspace::sRGB(255, 255, 255),
                                                      colorspace::sRGB(colour.rgb[1], 
                                                                       colour.rgb[2], 
                                                                       colour.rgb[3])))
  mixedcolhex <- rgb(mixedcol[1], mixedcol[2], mixedcol[3],
                     maxColorValue = 255)
}
green.light <- colalpha2hex("springgreen3", 0.125)
red.light <- colalpha2hex("red3", 0.125)
seasoncols <- c("springgreen3", "red3")
seasoncols.light <- c(green.light, red.light)
seasonnames <- c("Low", "High")

#COLOURS
mycols <- brewer.pal(3, "Greens")[-1]
mycols <- c(mycols, brewer.pal(3, "Reds")[-1])
mycols <- c(mycols[2], mycols[4], mycols[1], mycols[3], "grey60")
######################################################

#######################SAMPLING STATISTICS######################

#Sequencing
wgsdata <- read.csv("out/WGS-SNPs.csv")
wgssnps <- nchar(wgsdata$Barcode[1])
wgsdata <- wgsdata[,c(1,4,5)]
wgsdata <- wgsdata[(wgsdata$Xs) < (wgssnps-minwgssnps),]
wgsdata$date <- as.Date(paste(sub("\\S+_", "", wgsdata$Genome_ID),"01"), format = "%y%m%d")
wgsdata$type <- "genome"
barcodedata <- read.csv("out/barcodes-consensus.csv")
barcodesnps <- nchar(barcodedata$Barcode[1])
barcodedata <- barcodedata[,c(1,3,5,6)]
barcodedata <- barcodedata[(barcodedata$Ns+barcodedata$Xs) <= (barcodesnps - minbarcodesnps),]
barcodedata$type <- "barcode"
barcodedata$date <- as.Date(paste(sub("\\S+_", "", barcodedata$ID),"01"), format = "%y%m%d")
barcodedata$Study[barcodedata$Study == "WGS"] <- "Consensus"
study.sizes <- barcodedata %>%
  group_by(Study)%>%
  reframe(size = n())
study.sizes <- rbind(study.sizes,
                     data.frame("Study" = "Consensus and molecular", size = nrow(barcodedata)))
summary(barcodesnps - barcodedata$Ns - barcodedata$Xs)
summary(barcodesnps - barcodedata$Ns[barcodedata$Study == "Molecular"] - barcodedata$Xs[barcodedata$Study == "Molecular"])
summary(barcodesnps - barcodedata$Ns[barcodedata$Study == "Consensus"] - barcodedata$Xs[barcodedata$Study == "Consensus"])

gg <- ggplot(wgsdata)+
  geom_boxplot(aes(x = factor(1), y = wgssnps - Ns - Xs))+
  scale_x_discrete(labels = paste0("Genomes (", nrow(wgsdata), ")"))+
  scale_y_continuous("Number of homozygous SNPs",
                     limits = c(0,wgssnps), expand = c(0,0), n.breaks = 20)+
  theme(text=element_text(size=17),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave("out/homozygous-wgs.png", width = 5, height = 6, dpi = 400)

gg <- ggplot(barcodedata)+
  geom_boxplot(aes(x = factor(Study), y = 89 - Ns - Xs))+
  geom_boxplot(aes(x = factor("Consensus and molecular", levels = c(unique(Study), "Consensus and molecular")), y = 89 - Ns - Xs))+
  scale_x_discrete(labels = function(x){paste0(x, " barcodes (",study.sizes$size[study.sizes$Study == x], ")")})+
  scale_y_continuous("Number of homozygous SNPs",
                     limits = c(0,90), expand = c(0,0), n.breaks = 20)+
  theme(text=element_text(size=17),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave("out/homozygous-barcodes.png", width = 12, height = 6, dpi = 400)

barcodedata <- barcodedata[,-c(which(colnames(barcodedata) == "Study"))]
malariagendata <- merge(barcodedata, wgsdata, by.x = "ID", by.y = "Genome_ID", all = T)
malariagendata$date <- pmin(malariagendata$date.x, malariagendata$date.y, na.rm = T)
malariagendata$type <- pmin(malariagendata$type.x, malariagendata$type.y)
malariagendata$type[!is.na(malariagendata$type)] <- "both"
malariagendata$type[is.na(malariagendata$type)] <- pmin(malariagendata$type.x[is.na(malariagendata$type)], 
                                                        malariagendata$type.y[is.na(malariagendata$type)], na.rm = T)
malariagendata <- malariagendata[,c("ID", "date", "type")]
malariagendata$village <- substr(malariagendata$ID, 1, 1)

axisshift = 8.5
cols <- brewer.pal(9, "Set1")[2:3]
gg <- ggplot(malariagendata)+
  geom_rect(data=seasons, aes(xmin=pmax(date1, as.Date("2014-10-01")), 
                              xmax=pmin(date2, as.Date("2017-06-01")), 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.125)+
  scale_fill_manual("Transmission season", labels=seasonnames, 
                    values=seasoncols,
                    guide = guide_legend(ncol = 1, order = 1, 
                                         override.aes = c(color = "black")))+
  new_scale("fill")+
  geom_vline(data=yearspan,
             aes(xintercept = year), linewidth=1, linetype=2)+
  geom_bar_pattern(aes(x = as.Date(date), 
               group = factor(type, levels = c("barcode", "genome", "both")),
               pattern = factor(type, levels = c("barcode", "genome", "both")),
               pattern_fill = factor(type, levels = c("barcode", "genome", "both")),
               fill = type), linewidth = 0.5, color = "black", pattern_color = NA,
               pattern_density = 0.5,pattern_key_scale_factor = 0.25,
               show.legend = c(pattern_density = T, pattern = T, fill = T))+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"),
               limits = c(as.Date("2014-10-01"),as.Date("2017-06-01")))+
  ylab("Isolates")+
  xlab("Date")+
  scale_pattern_manual("Genetic data",
                       breaks = c("barcode", "genome", "both"),
                       labels=c("Barcode", "Genome", "Barcode and genome"), 
                       values = c("none", "none", "stripe"),
                       guide = guide_legend(order = 2, nrow = 2))+
  scale_fill_manual("Genetic data",
                    breaks = c("barcode", "genome", "both"),
                    labels=c("Barcode", "Genome", "Barcode and genome"), 
                    values=c(cols[2], cols[1], cols[1]),
                    guide = guide_legend(order = 2, nrow = 2,
                                         override.aes = c(linewidth = 1.5)))+
  scale_pattern_fill_manual("Genetic data",
                            breaks = c("barcode", "genome", "both"),
                            labels=c("Barcode", "Genome", "Barcode and genome"),
                            values=c(NA, NA, cols[2]),
                            guide = guide_legend(order = 2, nrow = 2))+
  scale_y_continuous(expand = c(0,0), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  annotate("text",
           label = "2014",
           x=as.Date("2014-11-17"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2014-10-05"),
           xend=as.Date("2014-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2015",
           x=as.Date("2015-06-22"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2015-01-03"),
           xend=as.Date("2015-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2016",
           x=as.Date("2016-06-22"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2016-01-03"),
           xend=as.Date("2016-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2017",
           x=as.Date("2017-03-15"),
           y=-axisshift, 
           size = 7) +
  annotate("segment",
           x=as.Date("2017-01-03"),
           xend=as.Date("2017-05-25"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  coord_cartesian(ylim=c(0,70), clip="off", expand = 0)+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        text=element_text(size=30),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(0,-20,0,0)),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'))
ggsave("out/high-quality-WGS-barcodes.png", width = 22, height = 7)

cols <- brewer.pal(4, "Set1")
gg <- ggplot(malariagendata)+
  geom_rect(data=seasons, aes(xmin=pmax(date1, as.Date("2014-10-01")), 
                              xmax=pmin(date2, as.Date("2017-06-01")), 
                              ymin=0, ymax=Inf, 
                              fill=type), alpha=0.125)+
  scale_fill_manual("Transmission season", labels=seasonnames, 
                    values=seasoncols,
                    guide = guide_legend(ncol = 1, order = 1, 
                                         override.aes = c(color = "black")))+
  new_scale("fill")+
  geom_vline(data=yearspan,
             aes(xintercept = year), linewidth=1, linetype=2)+
  geom_bar(aes(x = as.Date(date), 
                       group = factor(village, levels = rev(c("J", "K", "N", "P"))),
                       fill = village), linewidth = 0.5, color = "black")+
  scale_x_date(expand = c(0,0), 
               breaks = as.Date("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"),
               limits = c(as.Date("2014-10-01"),as.Date("2017-06-01")))+
  ylab("Isolates")+
  xlab("Date")+
  scale_fill_manual("Village",
                    values=c(cols),
                    guide = guide_legend(order = 2, nrow = 2,
                                         override.aes = c(linewidth = 0.5)))+
  scale_y_continuous(expand = c(0,0), n.breaks = 10,
                     labels = function(x){sprintf("%06s", x)})+
  annotate("text",
           label = "2014",
           x=as.Date("2014-11-17"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2014-10-05"),
           xend=as.Date("2014-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2015",
           x=as.Date("2015-06-22"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2015-01-03"),
           xend=as.Date("2015-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2016",
           x=as.Date("2016-06-22"),
           y=-axisshift, 
           size = 8) +
  annotate("segment",
           x=as.Date("2016-01-03"),
           xend=as.Date("2016-12-30"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  annotate("text",
           label = "2017",
           x=as.Date("2017-03-15"),
           y=-axisshift, 
           size = 7) +
  annotate("segment",
           x=as.Date("2017-01-03"),
           xend=as.Date("2017-05-25"),
           y=-axisshift+2.5,
           yend=-axisshift+2.5,
           color ="black") +
  coord_cartesian(ylim=c(0,70), clip="off", expand = 0)+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        text=element_text(size=30),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=margin(0,-20,0,0)),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'))
ggsave("out/high-quality-WGS-barcodes-villages.png", width = 22, height = 7)
############################################

######################FILTERING 21 BAD GENOTYPED SNPS######################

agreementplot <- ggplot()+
  geom_point(aes(x = Barcode.order, y = Agreement, color = !Barcode.order%in%badsnps),
             show.legend = F)+
  scale_x_continuous("Barcode locus", breaks = seq(5,100,5), minor_breaks = seq(1,101,1), 
                     limits = c(0.5,101.5), expand = c(0,0))+
  scale_y_continuous("Agreement", 
                     limits = c(0.2,1.01), n.breaks = 8, expand = c(0,0))+
  scale_color_manual(values = c("TRUE" = mycols[1], "FALSE" = mycols[2]))+
  theme(text=element_text(size=17),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.text = element_text(color = "black"))

#Agreement of each of the 101 SNPs between barcodes and WGS
gg <- agreementplot %+% snps.raw[!is.na(snps.raw$Agreement) & !is.nan(snps.raw$Agreement),]
ggsave("out/Agreement-WGS-molecular.png", width = 8, height = 5, dpi = 300)


sort(snps.raw.before1605[!is.na(snps.raw.before1605$Agreement) & !is.nan(snps.raw.before1605$Agreement),"Agreement"])
snps.raw.before1605[!is.na(snps.raw.before1605$Agreement) & !is.nan(snps.raw.before1605$Agreement),]
gg <- agreementplot %+% snps.raw.before1605[!is.na(snps.raw.before1605$Agreement) & !is.nan(snps.raw.before1605$Agreement),]
ggsave("out/Agreement-WGS-molecular-before1605.png", width = 8, height = 5, dpi = 300)

sort(snps.raw.after1605[!is.na(snps.raw.after1605$Agreement) & !is.nan(snps.raw.after1605$Agreement),"Agreement"])
gg <- agreementplot %+% snps.raw.after1605[!is.na(snps.raw.after1605$Agreement) & !is.nan(snps.raw.after1605$Agreement),]
ggsave("out/Agreement-WGS-molecular-after1605.png", width = 8, height = 5, dpi = 300)

#PLOT
heatmapplot <- ggplot()+
  geom_tile(aes(x=factor(SNP, levels=orderedsnps), y=factor(Sample, levels=orderedbarcodes),
                             fill=factor(Score, levels=c("-2","-1","0","1","2")),
                             color = factor(Score, levels=c("-2","-1","0","1","2"))), 
                         show.legend = F)+
  coord_fixed(expand = T)+
  scale_fill_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                    breaks = c("1", "-1", "2", "-2", "0"),
                    values=c(mycols,"grey60"),
                    drop = FALSE,
                    name = "")+
  scale_color_manual(labels=c("Match", "Mismatch", "N-N", "N-ATGC", "X-ATGCNX"),
                     breaks = c("1", "-1", "2", "-2", "0"),
                     values=c(mycols,"grey60"),
                     drop = FALSE,
                     name = "")+
  scale_x_discrete()+
  scale_y_discrete()+
  xlab("Barcode locus")+
  ylab("Sample")+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 30),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(),
        legend.position = "none")

#These 21 bad SNPs correlate well with WGS SNP before May 2016
bcodescores.before1605 <- bcodescores[,dateperbarcode <= 1605]
bcodescores.before1605.melt <- melt(bcodescores.before1605, 
                                    varnames = c("SNP", "Sample"), value.name = "Score")

gg <- heatmapplot %+% bcodescores.before1605.melt
ggsave("out/before-removal-before0516.png", plot = gg,
       width = 20, height = 14, dpi = 400)
#These 21 bad SNPs correlate badly with WGS SNP after May 2016
bcodescores.after1605 <- bcodescores[,dateperbarcode > 1605]
bcodescores.after1605.melt <- melt(bcodescores.after1605, 
                                   varnames = c("SNP", "Sample"), value.name = "Score")

gg <- heatmapplot %+% bcodescores.after1605.melt
ggsave("out/before-removal-after0516.png", plot = gg,
       width = 20, height = 32, dpi = 400)

#Before May 2016 (consensus barcodes)
bcodescores.cons.before1605 <- bcodescores.cons[,dateperbarcode <= 1605]
bcodescores.cons.before1605.melt <- melt(bcodescores.cons.before1605, 
                                             varnames = c("SNP", "Sample"), value.name = "Score")

gg <- heatmapplot %+% bcodescores.cons.before1605.melt
ggsave("out/after-removal-before0516.png", plot = gg,
       width = 20, height = 14, dpi = 400)

#After May 2016 (consensus barcodes)
bcodescores.cons.after1605 <- bcodescores.cons[,dateperbarcode > 1605]
bcodescores.cons.after1605.melt <- melt(bcodescores.cons.after1605, 
                                         varnames = c("SNP", "Sample"), value.name = "Score")

gg <- heatmapplot %+% bcodescores.cons.after1605.melt
ggsave("out/after-removal-after0516.png", plot = gg,
       width = 20, height = 32, dpi = 400)

#Stats about number of samples that recovered > 0 one of the 21 SNPs after consensus
badsnps.scores <- bcodescores[row.names(bcodescores)%in%badsnps,]
badsnps.scores.cons <- bcodescores.cons[row.names(bcodescores.cons)%in%badsnps,]
length(which(apply(badsnps.scores, 2, function(x){
  length(which(x == -1))})-
    apply(badsnps.scores.cons, 2, function(x){
      length(which(x == -1))})>0))
######################################



