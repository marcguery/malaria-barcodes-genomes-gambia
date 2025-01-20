#07/10/2021
#!/usr/bin/env Rscript

#This script will estimate the minimal allele frequency of the majority
#base called necessary to discard the minority base called as noise

#It will comparge WGS barcodes to genotyped barcodes and count the number of
#matches (mixed to mixed) and mismacthes (mixed to not mixed) resulting from
#the WGS barcodes generated using different values of cutoffs

###########################SCRIPT###########################
nnumbers <- data.frame(ID = NA, alignment = NA)
#All these cutoffs will be tested
minpropswgs <- seq(0.5,1,0.01)
wgsbarcodes.filename.prev <- wgsbarcodes.filename
consbarcodes.filename.prev <- consbarcodes.filename
for (minprop in minpropswgs){
  wgsbarcodes.filename <- "WGS-barcodes-tmp.csv"
  consbarcodes.filename <- "barcodes-consensus-tmp.csv"
  source("WGS-barcodes.R")
  source("merge-SNP-WGS.R")
  nnumbers <- merge(nnumbers, nnumber, all = T)
}

wgsbarcodes.filename <- wgsbarcodes.filename.prev
consbarcodes.filename <- consbarcodes.filename.prev

nnumbers <- nnumbers[!is.na(nnumbers$ID),]

###########################
###########################2. DETERMINE THE CUTOFF###########################

#We chose the cutoff maximizing N matches and minimizing N mismacthes
nnumbers.melted <- melt(nnumbers)
nnumbers.summary <- nnumbers.melted %>%
  group_by(alignment, variable) %>%
  summarise(Mean = mean(value))
#nn = N-N
#nl = N-ATGC
mycols <- brewer.pal(3, "Greens")[-c(1,3)]
mycols <- c("nn" = mycols, "nl" = brewer.pal(3, "Reds")[-c(1,3)])

gg <- ggplot(data = nnumbers.summary[nnumbers.summary$alignment%in%c("nn", "nl"),], 
             aes(x = 1-as.numeric(sub("ratio", "", variable)), y = Mean))+
  geom_path(aes(group = alignment, color = alignment), linewidth = 1, show.legend = F)+
  geom_point(data = nnumbers.summary[nnumbers.summary$alignment%in%c("nn", "nl") & 
                                       nnumbers.summary$variable=="ratio0.8",],
             aes(fill = alignment),
             size = 3, pch = 21, stroke = 0)+
  xlab("Threshold of within-sample MAF to call a genomic locus heterozygous")+
  ylab("Mean number of loci per pair")+
  scale_fill_manual(name = "", labels = c("nn" = "Heterozygous call match", 
                                          "nl" = "Heterozygous-homozygous call mismatch"),
                     values = mycols,
                     guide = guide_legend(reverse = T, nrow = 2, override.aes = c("stroke" = 1, 
                                                           "size" = 7)))+
  scale_color_manual(name = "", labels = c("nn" = "Heterozygous call match", 
                                           "nl" = "Heterozygous-homozygous call mismatch"),
                    values = mycols)+
  scale_y_continuous(n.breaks = 15, expand = c(0,0), limits = c(-0.05,5.5))+
  scale_x_continuous(expand = c(0,0), n.breaks = 15, limits = c(-0.01,0.51))+
  theme(axis.text = element_text(color="black"),
        legend.position = "top",
        legend.key = element_blank(),
        text = element_text(size = 17),
        panel.background = element_blank(),
        panel.grid = element_line(linetype = 3, color = "grey70"))
gg
#This corresponds to 0.8 of minimal majority allele called
ggsave("out/wgs-barcode-cutoff.png", width = 10, height = 5, dpi = 300)

rm(minprop)
######################################################