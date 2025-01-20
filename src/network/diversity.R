diversity <- edges.nocontinf[edges.nocontinf$date1 < "2017-01-01" & edges.nocontinf$date2 < "2017-01-01" & edges.nocontinf$fract_sites_IBD > -1 & !edges.nocontinf$Continf,]
diversity$kept <- diversity$fract_sites_IBD > 0.5
diversity$group <- sapply(diversity$fract_sites_IBD, FUN = function(x){
  if(x < 0.5){
    return(-1)
  }
  if(x >= 0.5 & x < 0.9){
    return(0)
  }
  if(x >= 0.9){
    return(1)
  }
})

covdensity <- data.frame(ecdf(diversity$fract_sites_IBD)(seq(0,1,0.001)))
covdensity <- cbind(row.names(covdensity), covdensity)
colnames(covdensity) <- c("IBD", "quantile")
covdensity$IBD <- seq(0,1,0.001) 

gg <- ggplot(diversity)+
  geom_hline(yintercept = c(1,10,100,1000,10000), 
             linetype = 2, color = "grey70", linewidth = 0.5)+
  geom_bar(aes(x = floor(fract_sites_IBD*20)/20,
               fill = factor(group)),
           show.legend = T, color = "black", linewidth = 0.25, width = 1/20)+
  scale_x_continuous(expand = c(0,0), n.breaks = 10)+
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(0+seq(1,10,1),
                                10+seq(10,90,10),
                                100+seq(100,900,100),
                                1000+seq(1000,9000,1000),
                                10000+seq(10000,20000,10000)),
                     labels = function(x){ifelse(log(x,10)%in%c(0:4),
                                                 round(x), "")},
                     trans = "log")+
  xlab("IBD")+
  ylab("Number of pairs")+
  scale_fill_manual(name = "",
                    values = c("white", "grey60", "black"),
                    labels = c("Unrelated", "Related", "Highly related"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 3, linewidth = 0.4),
        panel.grid.minor.y = element_blank(),
        panel.grid = element_line(color = "grey70"),
        text = element_text(size = 16),
        legend.position = "top",
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave("out/logdiversity-IBD.png", width = 10, height = 3.5)

gg <- ggplot(diversity)+
  geom_bar(aes(x = floor(fract_sites_IBD*20)/20,
               fill = factor(group)),
           show.legend = T, color = "black", linewidth = 0.25, width = 1/20)+
  geom_hline(yintercept = 500, linetype = 1, color = "red", linewidth = 0.4)+
  scale_x_continuous(expand = c(0,0), n.breaks = 10)+
  scale_y_continuous(expand = c(0,0), n.breaks = 10)+
  xlab("IBD")+
  ylab("Number of pairs")+
  scale_fill_manual(name = "",
                    values = c("white", "grey60", "black"),
                    labels = c("Unrelated", "Related", "Highly related"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 3, linewidth = 0.4),
        panel.grid.minor.y = element_blank(),
        panel.grid = element_line(color = "grey70"),
        text = element_text(size = 20),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.position.inside = c(0.5,0.85),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text = element_text(color = "black"))
gg

ggsave("out/diversity-IBD.png", width = 7, height = 3.5)

gg <- ggplot(diversity)+
  geom_bar(aes(x = floor(fract_sites_IBD*20)/20,
               y = after_stat(ifelse(count>=500, 501, count)),
               fill = factor(group)),
           show.legend = T, color = "black", linewidth = 0.25, width = 1/20)+
  geom_hline(yintercept = 500, linetype = 1, color = "red", linewidth = 0.4)+
  scale_x_continuous(expand = c(0,0), n.breaks = 10)+
  scale_y_continuous(expand = c(0,0), n.breaks = 5)+
  xlab("IBD")+
  ylab("Number of pairs")+
  scale_fill_manual(name = "",
                    values = c("white", "grey60", "black"),
                    labels = c("Unrelated", "Related", "Highly related"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 3, linewidth = 0.4),
        panel.grid.minor.y = element_blank(),
        panel.grid = element_line(color = "grey70"),
        text = element_text(size = 20),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave("out/diversity-IBD-zoom.png", width = 7, height = 3.5)

nrow(diversity[diversity$fract_sites_IBD > 0.5,])/
  nrow(diversity)
length(unique(c(diversity$ID2, diversity$ID1)))


#IBD between consensus/genotyped-only barcodes
bcodes <- read.csv("../read/out/barcodes-consensus.csv")
edges.consmol <- merge(edges, bcodes[,c("ID", "Study")],
                               by.x = "ID1", by.y = "ID", all.x = T)
edges.consmol <- merge(edges.consmol, bcodes[,c("ID", "Study")],
                               by.x = "ID2", by.y = "ID", suffixes = c("1", "2"),
                               all.x = T)
edges.consmol$Study1[is.na(edges.consmol$Study1)] <- "Consensus"
edges.consmol$Study2[is.na(edges.consmol$Study2)] <- "Consensus"
edges.consmol$Study1[edges.consmol$Study1 == "WGS"] <- "Consensus"
edges.consmol$Study2[edges.consmol$Study2 == "WGS"] <- "Consensus"
edges.consmol$commonStudy <- paste0(pmin(edges.consmol$Study1, edges.consmol$Study2), 
                                    pmax(edges.consmol$Study1,edges.consmol$Study2))

comp.consmol <- edges.consmol[edges.consmol$fract_sites_IBD>-1,] %>%
  group_by(commonStudy)%>%
  reframe(similar = as.numeric(length(which(fract_sites_IBD >= 0.5))),
          different = as.numeric(length(which(fract_sites_IBD < 0.5))))
comp.consmol <- data.frame(comp.consmol[-1], row.names = comp.consmol$commonStudy)
100*comp.consmol[,1]/(comp.consmol[,1]+comp.consmol[,2])
chisq.test(comp.consmol[-1,])
chisq.test(comp.consmol[-2,])
chisq.test(comp.consmol[-3,])

gg <- ggplot(edges.consmol[edges.consmol$fract_sites_IBD>-1,])+
  geom_bar(aes(x = commonStudy,
               fill = fract_sites_IBD >= 0.5, group = fract_sites_IBD >= 0.5), 
           position = "fill", color = NA)+
  scale_y_continuous(limits = c(0,0.02))
gg

edges.consmol.datevil <- edges.consmol[edges.consmol$fract_sites_IBD>-1,] %>%
  group_by(commonDate, commonStudy, commonVillage)%>%
  reframe(date1 = min(c(unique(date1), unique(date2))),
          date2 = max(c(unique(date1), unique(date2))),
          obs = n())
edges.consmol.datevil$commonYear <- paste0(year(edges.consmol.datevil$date1),
                                        year(edges.consmol.datevil$date2))

gg <- ggplot(edges.consmol.datevil)+
  geom_col(aes(x = commonYear,
               y = obs, fill = commonStudy, group = commonStudy), 
           color = NA, position = "fill")
gg

gg <- ggplot(edges.consmol.datevil)+
  geom_col(aes(x = commonVillage,
               y = obs, fill = commonStudy, group = commonStudy), 
           color = NA, position = "fill")
gg

