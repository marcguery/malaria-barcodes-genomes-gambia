###########################GLOBAL VARIABLES###########################
seasoncols <- c("springgreen3", "red3")
seasonnames <- c("Low", "High")
######################################################

######################READING######################
fws <- read.csv("out/fws-27577-5-3750-persample.csv")
fws.af <- read.table("out/maf.tsv", h = F)
nnum <- read.csv("out/Ns-prop.csv")
############################################

######################MERGING######################
coi <- merge(nnum, fws, by.x = "ID", by.y = "Sample", all = TRUE)

colnames(fws.af) <- c("name", "AF")
############################################

######################METHOD COMPARISON######################

mycols <- brewer.pal(3, "Greens")[3]
mycols <- c(mycols, brewer.pal(3, "Reds")[3])

gg <- ggplot(coi[!is.na(coi$prop.WGS) & !is.na(coi$Fws),])+
  geom_blank(aes(x = prop.WGS,
                 y = Fws))+
  # geom_hline(yintercept = 0.95, color = "blue", linetype = 2)+
  # geom_vline(xintercept = 0.005, color = "blue", linetype = 2)+
  geom_point(aes(x = prop.WGS,
                 y = Fws,
                 fill = (Fws > 0.95 & prop.WGS <= 0.005 | 
                           Fws < 0.95 & prop.WGS > 0.005)), 
             alpha = 0.9, 
             pch = 21, color = "black", show.legend = F)+
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("black", "black"))+
  scale_y_continuous(n.breaks = 10, limits = c(0.3,1.01),
                     expand = c(0,0))+
  scale_x_continuous(breaks = c(0.005,seq(0.05,0.5,0.05)),
                     labels = function(x){x*100},
                     expand = c(0.001,0.001))+
  xlab("Percentage of heterozygous loci (Genomes)")+
  ylab("Fws (Genomes)")+
  theme(panel.grid.minor.x = element_blank(),
        text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        axis.text = element_text(color = "black"),
        axis.ticks.x = element_blank())
gg
ggsave(paste0("out/NsGenomes-Fws.", imgfmt),
       width = 8, height = 5, units = "in", dpi = 500)

gg <- ggplot(coi[!is.na(coi$prop.SNP) & !is.na(coi$Fws),])+
  geom_blank(aes(x = prop.SNP,
                 y = Fws))+
  # geom_hline(yintercept = 0.95, color = "blue", linetype = 2)+
  # geom_vline(xintercept = 0.005, color = "blue", linetype = 2)+
  geom_point(aes(x = prop.SNP,
                 y = Fws,
                 fill = (Fws > 0.95 & prop.SNP <= 0.005 | 
                           Fws < 0.95 & prop.SNP > 0.005)), 
             alpha = 0.9, 
             pch = 21, color = "black", show.legend = F)+
  scale_fill_manual(breaks = c("TRUE", "FALSE"),
                    values = c("black", "black"))+
  scale_y_continuous(n.breaks = 10, limits = c(0.3,1.01),
                     expand = c(0,0))+
  scale_x_continuous(breaks = c(0.005,seq(0.05,0.5,0.05)),
                     labels = function(x){x*100},
                     expand = c(0.01,0.01))+
  xlab("Percentage of heterozygous loci (Barcodes)")+
  ylab("Fws (Genomes)")+
  theme(text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks.x = element_blank())
gg
ggsave(paste0("out/NsBarcodes-Fws.", imgfmt), 
       width = 8, height = 5, units = "in", dpi = 500)

coi$date <- paste0(substr(coi$ID, 10, 13),"01")

coiprop <- coi%>%
  group_by(date)%>%
  summarise(fwscoi = round(length(which(Fws<0.95))/length(which(!is.na(Fws))),2),
            fwsnum = length(which(!is.na(Fws))),
            nsbarcodecoi = round(length(which(prop.SNP>0.005))/length(which(!is.na(prop.SNP))),2),
            nsbarcodenum = length(which(!is.na(prop.SNP))),
            nswgscoi = round(length(which(prop.WGS>0.005))/length(which(!is.na(prop.WGS))),2),
            nswgsnum = length(which(!is.na(prop.WGS)))
  )

coiprop.melted <- melt(coiprop, id.vars = "date", measure.vars = seq(2,ncol(coiprop),2))
coiprop.melted2 <- melt(coiprop, id.vars = "date", measure.vars = seq(3,ncol(coiprop),2))
coiprop.melted$variable <- sub("coi$", "", coiprop.melted$variable)
coiprop.melted2$variable <- sub("num$", "", coiprop.melted2$variable)
colnames(coiprop.melted) <- c("Date", "Method", "COI")
colnames(coiprop.melted2) <- c("Date", "Method", "Numobs")
coiprop.melted <- merge(coiprop.melted, coiprop.melted2)
coiprop.melted <- coiprop.melted[with(coiprop.melted, order(Method, Date)),]
rm(coiprop.melted2)

coiprop.melted$Date <- as.Date(coiprop.melted$Date, format = "%y%m%d")

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(as.numeric(unique(format(coiprop.melted$Date,"%Y"))))
wetseason <- as.Date(as.character(as.Date(c(paste0(years[1], "-01-01"),
                                            unlist(lapply(years, function(x){
                                              paste(x, wetseason, "15", sep="-")
                                            })),
                                            paste0(years[length(years)], "-12-31")), "%Y-%m-%d")), "%Y-%m-%d")

seasons <- data.frame(wetseason[-length(wetseason)], wetseason[-1])
seasons$type <- rep(c("dry", "wet"), length.out=nrow(seasons))
seasons$cycle <- rep(c(-1,-1), length.out=nrow(seasons))
seasons$cycle <- ifelse(rep(seasons$type[1]=="wet", length.out=nrow(seasons)),
                        seasons$cycle+floor((0:(nrow(seasons)-1))/2)+1,
                        seasons$cycle+floor((1:nrow(seasons))/2)+1)
colnames(seasons) <- c("date1", "date2", "type", "cycle")
mindate <- min(coiprop.melted$Date)
maxdate <- max(coiprop.melted$Date)
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]
seasons$date2[which.max(seasons$date2)] <- maxdate+15
seasons$date1[which.min(seasons$date1)] <- mindate-15

##STANDARD ERROR####
coiprop.melted <- coiprop.melted[!is.nan(coiprop.melted$COI),]
coiprop.melted$err <- 1.96*sqrt(coiprop.melted$COI*(1-coiprop.melted$COI)/
                             coiprop.melted$Numobs)
coiprop.melted$errwilson <- binom.confint(round(coiprop.melted$COI*coiprop.melted$Numobs), 
                                          coiprop.melted$Numobs, method=c("wilson"))
########

axisshift <- 0.11
spread <- seq(-3,3, length.out = length(unique(coiprop.melted$Method)))

length(which(coiprop.melted$Numobs<5))
summary(coiprop.melted$Numobs[coiprop.melted$Numobs>=5])
cols <- brewer.pal(3, "Set1")

gg <- ggplot(data=coiprop.melted[coiprop.melted$Numobs>=5,])+
  geom_rect(data=seasons, aes(ymin=-Inf, ymax=Inf, 
                              xmin=date1, xmax=date2, 
                              fill=toupper(type)), alpha=0.125)+
  scale_fill_manual("Transmission season",
                    labels = seasonnames,
                    values=seasoncols,
                    guide = guide_legend(nrow = 2, order = 1, 
                                         override.aes = c(color = "black")))+
  new_scale("fill")+
  geom_segment(aes(x = Date+spread[as.numeric(factor(Method))], 
                   y = errwilson$lower+0.001,
                   xend = Date+spread[as.numeric(factor(Method))], 
                   yend = errwilson$upper-0.001,
                   color = Method,
                   group = Method),
               alpha = 0.3, show.legend = F, linewidth = 0.8)+
  geom_segment(aes(x = Date+spread[as.numeric(factor(Method))]-
                     2.5, 
                   y = errwilson$upper,
                   xend = Date+spread[as.numeric(factor(Method))]+
                     2.5, 
                   yend = errwilson$upper,
                   color = Method,
                   group = Method),
               alpha = 0.3, show.legend = F)+
  geom_segment(aes(x = Date+spread[as.numeric(factor(Method))]-
                     2.5, 
                   y = errwilson$lower,
                   xend = Date+spread[as.numeric(factor(Method))]+
                     2.5, 
                   yend = errwilson$lower,
                   color = Method,
                   group = Method),
               alpha = 0.3, show.legend = F)+
  geom_path(aes(x = Date+spread[as.numeric(factor(Method))], 
                y = COI, color = Method, group = Method), 
            show.legend = F)+
  geom_point(aes(x = Date+spread[as.numeric(factor(Method))], 
                 y = COI, fill = Method, pch = Method),
             stroke = 0.25, size = 2.25, color = "black",
             show.legend = c(pch = T, fill = T))+
  scale_shape_manual("COI",
                     breaks = c( "fws", 
                                 "nswgs", "nsbarcode"),
                     labels = c("Fws (Genomes)",
                                "Heterozygous loci (Genomes)", "Heterozygous loci (Barcodes)"),
                     values = c(21, 24, 21, 21, 24),
                     guide = guide_legend(order = 2, nrow = 2))+
  scale_color_manual("COI", values = rev(cols),
                     breaks = c( "fws", 
                                 "nswgs", "nsbarcode"),
                     labels = c("Fws (Genomes)",
                                "Heterozygous loci (Genomes)", "Heterozygous loci (Barcodes)"),
                    guide = guide_legend(nrow = 2, order = 2))+
  scale_fill_manual("COI", values = rev(cols),
                    breaks = c( "fws", 
                                "nswgs", "nsbarcode"),
                    labels = c("Fws (Genomes)",
                               "Heterozygous loci (Genomes)", "Heterozygous loci (Barcodes)"),
                     guide = guide_legend(nrow = 2, order = 2, 
                                          override.aes = c(size = 5, 
                                                           stroke = 1)))+
  scale_x_date("", expand = c(0,0), date_breaks = "2 months", date_minor_breaks = "1 month",
               date_labels = "%b")+
  scale_y_continuous(n.breaks = 10)+
  ylab("Proportion of mixed infections")+
  annotate("text",
           label = "2014",
           x=as.Date("2014-12-7"),
           y=-axisshift, 
           size = 5) +
  annotate("segment",
           x=as.Date("2014-11-18"),
           xend=as.Date("2014-12-30"),
           y=-axisshift+0.03,
           yend=-axisshift+0.03,
           color ="black") +
  annotate("text",
           label = "2015",
           x=as.Date("2015-06-22"),
           y=-axisshift, 
           size = 5) +
  annotate("segment",
           x=as.Date("2015-01-03"),
           xend=as.Date("2015-12-30"),
           y=-axisshift+0.03,
           yend=-axisshift+0.03,
           color ="black") +
  annotate("text",
           label = "2016",
           x=as.Date("2016-06-22"),
           y=-axisshift, 
           size = 5) +
  annotate("segment",
           x=as.Date("2016-01-03"),
           xend=as.Date("2016-12-30"),
           y=-axisshift+0.03,
           yend=-axisshift+0.03,
           color ="black") +
  annotate("text",
           label = "2017",
           x=as.Date("2017-03-15"),
           y=-axisshift, 
           size = 5) +
  annotate("segment",
           x=as.Date("2017-01-03"),
           xend=as.Date("2017-05-15"),
           y=-axisshift+0.03,
           yend=-axisshift+0.03,
           color ="black") +
  coord_cartesian(ylim=c(-0.01,1), clip="off", expand = 0)+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        legend.position = "top",
        panel.grid = element_line(color = "grey70", linetype = 3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank(),
        axis.text = element_text(color = "black"),
        text = element_text(size = 17),
        panel.background = element_blank())
gg
ggsave(paste0("out/coiprop-methods.", imgfmt), 
       width = 12, height = 6, units = "in", dpi = 400)
#Estimates of proportion of monoclonal infections from:
#FWS on WGS 27k SNPs
round(length(which(coi$Fws>=0.95))/nrow(coi[!is.na(coi$Fws),]),2)
#Hetero on 89 loci
round(length(which(coi$prop.SNP < 0.005))/nrow(coi[!is.na(coi$prop.SNP),]),2)
#Hetero on 27k loci
round(length(which(coi$prop.WGS < 0.005))/nrow(coi[!is.na(coi$prop.WGS),]),2)

#Hetero loci SNP and WGS vs Fws
length(which(coi$prop.SNP < 0.005 & coi$Fws < 0.95))/length(which(!is.na(coi$Fws) & coi$prop.SNP < 0.005))
length(which(coi$prop.WGS < 0.005 & coi$Fws < 0.95))/length(which(!is.na(coi$Fws) & coi$prop.WGS < 0.005))

length(which(coi$prop.SNP >= 0.005 & coi$Fws >= 0.95))/length(which(!is.na(coi$Fws) & coi$prop.SNP >= 0.005))
length(which(coi$prop.WGS >= 0.005 & coi$Fws >= 0.95))/length(which(!is.na(coi$Fws) & coi$prop.WGS > 0.005))
lm.barcode <- lm(coi$prop.SNP ~ 
             coi$Fws)
summary(lm.barcode)
lm.wgs <- lm(coi$prop.WGS ~ 
                   coi$Fws)
summary(lm.wgs)

coiprop.melted[coiprop.melted$Numobs>=5,] %>%
  group_by(Method) %>%
  reframe(mincoi = min(COI),
          maxcoi = max(COI))
############################################
