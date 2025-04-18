###########################GLOBAL VARIABLES###########################
seasoncols <- c("springgreen3", "red3")
seasonnames <- c("Low", "High")
######################################################

dim(drugres.cons.improvedbc[!drugres.cons.improvedbc$Comparison%in%c("Both unknown", "Unknown vs. called"),])
nrow(drugres.cons.improvedbc[grepl("^Match", drugres.cons.improvedbc$Comparison),])
483/545
mycols <- brewer.pal(4, "Greens")[-1]
mycols <- c(mycols, brewer.pal(5, "Blues")[-1])
mycols <- c(mycols, brewer.pal(3, "Reds")[-1])

gg <- ggplot(data = drugres.cons.improvedbc[drugres.cons.improvedbc$ShortMarker!="*" & !drugres.cons.improvedbc$Comparison%in%c("Both unknown", "Unknown vs. called"),])+
  geom_bar(aes(x = 1,
               fill = factor(Comparison, levels = rev(c("Match (Both sensitive)", "Match (Both resistant)",
                                                    "Match (Both mixed)",
                                                    "Partial Match (Genome mixed, barcode sensitive)",
                                                    "Partial Match (Genome mixed, barcode resistant)",
                                                    "Partial Match (Barcode mixed, genome sensitive)",
                                                    "Partial Match (Barcode mixed, genome resistant)",
                                                    "Mismatch (Barcode sensitive, genome resistant)",
                                                    "Mismatch (Genome sensitive, barcode resistant)")))), 
           stat = "count", position = "stack")+
  scale_fill_manual("Agreement",
                    breaks = rev(c("Match (Both sensitive)", "Match (Both resistant)",
                               "Match (Both mixed)",
                               "Partial Match (Genome mixed, barcode sensitive)",
                               "Partial Match (Genome mixed, barcode resistant)",
                               "Partial Match (Barcode mixed, genome sensitive)",
                               "Partial Match (Barcode mixed, genome resistant)",
                               "Mismatch (Barcode sensitive, genome resistant)",
                               "Mismatch (Genome sensitive, barcode resistant)")),
                    values = rev(mycols),
                    guide = guide_legend(ncol = 1,override.aes = c(color = "black"),
                                         byrow = TRUE))+
  scale_x_continuous(breaks = c(1),
                     expand = c(0,0))+
  scale_y_continuous("Pairs of markers",
                     n.break = 10, expand = c(0,0),
                     limits = c(0,550))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        legend.spacing.y = unit(0.5, "cm"))
gg
ggsave(paste0("out/DR_agreement.", imgfmt),
       width = 8, height = 6, dpi = 500)


drugres.cons.filtered <- drugres.cons.summary[!drugres.cons.summary$ShortMarker.Cons%in%c("X"),]
drugres.cons.filtered <- merge(drugres.cons.filtered, drugphenotype[,c("Name", "AA")])

drugres.cons.filtered.stats <- drugres.cons.summary[!drugres.cons.summary$ShortMarker.Cons%in%c("X"),]%>%
  group_by(Type)%>%
  summarise(ngeno = length(ShortMarker),
            nsamp = length(unique(Sample)))

drugres.cons.size <- drugres.cons.filtered%>%
  group_by(Name, AA, date)%>%
  summarise(ntot = length(AAsen))

drugres.cons.rank <- drugres.cons.filtered[drugres.cons.filtered$ShortMarker.Cons!=drugres.cons.filtered$AAsen,]%>%
  group_by(Name, AA, AAsen)%>%
  summarise(rank = as.numeric(factor(unique(ShortMarker.Cons))),
            ShortMarker.Cons = unique(ShortMarker.Cons))

drugres.cons.grouped <- drugres.cons.filtered%>%
  group_by(Name, date,AAsen,ShortMarker.Cons)%>%
  summarise(nmut = length(ShortMarker.Cons))

drugres.cons.grouped <- drugres.cons.grouped[drugres.cons.grouped$ShortMarker.Cons!=drugres.cons.grouped$AAsen,]
drugres.cons.grouped <- merge(drugres.cons.grouped,
                              drugres.cons.size)
drugres.cons.grouped <- merge(drugres.cons.grouped,
                              drugres.cons.rank)

c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)
  
  name <- paste(a$name, b$name, sep = "-")
  
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format, 
            domain = c(as.POSIXct(-Inf), as.POSIXct(Inf)))
  
}
rev_date <- c_trans("reverse", "time")
ori_date <- c_trans("identity", "time")
#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(as.numeric(unique(format(drugres.cons.grouped$date,"%Y"))))
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
mindate <- min(drugres.cons.grouped$date)
maxdate <- max(drugres.cons.grouped$date)
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]
seasons.narrow <- seasons
seasons.narrow$date1[seasons.narrow$date1==min(seasons.narrow$date1)] <- "2014-10-01"
seasons.narrow$date2[seasons.narrow$date2==max(seasons.narrow$date2)] <- "2017-07-01"

my_palette <- brewer.pal(7, "Set1")
drugres.cons.grouped.mean <- drugres.cons.grouped[drugres.cons.grouped$ntot>=5,]%>%
  group_by(Name,ShortMarker.Cons)%>%
  summarise(allnmut = sum(nmut),
            allntot = sum(ntot))
drugres.cons.grouped.mean
drugres.cons.grouped.mean$allnmut/drugres.cons.grouped.mean$allntot

drugres.cons.grouped.mean.int1 <- drugres.cons.grouped[drugres.cons.grouped$ntot>=5 & drugres.cons.grouped$date <= as.Date("2016-08-01"),]%>%
  group_by(Name,ShortMarker.Cons)%>%
  summarise(allnmut = sum(nmut),
            allntot = sum(ntot))
drugres.cons.grouped.mean.int1
drugres.cons.grouped.mean.int1$allnmut/drugres.cons.grouped.mean.int1$allntot

drugres.cons.grouped.mean.int2 <- drugres.cons.grouped[drugres.cons.grouped$ntot>=5 & drugres.cons.grouped$date > as.Date("2016-08-01"),]%>%
  group_by(Name,ShortMarker.Cons)%>%
  summarise(allnmut = sum(nmut),
            allntot = sum(ntot))
drugres.cons.grouped.mean.int2
drugres.cons.grouped.mean.int2$allnmut/drugres.cons.grouped.mean.int2$allntot

drugres.cons.grouped$err <- 1.96*sqrt((drugres.cons.grouped$nmut/drugres.cons.grouped$ntot)*
  (1-(drugres.cons.grouped$nmut/drugres.cons.grouped$ntot))/
  drugres.cons.grouped$ntot)

drugres.cons.grouped$errwilson <- binom.confint(drugres.cons.grouped$nmut, 
                                                    drugres.cons.grouped$ntot, method=c("wilson"))

###MEAN PREVALENCE OF MARKERS###
drugres.cons.pooledstats <- drugres.cons.grouped%>%
  group_by(Name)%>%
  summarise(binom.confint(sum(nmut), sum(ntot), method=c("wilson")))
###

dim(drugres.cons.grouped)
length(which(drugres.cons.grouped$ntot<5))
summary(drugres.cons.grouped$ntot[drugres.cons.grouped$ntot>=5])
axisshift <- 1/14
gg <- ggplot(drugres.cons.grouped[drugres.cons.grouped$ntot>=5 & drugres.cons.grouped$Name!="K13",]) +
  geom_rect(data=seasons.narrow, aes(ymin=-Inf, ymax=Inf, 
                                     xmin=as.POSIXct(date1), xmax=as.POSIXct(date2), 
                                     fill=toupper(type)), alpha=0.125)+
  scale_fill_manual("Transmission season",
                    labels = seasonnames,
                    values=seasoncols,
                    guide = guide_legend(nrow = 2, order = 1, 
                                         override.aes = c(color = "black")))+
  new_scale("fill")+
  geom_segment(aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6), 
                   y = errwilson$lower+0.001,
                   xend = as.POSIXct(date+2*as.numeric(factor(Name))-6), 
                   yend = errwilson$upper-0.001,
                color = paste0(Name,",", rank)),
               alpha = 0.3, show.legend = F, linewidth = 0.8)+
  geom_segment(aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6)
                   -60*60*24*2.5, 
                   y = errwilson$upper,
                   xend = as.POSIXct(date+2*as.numeric(factor(Name))-6)
                   +60*60*24*2.5, 
                   yend = errwilson$upper,
                   color = paste0(Name,",", rank)),
               alpha = 0.3, show.legend = F)+
  geom_segment(aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6)-
                     60*60*24*2.5, 
                   y = errwilson$lower,
                   xend = as.POSIXct(date+2*as.numeric(factor(Name))-6)+
                     60*60*24*2.5, 
                   yend = errwilson$lower,
                   color = paste0(Name,",", rank)),
               alpha = 0.3, show.legend = F)+
  geom_point(aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6),
                 y = nmut/ntot,
                fill = paste0(Name,",", rank)),
             pch = 21, stroke = NA, size = 2.5)+
  geom_line(aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6), 
                y = nmut/ntot,
                color = paste0(Name,",", rank),
                linetype = paste0(Name,",", rank)), 
            show.legend = F)+
  geom_text(data = drugres.cons.grouped[drugres.cons.grouped$ntot>=5 & drugres.cons.grouped$date=="2014-12-01" & drugres.cons.grouped$Name!="K13",],
            aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6)-60*60*24*8,
                y = nmut/ntot,
                color = paste0(Name,",", rank),
                label = Name),
            fontface = "bold", show.legend = F, hjust = "right")+
  geom_text(data = drugres.cons.grouped[drugres.cons.grouped$ntot>=5 & drugres.cons.grouped$date=="2017-05-01" & drugres.cons.grouped$Name!="K13",],
            aes(x = as.POSIXct(date+2*as.numeric(factor(Name))-6)+60*60*24*8,
                y = nmut/ntot,
                color = paste0(Name,",", rank),
                label = Name),
            fontface = "bold", show.legend = F, hjust = "left")+
  scale_color_manual("Resistance genotype",
                    breaks = c("AAT1,1", "CRT,1", "DHFR,1", "DHFR,2",
                               "DHPS,1", "K13,2", "MDR1,1"),
                    labels = c("AAT1 S258L", "CRT K76T", "DHFR S108N", "DHFR S108T",
                               "DHPS A437G", "K13 C580Y", "MDR1 N86Y"),
                    values = c(my_palette[1],my_palette[2],my_palette[3],my_palette[3],
                               my_palette[4],my_palette[5],my_palette[7]),
                    guide = guide_legend(nrow = 2))+
  scale_fill_manual("Resistance genotype",
                     breaks = c("AAT1,1", "CRT,1", "DHFR,1", "DHFR,2",
                                "DHPS,1", "K13,2", "MDR1,1"),
                     labels = c("AAT1 S258L", "CRT K76T", "DHFR S108N", "DHFR S108T",
                                "DHPS A437G", "K13 C580Y", "MDR1 N86Y"),
                     values = c(my_palette[1],my_palette[2],my_palette[3],my_palette[3],
                                my_palette[4],my_palette[5],my_palette[7]),
                     guide = guide_legend(nrow = 2, override.aes = c(size = 5, stroke = 1)))+
  scale_linetype_manual("Resistance genotype",
                        breaks = c("AAT1,1", "CRT,1", "DHFR,1", "DHFR,2",
                                   "DHPS,1", "K13,2", "MDR1,1"),
                        labels = c("AAT1 S258L", "CRT K76T", "DHFR S108N", "DHFR S108T",
                                   "DHPS A437G", "K13 C580Y", "MDR1 N86Y"),
                     values = c(1,1,1,2,1,1,1),
                     guide = guide_legend(nrow = 2))+
  scale_x_continuous(trans = ori_date,
                     expand = c(0,0), 
               breaks = as.POSIXct("2014-09-01")+months(seq(2,33,2)), 
               labels = c("Nov", 
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May", "Jul", "Sep", "Nov",
                          "Jan", "Mar", "May"),
               limits = c(as.POSIXct("2014-10-01"),as.POSIXct("2017-07-02")))+
  scale_y_continuous("Resistance marker prevalence",
                     n.breaks = 10,
                     expand = c(0,0))+
  ylab("")+
  xlab("")+
  coord_cartesian(ylim=c(-0.01,1.01), 
                  clip="off")+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        text=element_text(size=15),
        legend.position = "top",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank())
gg <- gg+
  annotate("text",
           label = "2014",
           size = 4,
           x=as.POSIXct("2014-12-7"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2014-11-18"),
           xend=as.POSIXct("2014-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2015",
           size = 4,
           x=as.POSIXct("2015-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2015-01-03"),
           xend=as.POSIXct("2015-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2016",
           size = 4,
           x=as.POSIXct("2016-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2016-01-03"),
           xend=as.POSIXct("2016-12-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60") +
  annotate("text",
           label = "2017",
           size = 4,
           x=as.POSIXct("2017-03-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2017-01-03"),
           xend=as.POSIXct("2017-05-30"),
           y=-axisshift+0.02,
           yend=-axisshift+0.02,
           color ="grey60")
gg
  
ggsave(paste0("out/DR-consensus-prevalence.", imgfmt), 
       width = 13, height = 8, dpi = 400)

drugres.cons.grouped.DHPS.1 <- drugres.cons.grouped %>%
  filter(ntot>=5 & date >= as.Date("2014-12-01") & date <= as.Date("2016-07-01") & Name == "DHPS")
lm.DHPS.1 <- lm(100*drugres.cons.grouped.DHPS.1$nmut/drugres.cons.grouped.DHPS.1$ntot  ~ as.numeric(drugres.cons.grouped.DHPS.1$date))

drugres.cons.grouped.DHPS.2 <- drugres.cons.grouped %>%
  filter(ntot>=5 & date >= as.Date("2016-10-01") & Name == "DHPS")
lm.DHPS.2 <- lm(100*drugres.cons.grouped.DHPS.2$nmut/drugres.cons.grouped.DHPS.2$ntot  ~ as.numeric(drugres.cons.grouped.DHPS.2$date))
summary(lm.DHPS.2)
lm.DHPS.df <- data.frame("starting" = c(as.numeric(as.Date("2014-12-01")), as.Date("2016-10-01")), 
                         "ending" = c(as.numeric(as.Date("2016-07-01")), as.Date("2017-05-01")),
                         "intercept" = c(summary(lm.DHPS.1)$coefficients[1,1], summary(lm.DHPS.2)$coefficients[1,1]),
                         "regression" = c(summary(lm.DHPS.1)$coefficients[2,1], summary(lm.DHPS.2)$coefficients[2,1]))

gg2 <- gg+
  geom_segment(data = lm.DHPS.df,
               mapping = aes(x = as.POSIXct(as.Date(starting)), xend = as.POSIXct(as.Date(ending)),
                   y = (intercept + regression*starting)/100, 
                   yend = (intercept + regression*ending)/100),
               linetype = 5, linewidth = 1,alpha = 0.7, color = my_palette[4])

gg2
##############################

