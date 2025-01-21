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
green.midlight <- colalpha2hex("springgreen3", 0.5)
red.midlight <- colalpha2hex("red3", 0.5)
seasoncols.light <- c(green.light, red.light)
seasoncols.midlight <- c(green.midlight, red.midlight)
seasoncols <- c("springgreen3", "red3")
seasonnames <- c("Low", "High")
######################################################

###########################EPIDEMIO###########################

#Epidemiologic data
posneg <- read.csv("rawdata/pf-test.csv")
colnames(posneg) <- c("sampleID", "ParticipantID", "Cohort",
                      "date", "collection", "infectivity", "treatment")
posneg$date <- as.Date(paste0(posneg$date, "-01"))
posneg$infectivity[posneg$sampleID%in%nodes$ID] <- "POS"
######################################################

###########################DURATION OF INFECTION###########################

ibdcutoff <- 0.9
doiclusters.edges <- edges[edges$fract_sites_IBD >= ibdcutoff & edges$sameIndividual == T,]

# Get rid of loops and ensure right naming of vertices
topology <- doiclusters.edges[,c("ID1", "ID2")]
topology <- simplify(graph_from_data_frame(topology[order(topology[[1]]),],directed = FALSE))

# Find all components
comps <- components(topology)
cluster.df <- data.frame(comps$membership)
cluster.df <- cbind(row.names(cluster.df), cluster.df)
colnames(cluster.df) <- c("ID", "clusterind")
doiclusters.nodes <- merge(nodes, cluster.df, all.x = T)
doiclusters.nodes$clusterind[is.na(doiclusters.nodes$clusterind)] <- 0
doiclusters.nodes$clusterind <- sprintf(paste0("%0",max(nchar(doiclusters.nodes$clusterind)),"d"), doiclusters.nodes$clusterind)

lineages <- doiclusters.nodes[doiclusters.nodes$clusterind != "00",] %>%
  group_by(clusterind, individual)%>%
  reframe(start = min(date),
          end = max(date))

studyParticipants <- read.csv("rawdata/study-participants.csv", header = T)
studyParticipants$DOB <- 2014-studyParticipants$Age..2014
lineages <- merge(lineages,
                       studyParticipants,
                       by.x = "individual", by.y = "Subject_ID", all.x = T)
lineages <- lineages[!is.na(lineages$DOB),]
lineages$startage <- time_length(difftime(as.Date(lineages$start), 
                                               as.Date(paste0(lineages$DOB, "-01-01"))), "years")
lineages <- lineages[lineages$startage > 5,]

lineages$individual <- factor(lineages$individual,
                              levels = unique(lineages$individual[with(lineages, order(start - end))]))
lineages$start <- as.POSIXct(lineages$start)
lineages$end <- as.POSIXct(lineages$end)

lineages$relpos <- 0
lineages$relpos[which(duplicated(lineages$individual, fromLast = T))] <- 0.15
lineages$relpos[which(duplicated(lineages$individual))] <- -0.15

posneg.subset <- merge(posneg, lineages, by.x = "ParticipantID", by.y = "individual")
posneg.subset <- posneg.subset[!is.na(posneg.subset$date) & (posneg.subset$date + 90 > posneg.subset$start) &
                                 (posneg.subset$date - 90 < posneg.subset$end),]
posneg.subset <- posneg.subset[!duplicated(posneg.subset[,c("sampleID", "date")]),]
posneg.subset <- posneg.subset[!is.na(posneg.subset$infectivity),]

nodes.subset <- merge(doiclusters.nodes, lineages[,c("individual", "start", "end", "clusterind")],
                      by = c("individual"), suffixes = c("", ".main"))
nodes.subset <- nodes.subset[(nodes.subset$date + 90 > nodes.subset$start) &
                               (nodes.subset$date - 90 < nodes.subset$end),]
nodes.subset <- nodes.subset[!paste0(nodes.subset$individual, nodes.subset$date)%in%paste0(lineages$individual, lineages$start),]
nodes.subset <- nodes.subset[!paste0(nodes.subset$individual, nodes.subset$date)%in%paste0(lineages$individual, lineages$end),]
nodes.subset$type <- ifelse(nodes.subset$clusterind == nodes.subset$clusterind.main,
                            "main", "add")

seasons.narrow <- seasons
seasons.narrow$date1[seasons.narrow$date1==min(seasons.narrow$date1)] <- "2014-11-15"
seasons.narrow$date2[seasons.narrow$date2==max(seasons.narrow$date2)] <- "2017-05-30"

table(lineages[!duplicated(lineages$individual),"Gender"])
as.numeric(lineages[duplicated(lineages$individual),"individual"])

sort(lineages$individual)[19]
axisshift <- length(unique(lineages$individual))/20

gg <- ggplot(lineages)+
  geom_rect(data=seasons.narrow, aes(ymin=-Inf, ymax=Inf, 
                              xmin=as.POSIXct(date1), xmax=as.POSIXct(date2), 
                              fill=toupper(type)))+
  scale_fill_manual("Transmission season",
                    values=seasoncols.light,
                    labels = seasonnames,
                    guide = guide_legend(nrow = 2, order = 1,
                                         override.aes = list(color = "black")))+
  new_scale("fill")+
  geom_rect(data = posneg.subset,
             aes(ymin = as.numeric(factor(ParticipantID, 
                                        levels = levels(lineages$individual)))-0.325,
                 ymax = as.numeric(factor(ParticipantID, 
                                          levels = levels(lineages$individual)))+0.325,
                 xmin = as.POSIXct(date-4.5), xmax = as.POSIXct(date+4.5),
                 fill = infectivity))+
  scale_fill_manual("P.f. infection",
                     breaks = c("NEG", "POS"),
                     labels = c("Negative", "Positive"),
                     values = seasoncols.midlight,
                     guide = guide_legend(nrow = 2, order = 2, override.aes = list(color = "black")))+
  geom_segment(aes(y=as.numeric(individual)-0.5, 
                   yend=as.numeric(individual)-0.5,
                   x=as.POSIXct("2014-11-15"), 
                   xend=as.POSIXct("2017-05-31")), 
               linewidth = 0.25, color = "grey40", linetype = 3)+
  geom_segment(data=data.frame(1),
               y=max(as.numeric(lineages$individual))+0.5, 
               yend=max(as.numeric(lineages$individual))+0.5,
               x=as.POSIXct("2014-11-15"), 
               xend=as.POSIXct("2017-05-31"), 
               linewidth = 0.25, color = "grey40", linetype = 3)+
  new_scale("fill")+
  geom_point(data = nodes.subset,
             aes(x = as.POSIXct(date), y = as.numeric(factor(individual, levels = levels(lineages$individual))),
                 fill = type),
             color = "black", pch = 21, size = 2, show.legend = c("fill" = T, "color" = F))+
  geom_segment(data = lineages[lineages$relpos == 0,],
               aes(y=as.numeric(individual), 
                 yend=as.numeric(individual),
                 x=start, xend=end, group = individual), linewidth = 0.5,
             show.legend = F)+
  geom_curve(data = lineages[lineages$relpos < 0,],
             aes(y=as.numeric(individual)+relpos, 
                 yend=as.numeric(individual)+relpos,
                 x=start, xend=end, group = individual), linewidth = 0.5,
             show.legend = F, arrow = arrow(length = unit(0, "inches"), type = "closed"),
             curvature = 0.12)+
  geom_curve(data = lineages[lineages$relpos > 0,],
             aes(y=as.numeric(individual)+relpos, 
                 yend=as.numeric(individual)+relpos,
                 x=start, xend=end, group = individual), linewidth = 0.5,
             show.legend = F, arrow = arrow(length = unit(0, "inches"), type = "closed"),
             curvature = -0.12)+
  geom_point(aes(y = as.numeric(individual)+relpos, 
                 x = start, 
                 group = individual, fill = "main"), size = 2,
             pch = 21, show.legend = c("fill" = T, "color" = F))+
  geom_point(aes(y = as.numeric(individual)+relpos, 
                 x = end, 
                 group = individual, fill = "main"), size = 2,
             pch = 21, show.legend = c("fill" = T, "color" = F))+
  scale_fill_manual("Barcode", breaks = c("main", "add"),
                     values = c("main" = "black", "add" = "grey90"),
                     labels = c("main" = "Dominant genotype", "add" = "Additional genotype"),
                     guide = guide_legend(order = 3, nrow = 2, override.aes = c(size = 4)))+
  scale_x_continuous(trans = ori_date, breaks = as.POSIXct(unique(c(nodes$date))), 
                     labels = month(unique(c(as.Date(nodes$date))), 
                                    label = T, abbr = T),
                     expand = c(0,0))+
  scale_y_continuous(breaks = c(1:length(unique(lineages$individual))),
                     labels = c(1:length(unique(lineages$individual))),
                     expand = c(0.017,0))+
  ylab("Participant rank")+
  xlab("")+
  annotate("text",
           label = "2014",
           size = 6,
           x=as.POSIXct("2014-12-7"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2014-11-18"),
           xend=as.POSIXct("2014-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           linewidth = 0.75,
           color ="black") +
  annotate("text",
           label = "2015",
           size = 6,
           x=as.POSIXct("2015-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2015-01-03"),
           xend=as.POSIXct("2015-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           linewidth = 0.75,
           color ="black") +
  annotate("text",
           label = "2016",
           size = 6,
           x=as.POSIXct("2016-06-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2016-01-03"),
           xend=as.POSIXct("2016-12-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           linewidth = 0.75,
           color ="black") +
  annotate("text",
           label = "2017",
           size = 6,
           x=as.POSIXct("2017-03-22"),
           y=-axisshift) +
  annotate("segment",
           x=as.POSIXct("2017-01-03"),
           xend=as.POSIXct("2017-05-30"),
           y=-axisshift+0.75,
           yend=-axisshift+0.75,
           linewidth = 0.75,
           color ="black") +
  coord_cartesian(ylim=c(1,length(unique(lineages$individual))), 
                  clip="off")+
  theme(plot.margin = unit(c(1,1,2,1), "lines"),
        text=element_text(size=22),
        legend.position = "top",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text = element_text(color = "black", size = 15),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.key = element_blank())
  
gg
ggsave(paste0("out/doi-", sub("\\.", "", as.character(ibdcutoff), perl = F), ".png"),
       width = 15, height = 9, dpi = 600)

#Ranks of individuals by duration of infection
lineages$days <- as.numeric(lineages$end - lineages$start)
doidata <- lineages
doidata <- doidata[,c("individual", "days", "start", "end")]
doidata$start <- substr(doidata$start, 1, 7)
doidata$end <- substr(doidata$end, 1, 7)
doidata$rank <- as.numeric(doidata$individual)
colnames(doidata) <- c("Individual", "Estimated duration of continuous infection (days)",
                       "Estimated start of infection (YYYY-MM)", "Estimated end of infection (YYYY-MM)", "Rank in figure 6")
write.csv(doidata, "out/doi-metadata.csv", quote = F, row.names = F)

femaledoi <- lineages$days[lineages$Gender=="F"]
maledoi <- lineages$days[lineages$Gender=="M"]


shortvlongdoi.gender <- data.frame(list("F"=c(length(femaledoi[femaledoi<90]), length(femaledoi[femaledoi>=90])), 
             "M" = c(length(maledoi[maledoi<90]), length(maledoi[maledoi>=90]))))
row.names(shortvlongdoi.gender) <- c("short", "long")
shortvlongdoi.gender[1,]/(shortvlongdoi.gender[1,]+shortvlongdoi.gender[2,])
chisq.test(shortvlongdoi.gender, correct = F)
fisher.test(shortvlongdoi.gender)

young <- lineages$days[lineages$startage < 15]
old <- lineages$days[lineages$startage >= 15]

shortvlongdoi.age <- data.frame(list("young" = c(length(young[young<90]), length(young[young>=90])),
                                     "old" = c(length(old[old<90]), length(old[old>=90]))))
row.names(shortvlongdoi.age) <- c("short", "long")

1-shortvlongdoi.age[1,]/(shortvlongdoi.age[1,]+shortvlongdoi.age[2,])

chisq.test(shortvlongdoi.age, correct = F)
t.test(young, old, alternative = "two.sided", var.equal = F)

mycols <- brewer.pal(3, "Set1")
gg <- ggplot(lineages)+
  geom_vline(xintercept = c(5.1,15,24.9),
             linetype = 1)+
  geom_violin(aes(x = pmin(round(startage/10)*10,20)+5*(as.numeric(factor(Gender))-1.5),
                  group = pmin(round(startage/10)*10,20)+5*(as.numeric(factor(Gender))-1.5),
                  y = round(days/30)),
              scale = "width")+
  geom_hline(yintercept = 3, linetype = 2, color = "grey40", linewidth = 1)+
  geom_point(aes(x = pmin(round(startage/10)*10,20)+5*(as.numeric(factor(Gender))-1.5),
                 group = factor(pmin(round(startage/10)*10,20)+5*(as.numeric(factor(Gender))-1.5)),
                 y = round(days/30),
                 fill = Gender),
             alpha = 1, pch = 21, size = 2, 
             position = position_jitter(width = 1.1, seed = 789, height = 0))+
  scale_fill_manual("",
                    breaks = c("F", "M"),
                    values = c(mycols[1], mycols[2]),
                    labels = c("Female", "Male"),
                    guide = guide_legend(override.aes = c(size = 4)))+
  xlab("Age group")+
  ylab("Duration of infection (months)")+
  scale_x_continuous(breaks = c(10, 20),
                     labels = c("5-15 years old", "Older than 15 years"),
                     expand = c(0,0), limits = c(4.99,25.01))+
  scale_y_continuous(n.breaks = 15, expand = c(0,0), limits = c(0,17.5))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        text = element_text(size = 18),
        panel.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_line(linetype = 3, color = "grey70"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave("out/doi-age.png", width = 8, height = 6, dpi = 350)

######################################################
