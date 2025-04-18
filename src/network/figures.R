###########################GLOBAL VARIABLES###########################
seasoncols <- c("springgreen3", "red3")
seasonnames <- c("Low", "High")
######################################################

###########STATS##############
if(removecontinf){
  outfolder="out/nocontinf/"
}else
{
  outfolder="out/"
}
############################

##############FIGURES##############

##############DATES################
gg <- ggplot(edges.grouped.datecompound)+
  geom_boxplot(aes(x = group,
                   group = group,
               y = horizontalIBD/horizontalMax),
               color = "black")+
  geom_segment(data=data.frame(1),
               aes(y = 0.2325, yend = 0.2325,
                   x = 1,
                   xend = 5), 
               color = "red3", linetype = 2, show.legend = F, linewidth = 0.5)+
  geom_segment(data=data.frame(c(1,2)),
               aes(y = c(0.2285, 0.2285), yend = c(0.2325, 0.2325),
                   x = c(1,5),
                   xend = c(1,5)), 
               color = "red3", linetype = 1, show.legend = F, linewidth = 0.5)+
  geom_text(data=data.frame(1),
            aes(y = 0.2375,
                x = 3,
                label = strrep("\U2731",1+floor(log(0.05/less2_more12.test$p.value, 
                                               base = 10)))),
            fontface = "bold", size = 4, color = "red3")+
  scale_x_continuous("Months between sampled pairs",
                     breaks = dfgrouplocation.datecompound$group,
                     labels=paste0(dfgrouplocation.datecompound$mingrouplocation, 
                                   "-",
                                   dfgrouplocation.datecompound$maxgrouplocation))+
  scale_y_continuous("Percentage of related isolates (IBD > 0.5)",
                     label = function(x){x*100},
                     breaks=seq(0,0.23, 0.01), limits = c(0,0.2385),expand = c(0.01,0))+
  theme(text = element_text(size=17), panel.grid.minor = element_blank(),
        legend.position="inside", legend.position.inside = c(0.8,0.8),
        legend.key.size = unit(2, "lines"), legend.spacing.x = unit(1, "lines"),
        legend.background = element_rect(fill = NA), legend.text = element_text(size=11),
        legend.box.background = element_rect(fill = "white", color ="black"),
        legend.key=element_blank(), legend.title = element_blank(),
        legend.spacing.y = unit(0, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(r = 0.1, l = 0.1, t = 0.1, b = 0.1, unit = "lines"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))
gg
ggsave(paste0(outfolder, "IBD-decrease-time.", imgfmt), gg, 
       width = 12, height = 6, dpi = 500)

gg <- ggplot(edges.grouped.datecompound)+
  geom_boxplot(aes(x=group, 
                   y=horizontalIBD/horizontalMax, 
                   fill=paste0(sameVillage, sameCompound),
                   group=factor(paste0(group, sameVillage, sameCompound),
                                levels = as.vector(outer(unique(sort(group)), 
                                                         c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"), 
                                                         paste0)))), 
               outlier.shape = NA)+
  geom_point(aes(x=group, 
                 y=horizontalIBD/horizontalMax, 
                 fill=paste0(sameVillage, sameCompound),
                 group=factor(paste0(group, sameVillage, sameCompound),
                              levels =as.vector(outer(unique(sort(group)), 
                                                      c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"), 
                                                      paste0)))), 
             pch=21, show.legend = F, position=position_jitterdodge(jitter.width = 0.1, seed = 123))+
  geom_segment(data=geotests,
               aes(y = c(0.235,0.2425,0.195), yend = c(0.235,0.2425,0.195),
                   x = fromgroup,
                   xend = togroup), 
               color = "red3", linetype = 2, show.legend = F, linewidth = 0.5)+
  geom_segment(data=geotests,
               aes(y = c(0.23,0.2375,0.19), yend = c(0.235,0.2425,0.195),
                   x = fromgroup,
                   xend = fromgroup), 
               color = "red3", linetype = 1, show.legend = F, linewidth = 0.5)+
  geom_segment(data=geotests,
               aes(y = c(0.23,0.2375,0.19), yend = c(0.23,0.2375,0.19)+0.005,
                   x = togroup,
                   xend = togroup), 
               color = "red3", linetype = 1, show.legend = F, linewidth = 0.5)+
  geom_text(data=geotests,
            aes(y = c(0.235,0.2425,0.195)+0.004,
                x = c((fromgroup+togroup)/2),
                label = ifelse(pvalues > 0.05, "ns", strrep("\U2731",1+floor(log(0.05/pvalues, 
                                                    base = 10))))),
            fontface = "bold", size = 4, color = "red3")+
  # stat_summary(aes(x = group,
  #                  y = horizontalIBD/horizontalMax,
  #                  group = paste0(sameVillage, sameCompound),
  #                  color = paste0(sameVillage, sameCompound)),
  #              fun = median, geom = "line", position = position_dodge(width = 0.75),
  #              show.legend = F, linetype = 2)+
  scale_x_continuous("Months between sampled pairs",
                     breaks = dfgrouplocation.datecompound$group,
                     labels=paste0(dfgrouplocation.datecompound$mingrouplocation, 
                                   "-",
                                   dfgrouplocation.datecompound$maxgrouplocation))+
  scale_y_continuous("Percentage of related isolates (IBD > 0.5)",
                     label = function(x){x*100},
                     breaks=seq(0,0.23, 0.01), limits = c(0,0.2485), expand = c(0.01,0))+
  scale_fill_brewer(name = "",
                    breaks = c("TRUETRUE", "TRUEFALSE", "FALSEFALSE"),
                    labels = list("TRUETRUE"= "Same household", 
                                  "TRUEFALSE" = "Different households (same village)", 
                                  "FALSEFALSE" = "Different households (different villages)"),
                    palette = 9, type = "seq")+
  theme(text = element_text(size=17), panel.grid.minor = element_blank(),
        legend.position = "inside", legend.position.inside = c(.8,.8), legend.key.size = unit(2, "lines"), 
        legend.spacing.x = unit(1, "lines"),
        legend.background = element_rect(fill = NA), legend.text = element_text(size=12),
        legend.box.background = element_rect(fill = "white", color ="black"),
        legend.key=element_blank(), legend.title = element_blank(),
        legend.spacing.y = unit(0, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(r = 0.1, l = 0.1, t = 0.1, b = 0.1, unit = "lines"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid = element_line(color = "grey70", linetype = 3),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))

gg
ggsave(paste0(outfolder, "compound-time.", imgfmt), gg, width = 12, height = 6, dpi = 500)

my_gradient <- brewer.pal(n = 9, "PuBu")
starcol <- my_gradient[7]

gg <- ggplot(edges.grouped.date.1season)+
  geom_blank(aes(color=paste0(commonSeason, cycleElapsed)))+
  geom_boxplot(aes(x=factor(paste0(commonSeason, cycleElapsed),
                           levels = c("wet0", "drywet0", "dry0", "drywet1")),
                  y=horizontalIBD/horizontalMax),
               outlier.shape = NA)+
  geom_point(aes(x=paste0(commonSeason, cycleElapsed),
                 y=horizontalIBD/horizontalMax, 
                 fill=paste0(commonSeason, cycleElapsed)), 
             pch=21, size=3, position=position_jitter(height=0, seed=100, 
                                                      width = 0.15), show.legend = F)+
  geom_segment(data=drywettests,
               aes(x = from, xend = to,
                   y = c(0.065,0.075,0.07),
                   yend = c(0.065,0.075,0.07)), 
               linetype = 5, show.legend = F, color = starcol)+
  geom_segment(data=drywettests,
               aes(x = to, xend = to,
                   y = c(0.065,0.075,0.07)-0.0001,
                   yend = c(0.0625,0.0725,0.0675)+0.0001), 
               linetype = 1, show.legend = F, color = starcol)+
  geom_segment(data=drywettests,
               aes(x = from, xend = from,
                   y = c(0.065,0.075,0.07)-0.0001,
                   yend = c(0.0625,0.0725,0.0675)+0.0001), 
               linetype = 1, show.legend = F, color = starcol)+
  geom_text(data=drywettests,
            aes(x = c(3.5,2.5,3),
                y = c(0.0665,0.0765,0.0715),
                label = strrep("\U2731",1+floor(log(0.05/pvalues, 
                                               base = 10)))),
            fontface = "bold", size = 4,color = starcol)+
  scale_color_manual(breaks = c("dry0", "wet0", "drywet0", "drywet1"),
                     values = c(seasoncols, "orange", "gold"))+
  scale_fill_manual(breaks = c("dry0", "wet0", "drywet0", "drywet1"),
                    values = c(seasoncols, "orange", "gold"))+
  scale_x_discrete("Transmission seasons",
                   breaks = c("wet0", "drywet0", "dry0", "drywet1"),
                   labels = c("Within high",
                              "High to low", 
                              "Within low",
                              "Low to high"), expand = c(0.125,0.125))+
  scale_y_continuous("Percentage of related isolates (IBD > 0.5)",
                     labels = function(x){x*100},
                     breaks = c(seq(0,1,0.01)),
                     limits = c(-0.00125,0.0795), expand = c(0,0))+
  theme(text=element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        axis.text = element_text(color = "black", size = 15),
        panel.grid.major.y = element_line(color = "grey70", linetype = 3))
gg
ggsave(paste0(outfolder, "wetdry-comparison.", imgfmt), plot = gg,
       width = 9, height = 5, dpi = 500)

seasons.combination.nocoh$groups <- factor(paste0(seasons.combination.nocoh$season1, 
                                                  seasons.combination.nocoh$season2),
                                     levels = c("drydry", "wetdry",
                                                "drywet", "wetwet"))

seasons.1season <- seasons.combination.nocoh[(seasons.combination.nocoh$groups == "wetwet" | 
                                                seasons.combination.nocoh$groups == "drydry" | 
                                                seasons.combination.nocoh$groups == "wetdry") &
                                               seasons.combination.nocoh$cycle1 == seasons.combination.nocoh$cycle2,]
seasons.1season <- rbind(seasons.1season,
                         seasons.combination.nocoh[seasons.combination.nocoh$groups == "drywet" &
                                                     seasons.combination.nocoh$cycle2 -
                                                     seasons.combination.nocoh$cycle1 == 1,])
#my_gradient <- brewer.pal(n = 9, "YlOrRd")

gg <- ggplot(data = edges.grouped.date)+
  geom_segment(data = seasons.combination.nocoh[seasons.combination.nocoh$date1end >= seasons.combination.nocoh$date2start,],
               aes(y = as.POSIXct(date1start),
                   yend = as.POSIXct(date1end),
                   color = season1,
                   x = as.POSIXct(date2start),
                   xend = as.POSIXct(date2start)),
               show.legend = NA, linetype = 1, linewidth = 1.25, alpha = 0.5)+
  geom_segment(data = seasons.combination.nocoh[seasons.combination.nocoh$date2start >= seasons.combination.nocoh$date1start,],
               aes(x = as.POSIXct(date1start),
                   xend = as.POSIXct(date1end),
                   color = season1,
                   y = as.POSIXct(date2start),
                   yend = as.POSIXct(date2start)),
               show.legend = F, linetype = 1, linewidth = 1.25, alpha = 0.5)+
  geom_segment(data = seasons.nocoh,
               aes(x = as.POSIXct(date1),
                   xend = as.POSIXct(date2),
                   color = type,
                   y = as.POSIXct(max(date2)),
                   yend = as.POSIXct(max(date2))), linetype = 1,
               show.legend = F, linewidth = 1.25, alpha = 0.5)+
  geom_segment(data = seasons.nocoh[which.max(seasons.nocoh$date1),],
               aes(x = as.POSIXct(max(date2)),
                   xend = as.POSIXct(max(date2)),
                   color = type,
                   y = as.POSIXct(date1),
                   yend = as.POSIXct(date2)), linetype = 1,
               show.legend = F, linewidth = 1.25, alpha = 0.5)+
  geom_point(data = seasons.combination.nocoh,
             aes(y = as.POSIXct(date1start),
                 x = as.POSIXct(date2start)),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  geom_point(data = seasons.nocoh,
             aes(x = as.POSIXct(date1),
                 y = as.POSIXct(max(date2))),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  geom_point(data = seasons.nocoh[which.max(seasons.nocoh$date1),],
             aes(x = as.POSIXct(max(date2)),
                 y = as.POSIXct(date1)),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  geom_point(data = seasons.nocoh[which.max(seasons.nocoh$date1),],
             aes(x = as.POSIXct(max(date2)),
                 y = as.POSIXct(date2)),
             show.legend = NA, size = 1.25, color = "white", pch = 15)+
  scale_color_manual("Transmission season", values = seasoncols, labels = seasonnames,
                     guide = guide_legend(order = 1, ncol = 2, override.aes = c(linewidth = 2)))+
  new_scale("color")+
  geom_tile(data = seasons.1season,
            aes(x = as.POSIXct((as.numeric(date1start)+as.numeric(date1end))/2),
                y = date2start+60*60*24*22.5,
                color = groups),
            width = 60*60*24*95, height = 60*60*24*20, 
            fill = NA, linewidth = 1, linetype = 1, show.legend = F)+
  geom_text(data = seasons.1season,
            aes(x = as.POSIXct((as.numeric(date1start)+as.numeric(date1end))/2),
                y = date2start+60*60*24*22.5,
                label = factor(groups, levels = c("wetwet", "drydry",
                                                  "drywet", "wetdry"),
                               labels = c("Within high", "Within low", 
                                          "Low to high", "High to low"))),
            show.legend = F, fontface = "bold")+
  scale_color_manual("Seasonal groups", breaks = c("drydry", "wetwet",
                                                   "wetdry", "drywet"), 
                     values = c(seasoncols, "orange", "gold"))+
  geom_point(aes(x = as.POSIXct(datemin),
                 y = as.POSIXct(datemax),
                 size = horizontalMax),
             show.legend = NA, color = "black")+
  scale_size_continuous("Number of comparisons", breaks = c(5,10,100,500,1000,1500),
                        range = c(2,10),
                        guide = guide_legend(nrow = 2, order = 2))+
  new_scale("color")+
  new_scale("size")+
  geom_point(aes(x = as.POSIXct(datemin),
                 y = as.POSIXct(datemax),
                 size = horizontalMax,
                 color = horizontalIBD/horizontalMax),
             show.legend = NA)+
  scale_color_stepsn("Percentage of related isolates", colours = my_gradient, 
                       breaks = seq(0,1,0.01),
                     label = function(x){x*100},
                       guide = guide_colorbar(order = 3, title.position = "top",
                                              direction = "horizontal", barwidth = 17,
                                              label.hjust = -0.001))+
  scale_size_continuous("Number of comparisons", breaks = c(5,10,100,500,1000,1500),
                        range = c(1.5,9.5),
                        guide = guide_legend(nrow = 2, order = 2))+
  scale_y_continuous(trans = rev_date, 
                     breaks = unique(as.POSIXct(c(edges.grouped.date.not1season$datemin, edges.grouped.date.not1season$datemax))), 
                     labels = month(unique(as.Date(c(edges.grouped.date.not1season$datemin, 
                                                     edges.grouped.date.not1season$datemax))), label = T),
                     expand = c(0.01,0))+
  scale_x_continuous(trans = ori_date, 
                     breaks = unique(as.POSIXct(c(edges.grouped.date.not1season$datemin, edges.grouped.date.not1season$datemax))), 
                     labels = month(unique(as.Date(c(edges.grouped.date.not1season$datemin, 
                                                     edges.grouped.date.not1season$datemax))), label = T),
                     expand = c(0.01,0))+
  annotate("text",
           label = "2014",
           x=as.POSIXct("2014-10-23"),
           y=as.POSIXct("2017-01-15"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2014-08-18"),
           xend=as.POSIXct("2014-12-30"),
           y=as.POSIXct("2017-01-08"),
           yend=as.POSIXct("2017-01-08"),
           color ="black") +
  annotate("text",
           label = "2015",
           x=as.POSIXct("2015-06-15"),
           y=as.POSIXct("2017-01-15"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2015-01-03"),
           xend=as.POSIXct("2015-12-30"),
           y=as.POSIXct("2017-01-08"),
           yend=as.POSIXct("2017-01-08"),
           color ="black") +
  annotate("text",
           label = "2016",
           x=as.POSIXct("2016-06-08"),
           y=as.POSIXct("2017-01-15"), 
           size = 3.5) +
  annotate("segment",
           x=as.POSIXct("2016-01-03"),
           xend=as.POSIXct("2016-12-15"),
           y=as.POSIXct("2017-01-08"),
           yend=as.POSIXct("2017-01-08"),
           color ="black") +
  annotate("text",
           label = "2014",
           y=as.POSIXct("2014-10-23"),
           x=as.POSIXct("2014-07-13"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2014-08-18"),
           yend=as.POSIXct("2014-12-30"),
           x=as.POSIXct("2014-07-20"),
           xend=as.POSIXct("2014-07-20"),
           color ="black") +
  annotate("text",
           label = "2015",
           y=as.POSIXct("2015-06-15"),
           x=as.POSIXct("2014-07-13"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2015-01-03"),
           yend=as.POSIXct("2015-12-30"),
           x=as.POSIXct("2014-07-20"),
           xend=as.POSIXct("2014-07-20"),
           color ="black") +
  annotate("text",
           label = "2016",
           y=as.POSIXct("2016-06-08"),
           x=as.POSIXct("2014-07-13"), 
           size = 3.5, angle = -90, hjust = 0.5) +
  annotate("segment",
           y=as.POSIXct("2016-01-03"),
           yend=as.POSIXct("2016-12-15"),
           x=as.POSIXct("2014-07-20"),
           xend=as.POSIXct("2014-07-20"),
           color ="black") +
  theme_classic()+
  coord_cartesian(ylim = c(as.POSIXct("2016-12-9"),
                           as.POSIXct("2014-08-20")),
                  xlim = c(as.POSIXct("2014-08-20"),
                           as.POSIXct("2016-12-9")),
                  clip="off")+
  theme(plot.margin = unit(c(0,1,1,1), "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        axis.text.y = element_text(angle = -90, hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.ticks = element_blank(),
        text = element_text(size = 14),
        axis.line = element_blank(),
        legend.background = element_rect(color = "grey20", linewidth = 0.5))
gg
ggsave(paste0(outfolder, "seasons-combinations.", imgfmt), plot = gg,
       width = 10, height = 10, dpi = 500)
##################################
