###################DATA###################
#get clusters within individuals from network-doi-strains
edges.identical <- edges[edges$fract_sites_IBD>=ibdcutoff,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.identical[,c("ID1", "ID2")]
edges.topology <- simplify(graph_from_data_frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
comps <- components(edges.topology)
cluster.df <- data.frame(comps$membership)
cluster.df <- cbind(row.names(cluster.df), cluster.df)
colnames(cluster.df) <- c("ID", "clusteridentical")
nodes.identical <- merge(nodes, cluster.df, all.x = T)
nodes.identical$clusteridentical[is.na(nodes.identical$clusteridentical)] <- 0

nodes.identical$clusteridentical <- sprintf(paste0("%0",max(nchar(nodes.identical$clusteridentical)),"d"), nodes.identical$clusteridentical)

continf <- nodes.identical
continf.summary <- continf %>%
  group_by(individual, clusteridentical) %>%
  summarise(earlierind = min(ID), numind = length(ID))
continf.summary <- continf.summary[continf.summary$clusteridentical != "00" & continf.summary$numind>1,]
continf <- merge(continf, continf.summary)
edges.nocontinf <- merge(edges, continf[,c("ID", "earlierind")],
                         by.x = "ID1", by.y = "ID", all.x = TRUE)
edges.nocontinf <- merge(edges.nocontinf, continf[,c("ID", "earlierind")],
                                 by.x = "ID2", by.y = "ID", all.x = TRUE, 
                         suffixes = c("1","2"))
edges.nocontinf$earlierind1[is.na(edges.nocontinf$earlierind1)] <- edges.nocontinf$ID1[is.na(edges.nocontinf$earlierind1)]
edges.nocontinf$earlierind2[is.na(edges.nocontinf$earlierind2)] <- edges.nocontinf$ID2[is.na(edges.nocontinf$earlierind2)]
edges.nocontinf$Continf <- edges.nocontinf$earlierind1!=edges.nocontinf$ID1 | edges.nocontinf$earlierind2!=edges.nocontinf$ID2
nodes.nocontinf <- nodes[nodes$ID%in%c(edges.nocontinf$ID1, edges.nocontinf$ID2),]

write.table(edges.nocontinf,
            file = "out/nocontinf/edges-nocontinf.tsv",
            quote = F, row.names = F)
write.table(nodes.nocontinf,
            file = "out/nocontinf/nodes-nocontinf.tsv",
            quote = F, row.names = F)
edgesfiltered <- edges[edges$fract_sites_IBD >= clusterIBDmin,]
edgesfiltered.nocontinf <- merge(edgesfiltered, continf[,c("ID", "earlierind")],
                                 by.x = "ID1", by.y = "ID", all.x = TRUE)
edgesfiltered.nocontinf <- merge(edgesfiltered.nocontinf, continf[,c("ID", "earlierind")],
                                 by.x = "ID2", by.y = "ID", all.x = TRUE, suffixes = c("1","2"))
edgesfiltered.nocontinf$earlierind1[is.na(edgesfiltered.nocontinf$earlierind1)] <- edgesfiltered.nocontinf$ID1[is.na(edgesfiltered.nocontinf$earlierind1)]
edgesfiltered.nocontinf$earlierind2[is.na(edgesfiltered.nocontinf$earlierind2)] <- edgesfiltered.nocontinf$ID2[is.na(edgesfiltered.nocontinf$earlierind2)]
edgesfiltered.nocontinf$Continf <- edgesfiltered.nocontinf$earlierind1!=edgesfiltered.nocontinf$ID1 | edgesfiltered.nocontinf$earlierind2!=edgesfiltered.nocontinf$ID2

nodesfiltered.nocontinf <- nodes[nodes$ID%in%c(edgesfiltered.nocontinf$ID1, edgesfiltered.nocontinf$ID2),]

write.table(edgesfiltered.nocontinf,
            file = "out/nocontinf/edges-related-nocontinf.tsv",
            quote = F, row.names = F)
# write.table(nodesfiltered.nocontinf,
#             file = "out/nocontinf/nodes-related-nocontinf.tsv",
#             quote = F, row.names = F)

relatednessscore1 <- edges.nocontinf[edges.nocontinf$date1 < as.Date("2017-01-01") & 
                                       edges.nocontinf$date2 < as.Date("2017-01-01") &
                                       edges.nocontinf$Continf==FALSE,]%>%
  group_by(ID1)%>%
  summarize(score1 = max(fract_sites_IBD))
relatednessscore2 <- edges.nocontinf[edges.nocontinf$date1 < as.Date("2017-01-01") & 
                                       edges.nocontinf$date2 < as.Date("2017-01-01") &
                                       edges.nocontinf$Continf==FALSE,]%>%
  group_by(ID2)%>%
  summarize(score2 = max(fract_sites_IBD))
relatednessscore <- merge(relatednessscore1, relatednessscore2, 
                          by.x = "ID1",by.y = "ID2", 
                          all = T)
relatednessscore$score <- pmax(relatednessscore$score1, 
                               relatednessscore$score2, na.rm = T)
relatednessscore <- merge(relatednessscore,  nodes, by.x = "ID1", by.y = "ID")


edgesfiltered.nocontinf.all <- edgesfiltered.nocontinf
unrelnodes <- unique(relatednessscore$ID1[relatednessscore$score<0.5])
unrelnodes.df <- data.frame(matrix(ncol = ncol(edgesfiltered.nocontinf.all),
                                nrow = length(unrelnodes)))
colnames(unrelnodes.df) <- colnames(edgesfiltered.nocontinf.all)
unrelnodes.df$ID2 <- edgesfiltered.nocontinf.all$ID2[1]
unrelnodes.df$ID1 <- unrelnodes
unrelnodes.df$fract_sites_IBD <- 0
unrelnodes.df$Continf <- F
edgesfiltered.nocontinf.all <- rbind(edgesfiltered.nocontinf.all,
                                     unrelnodes.df)
nodesfiltered.nocontinf.all <- nodes[nodes$ID%in%c(edgesfiltered.nocontinf.all$ID1, edgesfiltered.nocontinf.all$ID2),]
write.table(edgesfiltered.nocontinf.all,
            file = "out/nocontinf/edges-related-missingpairs-nocontinf.tsv",
            quote = F, row.names = F)
# write.table(nodesfiltered.nocontinf.all,
#             file = "out/nocontinf/nodes-related-missingpairs-nocontinf.tsv",
#             quote = F, row.names = F)
write.csv(relatednessscore, "out/unique-barcodes.csv", quote = F, row.names = F)
nrow(relatednessscore)
length(which(relatednessscore$score >= 0.5 & relatednessscore$score < 0.9))/nrow(relatednessscore)
length(which(relatednessscore$score >= 0.9))/nrow(relatednessscore)
summary(relatednessscore$score)

relatednessscore.date <- relatednessscore%>%
  group_by(date)%>%
  summarize(numobs = length(score),
            propshared = length(which(score > 0.9))/length(score))
