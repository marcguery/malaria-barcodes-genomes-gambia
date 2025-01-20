###########################HMM IBD FILE###########################
databarcodes <- read_delim(file = paste0("../hmmibd/out/", barcode.hmmfile), delim = "\t")
barcodes <- read.csv("../read/out/barcodes-consensus.csv")
barcodes.comp.num <- read.csv("../read/out/barcodes-comp-sites.csv")
wgs.nodes <- read.table("out/WGS-nodes.tsv", header = TRUE)
wgs.edges <- read.table("out/WGS-edges.tsv", header = TRUE)
rescue <- T
######################################################

###########################ATTRIBUTING IDs###########################
databarcodes$Var1 <- pmin(databarcodes$sample1, databarcodes$sample2)
databarcodes$Var2 <- pmax(databarcodes$sample1, databarcodes$sample2)
databarcodes$sample1 <- databarcodes$Var1
databarcodes$sample2 <- databarcodes$Var2
databarcodes <- databarcodes[,c(1:(ncol(databarcodes)-2))]
colnames(databarcodes)[1:2] <- c("ID1", "ID2")

all(table(c(databarcodes$ID1, databarcodes$ID2))==424)

databarcodes <- merge(databarcodes,
                      barcodes.comp.num, 
                      by = c("ID1", "ID2"),
                      all.x = T)
prevcolnames <- colnames(databarcodes)
###########################ADD WGS DATA###########################
######To test different cutoff of minsites and mininfsites:######
source("pairwise-sim.R") #Calculate identity between barcodes instead of IBD
prevminsites <- minsites
prevmininfsites <- mininfsites
minsites.test <- seq(0,80,10)
mininfsites.test <- seq(0,15,5)
agreement.all <- data.frame(minsites = c(),
                       mininfsites = c(),
                       bcodepairs = c(),
                       rsquared05 = c(),
                       sens05 = c(),
                       spec05 = c(),
                       prec05 = c(),
                       agreement.kappa05 = c(),
                       rmse = c(),
                       rmsesup05 = c(),
                       rmseinf05 = c())
for (minsites.ctoff in minsites.test){
  for (mininfsites.ctoff in mininfsites.test){
    minsites <- minsites.ctoff
    mininfsites <- mininfsites.ctoff
    source("WGS-barcodes-correlation.R") #Compare IBDs between genotyped and WGS samples
    agreement.all <- rbind(agreement.all, agreement.df)
  }

}
write.table(agreement.all,
            file = "out/IBD-WGS_barcode-agreement.tsv",
            quote = F, row.names = F, sep = "\t")
minsites <- prevminsites
mininfsites <- prevmininfsites

############
source("pairwise-sim.R") #Calculate identity between barcodes instead of IBD
source("WGS-barcodes-correlation.R")

if (rescue==TRUE){
  databarcodes <- merge(databarcodes,
                        edges.both[,c("ID1", "ID2", "fract_sites_IBD.wgs", "accuracy.IBD")],
                        by = c("ID1", "ID2"),
                        all = T)
  databarcodes$N_informative_sites[is.na(databarcodes$N_comp_sites)] <- mininfsites+1
  databarcodes$N_comp_sites[is.na(databarcodes$N_comp_sites)] <- minsites+1
  
  allsamples <- expand.grid(unique(c(databarcodes$ID1,databarcodes$ID2)),
                            unique(c(databarcodes$ID1,databarcodes$ID2)), stringsAsFactors = F)
  allsamples <- allsamples[allsamples$Var1 < allsamples$Var2,]
  colnames(allsamples) <- c("ID1", "ID2")
  databarcodes <- merge(databarcodes,
                        allsamples, by = c("ID1", "ID2"),
                        all.y = T)
  databarcodes$N_informative_sites[is.na(databarcodes$N_comp_sites)] <- -1
  databarcodes$N_comp_sites[is.na(databarcodes$N_comp_sites)] <- -1
  
  databarcodes$fract_sites_IBD[databarcodes$N_comp_sites<minsites] <- -1
  databarcodes$fract_sites_IBD[databarcodes$N_informative_sites<mininfsites] <- -1
  
  databarcodes$accuracy.IBD <- as.character(databarcodes$accuracy.IBD)
  databarcodes$accuracy.IBD[is.na(databarcodes$accuracy.IBD)] <- "Unknown"
  databarcodes$fract_sites_IBD.wgs[is.na(databarcodes$fract_sites_IBD.wgs)] <- -1
  
  databarcodes$fract_sites_IBD[databarcodes$accuracy.IBD=="Unknown" & databarcodes$fract_sites_IBD.wgs > -1] <- databarcodes$fract_sites_IBD.wgs[databarcodes$accuracy.IBD=="Unknown" & databarcodes$fract_sites_IBD.wgs > -1]
  databarcodes$fract_sites_IBD[is.na(databarcodes$fract_sites_IBD)] <- databarcodes$fract_sites_IBD.wgs[is.na(databarcodes$fract_sites_IBD)]
  
  databarcodes$fract_sites_IBD[databarcodes$accuracy.IBD=="False Positive"] <- databarcodes$fract_sites_IBD.wgs[databarcodes$accuracy.IBD=="False Positive"]
  databarcodes$fract_sites_IBD[databarcodes$accuracy.IBD=="False Negative"] <- databarcodes$fract_sites_IBD.wgs[databarcodes$accuracy.IBD=="False Negative"]
  
  databarcodes <- databarcodes[,c(1:(ncol(databarcodes)-2))]
  databarcodes <- databarcodes[,order(factor(colnames(databarcodes), levels = prevcolnames))]
}else{
  databarcodes$fract_sites_IBD[databarcodes$N_comp_sites<minsites] <- -1
  databarcodes$fract_sites_IBD[databarcodes$N_informative_sites<mininfsites] <- -1
}


######################################################

###########################SPACE/TIME INFORMATION###########################
##DATES
databarcodes$date1 <- str_extract(databarcodes$ID1, pattern = "[0-9]{4}$")
databarcodes$date2 <- str_extract(databarcodes$ID2, pattern = "[0-9]{4}$")
date1 <- databarcodes$date1
date2 <- databarcodes$date2

date1[!is.na(date1)] <- paste(date1[!is.na(date1)], "01", sep="")
date2[!is.na(date2)] <- paste(date2[!is.na(date2)], "01", sep="")
date1.date <- as.Date(date1, "%y%m%d")
date2.date <- as.Date(date2, "%y%m%d")
diffdate <- round(abs(difftime(date1.date,date2.date, units="days")), digits = 2)
dateorder <- abs(as.numeric(as.factor(date1))-as.numeric(as.factor(date2)))
databarcodes$dayElapsed <- diffdate
databarcodes$groupElapsed <- dateorder
databarcodes$sameDate <- date1==date2
databarcodes$commonDate <- NA
databarcodes$commonDate[databarcodes$sameDate] <- databarcodes$date1[databarcodes$sameDate]

##LOCATIONS
databarcodes$individual1 <- sub(pattern = "_\\S+", replacement = "", x = databarcodes$ID1)
databarcodes$individual2 <- sub(pattern = "_\\S+", replacement = "", x = databarcodes$ID2)
databarcodes$sameIndividual <- databarcodes$individual1==databarcodes$individual2
databarcodes$commonIndividual <- NA
databarcodes$commonIndividual[databarcodes$sameIndividual] <- databarcodes$individual1[databarcodes$sameIndividual]


databarcodes$village1 <- substr(databarcodes$ID1, 1, 1)
databarcodes$village2 <- substr(databarcodes$ID2, 1, 1)
databarcodes$sameVillage <- databarcodes$village1==databarcodes$village2
databarcodes$commonVillage <- NA
databarcodes$commonVillage[databarcodes$sameVillage] <- databarcodes$village1[databarcodes$sameVillage]

databarcodes$compound1 <- substr(databarcodes$ID1, 1, 4)
databarcodes$compound2 <- substr(databarcodes$ID2, 1, 4)
databarcodes$sameCompound <- databarcodes$compound1==databarcodes$compound2
databarcodes$commonCompound <- NA
databarcodes$commonCompound[databarcodes$sameCompound] <- databarcodes$compound1[databarcodes$sameCompound]

#Number of samples
length(unique(c(databarcodes$ID1, databarcodes$ID2)))
#Number of individuals
length(unique(c(databarcodes$individual1, databarcodes$individual2)))
######################################################

#############NODES############
nodes <- databarcodes[!duplicated(databarcodes$ID1),c("ID1", "individual1", "date1", "village1", "compound1")]
colnames(nodes) <- c("ID2", "individual2", "date2", "village2", "compound2")
nodes <- rbind(nodes, databarcodes[!duplicated(databarcodes$ID2),c("ID2", "individual2", "date2", "village2", "compound2")])
nodes <- nodes[!duplicated(nodes$ID2),]

colnames(nodes) <- c("ID", "individual", "date", "village", "compound")

#################EDGES#################
edges <- databarcodes[,c("ID1","ID2", 
                         "fract_sites_IBD", "dayElapsed",
                     "groupElapsed", "sameIndividual", "commonIndividual", 
                     "individual1", "individual2",
                     "sameDate", "commonDate", "date1", "date2",
                     "sameVillage", "commonVillage", "village1", "village2",
                     "sameCompound", "commonCompound", "compound1", "compound2")]
edges <- edges[with(edges, order(date1, date2)),]

combos <- strsplit(paste(edges$individual1[!edges$sameIndividual], edges$individual2[!edges$sameIndividual], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonIndividual[!edges$sameIndividual] <- combos

combos <- strsplit(paste0(edges$village1[!edges$sameVillage], edges$village2[!edges$sameVillage]), split="")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonVillage[!edges$sameVillage] <- combos

combos <- strsplit(paste(edges$compound1[!edges$sameCompound], edges$compound2[!edges$sameCompound], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonCompound[!edges$sameCompound] <- combos

combos <- strsplit(paste(edges$date1[!edges$sameDate], edges$date2[!edges$sameDate], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonDate[!edges$sameDate] <- combos
##################################

#############FILTER EDGES############
edges.filtered <- edges[edges$fract_sites_IBD>=clusterIBDmin,]
# Get rid of loops and ensure right naming of vertices
edges.topology <- edges.filtered[,c("ID1", "ID2")]
edges.topology <- simplify(graph_from_data_frame(edges.topology[order(edges.topology[[1]]),],directed = FALSE))

# Find all components
comps <- components(edges.topology)
cluster.df <- data.frame(comps$membership)
cluster.df <- cbind(row.names(cluster.df), cluster.df)
colnames(cluster.df) <- c("ID", "cluster")
nodes <- merge(nodes, cluster.df, all.x = T)
nodes$cluster[is.na(nodes$cluster)] <- 0

nodes$cluster <- sprintf(paste0("%0",max(nchar(nodes$cluster)),"d"), nodes$cluster)

#################EDGES#################
edges <- merge(edges, nodes[,c("ID", "cluster")],
               by.x = "ID1", by.y = "ID", all.x = TRUE)
edges <- merge(edges, nodes[,c("ID", "cluster")],
               by.x = "ID2", by.y = "ID", all.x = TRUE,
               suffixes = c("1", "2"))

edges$sameCluster <- edges$cluster1 == edges$cluster2
edges$commonCluster <- NA
edges$commonCluster[edges$sameCluster] <- edges$cluster1[edges$sameCluster]

combos <- strsplit(paste(edges$cluster1[!edges$sameCluster], edges$cluster2[!edges$sameCluster], sep =";"), split=";")
combos <- sapply(combos, sort)
if (length(combos) > 0){
  combos <- paste0(combos[1,], combos[2,])
  edges$commonCluster[!edges$sameCluster] <- combos
}
##################################

nodes$date <- as.Date(paste0(nodes$date, "01"), "%y%m%d")
edges$date1 <- as.Date(paste0(edges$date1, "01"), "%y%m%d")
edges$date2 <- as.Date(paste0(edges$date2, "01"), "%y%m%d")

nodes$season <- apply(nodes, 1, 
                      FUN = function(x){
                        seasons$type[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})
nodes$cycle <- apply(nodes, 1, 
                     FUN = function(x){
                       seasons$cycle[seasons$date1 <= x [3] & seasons$date2 >= x[3]]})


edges <- merge(edges, nodes[,c("ID", "season", "cycle")], all.x = T,
                 by.x = "ID1", by.y = "ID")
edges <- merge(edges, nodes[,c("ID", "season", "cycle")], all.x = T,
                 by.x = "ID2", by.y = "ID", suffixes = c("1", "2"))
  

##Seasons
edges$sameSeason <- edges$season1==edges$season2
edges$commonSeason <- NA
edges$commonSeason[edges$sameSeason] <- edges$season1[edges$sameSeason]

combos <- strsplit(paste(edges$season1[!edges$sameSeason], edges$season2[!edges$sameSeason], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonSeason[!edges$sameSeason] <- combos

##Cycles
edges$sameCycle <- edges$cycle1==edges$cycle2
edges$commonCycle <- NA
edges$commonCycle[edges$sameCycle] <- edges$cycle1[edges$sameCycle]

combos <- strsplit(paste(edges$cycle1[!edges$sameCycle], edges$cycle2[!edges$sameCycle], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
edges$commonCycle[!edges$sameCycle] <- combos

edges$cycleElapsed <- abs(edges$cycle2 - edges$cycle1)

edges.filtered <- edges[edges$fract_sites_IBD>=clusterIBDmin,]

###########################SAVE DATA###########################
write.table(nodes,
            file = "out/nodes.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(edges,
            file = "out/edges.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(edges.filtered,
            file = "out/edges-related.tsv",
            quote = F, row.names = F, sep = "\t")
##################################