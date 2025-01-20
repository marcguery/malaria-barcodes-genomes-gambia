###########################HMM IBD FILE###########################
databarcodes <- read_delim(file = paste0("../hmmibd/out/", wgs.hmmfile), delim = "\t")
databarcodes$Var1 <- pmin(databarcodes$sample1, databarcodes$sample2)
databarcodes$Var2 <- pmax(databarcodes$sample1, databarcodes$sample2)
databarcodes$sample1 <- databarcodes$Var1
databarcodes$sample2 <- databarcodes$Var2
databarcodes <- databarcodes[,c(1:(ncol(databarcodes)-2))]
colnames(databarcodes)[1:2] <- c("ID1", "ID2")
all(table(c(databarcodes$ID1, databarcodes$ID2))==198)

wgs.comp.num <- read.csv("../read/out/WGS-comp-sites.csv")
databarcodes <- merge(databarcodes,
                      wgs.comp.num, 
                      by = c("ID1", "ID2"),
                      all.x = T)
######################################################

###########################FILTERING###########################
databarcodes$fract_sites_IBD[databarcodes$N_informative_sites<mininfsites] <- -1
summary(databarcodes$N_comp_sites)
summary(databarcodes$N_comp_sites[databarcodes$fract_sites_IBD > -1])
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

databarcodes.all <- databarcodes
######################################################

#############NODES############
nodes <- databarcodes[!duplicated(databarcodes$ID1),c("ID1", "individual1", "date1", "village1", "compound1")]
colnames(nodes) <- c("ID2", "individual2", "date2", "village2", "compound2")
nodes <- rbind(nodes, databarcodes[!duplicated(databarcodes$ID2),c("ID2", "individual2", "date2", "village2", "compound2")])
nodes <- nodes[!duplicated(nodes$ID2),]

colnames(nodes) <- c("ID", "individual", "date", "village", "compound")
##########################
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

###########################SAVE DATA###########################
write.table(nodes,
            file = "out/WGS-nodes.tsv",
            quote = F, row.names = F)
write.table(edges,
            file = "out/WGS-edges.tsv",
            quote = F, row.names = F)
##################################
