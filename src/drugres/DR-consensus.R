#MDR1 Y184F is not related to drug susceptibility
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-020-03506-z
#K13 C580Y absent
############COMBINE WGS and BARCODE############
drugres <- rbind(barcode.drugres,wgs.drugres.phenotype.summary)
drugres$Name <- factor(drugres$Name)
########################
drugres.cons <- drugres

drugres.cons.summary <- drugres.cons %>%
  group_by(Type, Name, AAsen, Sample)%>%
  summarise(Combinedmarker =  paste0(unique(Drugbarcode), collapse = ""))
drugres.cons.summary$ShortMarker <- drugres.cons.summary$Combinedmarker
#Remove from the combined markers calls from sensitive alleles
drugres.cons.summary$ShortMarker <- apply(drugres.cons.summary[,c("AAsen", "Combinedmarker")],
                                          1,FUN = function(x){
                                            ifelse(nchar(x[2])>1,
                                                   unique(unlist(str_split(sub(x[1],"",x[2]), ""))),
                                                   x[2])
                                            
                                          })
#Check if there are never two different resistant alleles
all(nchar(drugres.cons.summary$ShortMarker)==1)

drugres.cons.summary$Phenotype <- NA
drugres.cons.summary$Phenotype[drugres.cons.summary$ShortMarker == drugres.cons.summary$AAsen] <- "Sensitive"
drugres.cons.summary$Phenotype[drugres.cons.summary$ShortMarker != drugres.cons.summary$AAsen] <- "Resistant"
drugres.cons.summary$Phenotype[nchar(drugres.cons.summary$Combinedmarker) > 1] <- "Mixed"
drugres.cons.summary$Phenotype[drugres.cons.summary$ShortMarker == "X"] <- "Unknown"

drugres.cons.summary.stats <- drugres.cons.summary[drugres.cons.summary$ShortMarker != "X",]%>%
  group_by(Type)%>%
  summarise(ngeno = length(ShortMarker),
            nsamp = length(unique(Sample)))
drugres.cons.summary.combined.stats <- drugres.cons.summary[drugres.cons.summary$ShortMarker != "X",]%>%
  group_by(Sample)%>%
  summarise(ngeno = length(unique(Name)))
nrow(drugres.cons.summary.combined.stats)
sum(drugres.cons.summary.combined.stats$ngeno)

drugres.cons.improvedbc <- drugres.cons.summary[drugres.cons.summary$Type == "Barcode",]
drugres.cons.improvedbc <- merge(drugres.cons.improvedbc,
                                 drugres.cons.summary[drugres.cons.summary$Type == "WGS",c("Sample", "Name", "ShortMarker", "Phenotype")],
                                 by = c("Sample", "Name"),
                                 suffixes = c("", ".WGS"),
                                 all.x = T)
drugres.cons.improvedbc.stats <- drugres.cons.improvedbc[!is.na(drugres.cons.improvedbc$Phenotype.WGS) & drugres.cons.improvedbc$ShortMarker.WGS != "X",]
nrow(drugres.cons.improvedbc.stats)

drugres.cons.improvedbc$ShortMarker.WGS[drugres.cons.improvedbc$ShortMarker.WGS == "X"] <- NA
mismatches <- !is.na(drugres.cons.improvedbc$ShortMarker.WGS) & drugres.cons.improvedbc$ShortMarker.WGS != "X" &
  drugres.cons.improvedbc$ShortMarker != drugres.cons.improvedbc$ShortMarker.WGS
drugres.cons.improvedbc$ShortMarker.Cons <- NA
drugres.cons.improvedbc$Phenotype.Cons <- NA
drugres.cons.improvedbc$ShortMarker.Cons[mismatches] <- drugres.cons.improvedbc$ShortMarker.WGS[mismatches]
drugres.cons.improvedbc$Phenotype.Cons[mismatches] <- drugres.cons.improvedbc$Phenotype.WGS[mismatches]

drugres.cons.improvedbc$ShortMarker.WGS[is.na(drugres.cons.improvedbc$ShortMarker.WGS)] <- "X"
drugres.cons.improvedbc$Comparison <- NA
drugres.cons.improvedbc$Comparison[paste0(drugres.cons.improvedbc$ShortMarker,drugres.cons.improvedbc$ShortMarker.WGS) =="XX"] <- "Both unknown"
drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) & 
                                     str_detect(paste0(drugres.cons.improvedbc$ShortMarker,drugres.cons.improvedbc$ShortMarker.WGS),"X")] <- "Unknown vs. called"
drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) & 
                                     drugres.cons.improvedbc$Phenotype == "Sensitive" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Sensitive" ] <- "Match (Both sensitive)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) & 
                                     drugres.cons.improvedbc$Phenotype == "Resistant" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Resistant"] <- "Match (Both resistant)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Mixed" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Mixed" ] <- "Match (Both mixed)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Sensitive" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Mixed"] <- "Partial Match (Genome mixed, barcode sensitive)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Resistant" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Mixed"] <- "Partial Match (Genome mixed, barcode resistant)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Mixed" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Sensitive"] <- "Partial Match (Barcode mixed, genome sensitive)"

drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Mixed" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Resistant"] <- "Partial Match (Barcode mixed, genome resistant)"


drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Sensitive" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Resistant"] <- "Mismatch (Barcode sensitive, genome resistant)"
drugres.cons.improvedbc$Comparison[is.na(drugres.cons.improvedbc$Comparison) &
                                     drugres.cons.improvedbc$Phenotype == "Resistant" & 
                                     drugres.cons.improvedbc$Phenotype.WGS == "Sensitive"] <- "Mismatch (Genome sensitive, barcode resistant)"


drugres.cons.summary <- merge(drugres.cons.summary,
                              drugres.cons.improvedbc[,c("Sample", "Type", "Name","ShortMarker.Cons", "Phenotype.Cons")],
                              all = T)
drugres.cons.summary$Type[!is.na(drugres.cons.summary$ShortMarker.Cons)] <- "Consensus"
drugres.cons.summary$ShortMarker.Cons[is.na(drugres.cons.summary$ShortMarker.Cons)] <- drugres.cons.summary$ShortMarker[is.na(drugres.cons.summary$ShortMarker.Cons)]
drugres.cons.summary$Phenotype.Cons[is.na(drugres.cons.summary$Phenotype.Cons)] <- drugres.cons.summary$Phenotype[is.na(drugres.cons.summary$Phenotype.Cons)]

drugres.cons.summary <- drugres.cons.summary[!duplicated(drugres.cons.summary[,c("Sample", "Name","ShortMarker.Cons")]),]

############ADD SAMPLE METADATA############
drugres.cons.summary$individual <- substr(drugres.cons.summary$Sample, 1, 8)
drugres.cons.summary$date <- as.Date(paste0(substr(drugres.cons.summary$Sample, 10, 13),"01"),format = "%y%m%d")
drugres.cons.summary$village <- substr(drugres.cons.summary$Sample,1,1)
drugres.cons.summary$compound <- substr(drugres.cons.summary$Sample,1,4)
wetseason <- c(8,12)
cycle1 <- c(as.Date("2014-12-01"),
            as.Date("2015-07-01"))
cycle2 <- c(as.Date("2015-08-01"),
            as.Date("2016-07-01"))
cycle3 <- c(as.Date("2016-08-01"),
            as.Date("2017-07-01"))
drugres.cons.summary$season <- ifelse(month(drugres.cons.summary$date)>=wetseason[1],
                                      "wet","dry")
drugres.cons.summary$cycle[drugres.cons.summary$date <= cycle1[[2]]] <- 1
drugres.cons.summary$cycle[drugres.cons.summary$date <= cycle2[[2]] & drugres.cons.summary$date >= cycle2[[1]] ] <- 2
drugres.cons.summary$cycle[drugres.cons.summary$date <= cycle3[[2]] & drugres.cons.summary$date >= cycle3[[1]] ] <- 3


drugres.cons.summary.stats <- drugres.cons.summary[!drugres.cons.summary$ShortMarker.Cons%in%c("X"),]%>%
  group_by(Type)%>%
  summarise(ngeno = length(Sample),
            nsamp = length(unique(Sample)))
nrow(drugres.cons.improvedbc[!is.na(drugres.cons.improvedbc$Comparison) & !grepl("known", drugres.cons.improvedbc$Comparison),])
########################

############WRITE############
write.csv(drugres.cons.summary, file = "out/DR_consensus_WGS-barcodes.csv", quote = F, row.names = F)
########################
