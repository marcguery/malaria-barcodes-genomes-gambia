############READING##############
barcode.drugres <- read.table("rawdata/Barcode-DR.tsv", h = T, sep = "\t")
############################

############DOUBLE CALLS############
#Split mixed calls
altbarcodes <- barcode.drugres[str_detect(barcode.drugres$Alleles,","),]
altbarcodes$Alleles <- sub("\\S+,","", altbarcodes$Alleles)
barcode.drugres$Alleles <- sub(",\\S+","", barcode.drugres$Alleles)
barcode.drugres <- rbind(barcode.drugres,altbarcodes)
########################

############ADD DR METADATA############

phenotype.summary <- drugphenotype[,c("Group", "Name", "AAsen", "AAres")] %>%
  group_by(Group,Name)%>%
  summarise_at(c(1:2),
               function(x){
                 paste0(x, collapse = "")
               })
barcode.drugres <- merge(barcode.drugres,phenotype.summary, by.x = "Protein", by.y = "Name")
barcode.drugres$Phenotype <- NA
barcode.drugres$Phenotype[barcode.drugres$Alleles==barcode.drugres$AAsen] <- "Sensitive"
barcode.drugres$Phenotype[barcode.drugres$Alleles==barcode.drugres$AAres] <- "Resistant"
barcode.drugres$Phenotype[str_detect(barcode.drugres$Alleles,"X")] <- NA
########################


############ADD SAMPLE METADATA############
barcode.drugres$individual <- substr(barcode.drugres$Barcode_ID, 1, 8)
barcode.drugres$date <- as.Date(paste0(substr(barcode.drugres$Barcode_ID, 10, 13),"01"),format = "%y%m%d")
barcode.drugres$village <- substr(barcode.drugres$Barcode_ID,1,1)
barcode.drugres$compound <- substr(barcode.drugres$Barcode_ID,1,4)
wetseason <- c(8,12)
cycle1 <- c(as.Date("2014-12-01"),
            as.Date("2015-07-01"))
cycle2 <- c(as.Date("2015-08-01"),
            as.Date("2016-07-01"))
cycle3 <- c(as.Date("2016-08-01"),
            as.Date("2017-07-01"))
barcode.drugres$season <- ifelse(month(barcode.drugres$date)>=wetseason[1],
                          "wet","dry")
barcode.drugres$cycle[barcode.drugres$date <= cycle1[[2]]] <- 1
barcode.drugres$cycle[barcode.drugres$date <= cycle2[[2]] & barcode.drugres$date >= cycle2[[1]] ] <- 2
barcode.drugres$cycle[barcode.drugres$date <= cycle3[[2]] & barcode.drugres$date >= cycle3[[1]] ] <- 3
barcode.drugres$Type <- "Barcode"
colnames(barcode.drugres)[1:4] <- c("Name", "Sample", "Position", "Drugbarcode")
barcode.drugres <- barcode.drugres[,-c(3)]
########################

