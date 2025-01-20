###############SNP and WGS BARCODES###############
WGSbarcodes <- read.csv(paste0("out/", wgsbarcodes.filename))
SNPbarcodes <- read.csv("rawdata/101-genotyped-SNPS.csv")
allSNPs <- read.csv("rawdata/101-SNP-position.csv")
##############################

###############MERGING BARCODES###############
SNPbarcodes$ID <- SNPbarcodes$Barcode_ID
SNPbarcodes$Study <- "Molecular"
SNP.missing <- apply(do.call("rbind", strsplit(SNPbarcodes$Barcode, split = "")), 2, FUN = function(x){
  length(which(x=="X"))
})/nrow(SNPbarcodes)
which(SNP.missing == 1) #Missing loci in molecular barcodes

WGSbarcodes$ID <- WGSbarcodes$Genome_ID
WGSbarcodes$Study <- "WGS"
WGS.missing <- apply(do.call("rbind", strsplit(WGSbarcodes$Barcode, split = "")), 2, FUN = function(x){
  length(which(x=="X"))
})/nrow(WGSbarcodes)
which(WGS.missing == 1) #Missing loci in WGS barcodes

barcodes <- rbind(WGSbarcodes[,c("ID", "Barcode", "Study")], SNPbarcodes[, c("ID", "Barcode", "Study")])
barcodes$Date <- as.numeric(sub("\\S+_", "", barcodes$ID))
##############################

##########CORRELATE WGS and SNP barcodes##########

#Samples having a SNPs and WGS barcode
dups <- names(which(table(barcodes$ID)==2))

#Code for each site comparison
#This will trigger a warning as there are expected NAs in the scores

bcodescores <- sapply(dups,
                      FUN = function(x) { 
                        snpbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study != "WGS" & barcodes$ID==x], ""))
                        wgsbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study == "WGS" & barcodes$ID==x], ""))
                        positionqual(snpbcode, wgsbcode) })
row.names(bcodescores) <- c(1:nrow(bcodescores))
#Removal of SNPs not present in WGS
bcodescores <- bcodescores[-which(WGS.missing==1),]
bcodescores.melted <- melt(bcodescores, varnames = c("SNP", "Sample"), value.name = "Score")

#Summary of the comparison of each SNP
bcodescores.summary <- statfromscore(bcodescores, 1)
dateperbarcode <- as.numeric(sub("\\S+_", "", colnames(bcodescores)))
bcodescores.summary.before1605 <- statfromscore(bcodescores[,dateperbarcode <= 1605], 1)
bcodescores.summary.after1605 <- statfromscore(bcodescores[,dateperbarcode > 1605], 1)

#Summary of the comparison of each SNP in each barcode
bcodescores.summaryperbarcode <- statfromscore(bcodescores, 2)

#PLOT ORDER based on SNP quality
#Order of SNPs for plotting: 
##least number of individuals with different SNPs --> most
##least number of individuals with unkowns SNPs --> most
orderedsnps <- rownames(bcodescores.summary[with(bcodescores.summary, order(-different, -unknown)),])
#Order of barcode comparison for plotting
##past --> present
##least number of different SNPs --> most
##least number of unkowns SNPs --> most
datetile <- sub("\\S+_", "", bcodescores.melted$Sample)
bcodescores.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.melted$date <- as.numeric(datetile)
orderedbarcodes <- colnames(bcodescores[,with(bcodescores.summaryperbarcode, order(-date, -different, -unknown))])

bcodescores.summary$Barcode.order <- row.names(bcodescores.summary)
bcodescores.summary.before1605$Barcode.order <- row.names(bcodescores.summary.before1605)
bcodescores.summary.after1605$Barcode.order <- row.names(bcodescores.summary.after1605)
snps.raw <- merge(allSNPs, bcodescores.summary, all.x=T, by = "Barcode.order")
snps.raw.before1605 <- merge(allSNPs, bcodescores.summary.before1605, all.x=T, by = "Barcode.order")
snps.raw.after1605 <- merge(allSNPs, bcodescores.summary.after1605, all.x=T, by = "Barcode.order")

#21 bad SNPs
badsnps <- snps.raw$Barcode.order[!is.na(snps.raw$Agreement) & snps.raw$Agreement < 0.8]

#Resetting 21 bad SNPs to 'X' in molecular barcodes
barcodes.raw <- barcodes
barcodes$Barcode[barcodes$Study != "WGS" & barcodes$Date > 1605] <- sapply(barcodes$Barcode[barcodes$Study != "WGS" & barcodes$Date > 1605],
                             FUN = function(x){
                               bcode <- unlist(strsplit(x, split = ""))
                               bcode[badsnps] <- "X"
                               return(paste0(bcode, collapse = ""))
                             })
####################

#############COMPARING BARCODES WITHOUT BAD SNPS#############
#ESTIMATE THE CUTOFF FOR N in WGS barcodes
#Barcodes without bad SNPs
#This will trigger a warning as there are expected NAs in the scores
bcodescores.goodsnps <- sapply(dups,
                               FUN = function(x) { 
                                 snpbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study != "WGS" & barcodes$ID==x], ""))
                                 wgsbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study == "WGS" & barcodes$ID==x], ""))
                                 positionqual(snpbcode, wgsbcode) })
row.names(bcodescores.goodsnps) <- c(1:nrow(bcodescores.goodsnps))
#Removal of SNPs not present in WGS
bcodescores.goodsnps <- bcodescores.goodsnps[-which(WGS.missing==1),]
bcodescores.goodsnps.melted <- melt(bcodescores.goodsnps, varnames = c("SNP", "Sample"), value.name = "Score")

#Summary of the comparison of each SNP
bcodescores.goodsnps.summary <- statfromscore(bcodescores.goodsnps, 1)

#Summary of the comparison of each SNP in each barcode
bcodescores.goodsnps.summaryperbarcode <- statfromscore(bcodescores.goodsnps, 2)

#Adding dates
dateperbarcode2 <- sub("\\S+_", "", colnames(bcodescores))
datetile <- sub("\\S+_", "", bcodescores.goodsnps.melted$Sample)
bcodescores.goodsnps.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.goodsnps.melted$date <- as.numeric(datetile)

bcodescores.goodsnps.summary$Barcode.order <- row.names(bcodescores.goodsnps.summary)
snps.goodsnps <- merge(allSNPs, bcodescores.goodsnps.summary, all.x=T, by = "Barcode.order")

#This dataframe can be used to determine the appropriate cutoff
#of the minimal allele frequency of the majority base called
#to discard the minority base called
#(i. e., to add mixed calls in WGS barcodes)
bcodescores.goodsnps.summaryperbarcode$ID <- colnames(bcodescores.goodsnps)
nnumber <- melt(bcodescores.goodsnps.summaryperbarcode[,c("ID", "same", "different", "unknown", "nn", "nl")])
if (!exists("minprop")){
  minprop <- ""
}
colnames(nnumber) <- c("ID", "alignment", paste0("ratio", minprop))
############################

##############MAKING CONSENSUS BARCODES###############
#Consensus barcodes are the same as original barcodes with:
#Barcodes from WGS replace all unknown calls of genotyped barcodes
#Sites with a discrepancy between SNP and WGS calls are set to unkown call

#Updating barcodes with consensus barcodes, removing duplicated samples from
#SNP genotyping and WGS
barcodes.goodsnps <- barcodes
consbarcodes <- sapply(dups,
                       FUN = function(x) { 
                         snpbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study != "WGS" & barcodes$ID==x], ""))
                         wgsbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study == "WGS" & barcodes$ID==x], ""))
                         make_consensus(snpbcode, wgsbcode)})
consbarcodes.df <- data.frame(ID = names(consbarcodes), 
                              Barcode = consbarcodes,
                              Study = "Consensus",
                              Date = as.numeric(sub("\\S+_", "", names(consbarcodes))))
barcodes <- barcodes[!(barcodes$ID%in%dups & barcodes$Study == "Molecular"),]
barcodes <- rbind(barcodes, consbarcodes.df)

#This will trigger a warning as there are expected NAs in the scores
bcodescores.cons <- sapply(dups,
                           FUN = function(x) { 
                             snpbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study != "WGS" & barcodes$ID==x], ""))
                             wgsbcode <- unlist(strsplit(barcodes$Barcode[barcodes$Study == "WGS" & barcodes$ID==x], ""))
                             positionqual(snpbcode, wgsbcode) })
row.names(bcodescores.cons) <- c(1:nrow(bcodescores.cons))
#Removal of SNPs not present in WGS
bcodescores.cons <- bcodescores.cons[-which(WGS.missing==1),]
bcodescores.cons.melted <- melt(bcodescores.cons, varnames = c("SNP", "Sample"), value.name = "Score")

#Summary of the comparison of each SNP
bcodescores.cons.summary <- statfromscore(bcodescores.cons, 1)

#Summary of the comparison of each SNP in each barcode
bcodescores.cons.summaryperbarcode <- statfromscore(bcodescores.cons, 2)

#Adding dates
datetile <- sub("\\S+_", "", bcodescores.cons.melted$Sample)
bcodescores.cons.summaryperbarcode$date <- as.numeric(dateperbarcode)
bcodescores.cons.melted$date <- as.numeric(datetile)
#Snp discordance between genotyped SNP and WGS SNP
bcodescores.cons.summary$Barcode.order <- row.names(bcodescores.cons.summary)
snps.cons <- merge(allSNPs, bcodescores.cons.summary, all.x=T, by = "Barcode.order")


#Removing WGS barcodes having a new consensus
barcodes <- barcodes[-which(barcodes$Study=="WGS" & barcodes$ID%in%dups),]
#Removing SNPs absent in WGS from consensus barcodes
barcodes$Barcode <- sapply(barcodes$Barcode,
                           FUN = function(x){
                             bcode <- unlist(strsplit(x, split = ""))
                             bcode <- bcode[-which(WGS.missing==1)]
                             bcode <- paste0(bcode, collapse = "")
                           })

##############FILTERING LOW QUALITY BARCODES##############
#Add Ns and Xs
n.num <- sapply(X = barcodes$Barcode, 
               FUN = function (x) { 
                 bcode <- unlist(strsplit(x, split = ""))
                 length(which(bcode=="N"))})
x.num <- sapply(X = barcodes$Barcode, 
               FUN = function (x) { 
                 bcode <- unlist(strsplit(x, split = ""))
                 length(which(bcode=="X"))})
barcodes$Ns <-  n.num
barcodes$Xs <-  x.num

goodsnps <- nchar(barcodes$Barcode) - barcodes$Xs - barcodes$Ns
#We need to get rid of barcodes having too much unknown SNPs
ctoffsnps <- 30 #Minimal number of SNPs for a barcode to be kept
#Distribution of number of good quality SNPs in barcodes
quantile(goodsnps)
#number of barcodes having no good SNPs at all
length(which(goodsnps<ctoffsnps))

barcodes <- barcodes[goodsnps >= ctoffsnps,]
############################
##############SEASONS##############

#dry (january to july) wet (beginning august to end december)
wetseason <- c(8,12)
wetseason <- as.character(wetseason)

wetseason[nchar(wetseason)==1] <- paste0("0", wetseason[nchar(wetseason)==1])

years <- sort(unique(year(as.Date(str_sub(barcodes$Date, 1, 2),"%y"))))
wetseason <- as.Date(c(paste0(years[1], "-01-01"),
                       unlist(lapply(years, function(x){
                         paste(x, wetseason, "15", sep="-")
                       })),
                       paste0(years[length(years)], "-12-31")), "%Y-%m-%d")

seasons <- data.frame(wetseason[-length(wetseason)], wetseason[-1])
seasons$type <- rep(c("dry", "wet"), length.out=nrow(seasons))
colnames(seasons) <- c("date1", "date2", "type")
mindate <- min(as.Date(paste0(barcodes$Date, "01"), format = "%y%m%d"))
#maxdate <- max(plotdataset$Date)
maxdate <- as.Date("2017-01-01")
#Time zone warning but it is ok
seasons <- seasons[seasons$date2>mindate & seasons$date1 < maxdate,]

yearspan <- data.frame(as.Date("2014-01-01")+months(c(0,12,24,36)))
colnames(yearspan) <- "year"
yearspan$manualposition <- c(9, 6, 6, 3)
############################

##############SAVING FILES##############
write.csv(barcodes, file = paste0("out/", consbarcodes.filename), quote = F, row.names = F)
write.csv(snps, file = "out/consensus-SNPs.csv", quote = F, row.names = F)
####################################
  