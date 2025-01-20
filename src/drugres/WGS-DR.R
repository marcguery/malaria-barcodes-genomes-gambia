#This script needs functions from ../read/functions.R

############READING##############
drugpositions <- read.table("rawdata/DR-positions.tsv", h = T)
drugphenotype <- read.table("rawdata/DR-phenotypes.tsv",h = T,sep="\t")
wgs.drugres <- read.table("rawdata/WGS-DR-annotated.csv", h = T, sep = ";", comment.char = "")
############################

############WGS DR SNPS##############
bigsnpsfile <- "rawdata/WGS-DR"
vcf_header <- system(str_glue('cat "{bigsnpsfile}".vcf | grep -m1 "#CHROM"'),
                     ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]

bigsnps.data <- read_delim(paste0(bigsnpsfile, ".vcf"), delim = "\t",
                           comment = "#", col_types = cols(.default = "c"), col_names = vcf_header)
colnames(bigsnps.data) <- make.names(vcf_header)

bigsnps.ad <- read_delim(paste0(bigsnpsfile, "_allele-depth.tsv"),
                         delim="\t", col_types = cols(.default = "c"))

##################SHAPE##################
bigsnps.data <- bigsnps.data[,c(1:5)]
colnames(bigsnps.data)[1] <- "CHROM"

bigsnps.merged <- merge(bigsnps.data, bigsnps.ad, by = c("CHROM", "POS"), sort = F)
bigsnps.merged <- bigsnps.merged[paste0(bigsnps.merged$CHROM, bigsnps.merged$POS)%in%paste0(drugpositions$Chromosome, drugpositions$Position),]
####################################



##################USING AD TO GET GENOTYPE##################
bigsnps <- bigsnps.merged
bigsnps[,-c(1:5)] <- apply(bigsnps.merged[,-c(1:5)], MARGIN = c(1,2),
                           FUN = getgenotype)
#These are not all bi-allelic calls, the max number of alts is 3
maxalt <- max(apply(bigsnps[,-c(1:5)],
                     MARGIN = c(1,2),
                     FUN = function(x){as.numeric(substr(x, 4,4))}), na.rm = T)
#Conversion
bigsnps[,-c(1:5)] <- apply(bigsnps[,-c(1:5)], MARGIN = 2,
                           FUN = convertbarcode, 
                           ref=as.vector(bigsnps.merged$REF), alt=as.vector(bigsnps.merged$ALT), maxaltn = maxalt)
####################################

#################TRANSFORMING DATAFRAME#################
sample.bigsnps <- data.frame(t(bigsnps))
colnames(sample.bigsnps) <- paste(sample.bigsnps[1,], 
                                  sub(pattern =" ", replacement = "", sample.bigsnps[2,]), 
                                  sep = ".")
sample.bigsnps <- as.data.frame(sample.bigsnps[-c(1,2,3),])
sample.bigsnps <- cbind(row.names(sample.bigsnps), sample.bigsnps)
colnames(sample.bigsnps)[1] = "Internal.Sample.ID"

sample.bigsnps$Barcode <- apply(X = sample.bigsnps[,-1], MARGIN = 1,
                                FUN =  function(x) {paste(x, sep="", collapse = "")})
sample.bigsnps$Study <- "BigWGS"

sample.bigsnps$Internal.Sample.ID <- make.names(sample.bigsnps$Internal.Sample.ID)
########################################################


############FORMAT DR POSITIONS############
wgs.drugres.summary <- wgs.drugres[,c("Chromosome", "Position", "Annotation", "Ref", "Alt", "AAref", "AAalt")]
wgs.drugres.summary$Annotation <- sub("\\)$", "" ,sub("(\\S|\\s)+\\(", "", wgs.drugres.summary$Annotation))
wgs.drugres.summary$AA <- gsub("[A-Z]", "", wgs.drugres.summary$AAref)
wgs.drugres.summary$AA <- as.numeric(sub("/(\\S|\\s)+", "", sub("(\\S|\\s)+\\(", "", wgs.drugres.summary$AA)))

wgs.drugres.summary <- merge(wgs.drugres.summary, drugphenotype, 
                             by.x = c("Annotation", "AA"), by.y = c("Gene", "AA"))

wgs.drugres.summary$Phenotype <- NA
wgs.drugres.summary$Phenotype[str_detect(wgs.drugres.summary$AAsen, substr(wgs.drugres.summary$AAalt,1,1))] <- "Sensitive"
wgs.drugres.summary$Phenotype[str_detect(wgs.drugres.summary$AAres,substr(wgs.drugres.summary$AAalt,1,1))] <- "Resistant"
wgs.drugres.summary$Phenotype[str_detect(wgs.drugres.summary$AAsen, substr(wgs.drugres.summary$AAalt,1,1)) & str_detect(wgs.drugres.summary$AAres,substr(wgs.drugres.summary$AAalt,1,1))] <- "Mixed"

wgs.drugres.summary$ID <- paste(wgs.drugres.summary$Chromosome, wgs.drugres.summary$Position, sep = ".")
########################

############CREATE DRUG RESISTANCE BARCODES############
#For each sample, create are as many barcodes as the 
# maximum number of alts
aabarcodes.bigsnps <- data.frame("Sample" = 0, "AAbarcode" = "", "Alt" = "")
for (sample in c(3:(nrow(sample.bigsnps)))){
  for (col in colnames(sample.bigsnps[,c(2:(ncol(sample.bigsnps)-2))])){
    base <- sample.bigsnps[[col]][[sample]]
    if (substr(base,1,1) == "N"){
      if (nchar(base) == 1){
        bases <- c(sample.bigsnps[[col]][1],
                   sample.bigsnps[[col]][2])
        alts <- c(1,2)
      }else{
        indexes <- as.numeric(unlist(str_split(sub("^N", "", base),",")))
        availbases <- c(sample.bigsnps[[col]][1],
                        unlist(str_split(sample.bigsnps[[col]][2],",")))
        bases <- availbases[indexes]
        alts <- c(1:length(bases))
      }
    }else{
      bases <- c(base)
      alts <- c(1)
    }
    for (i in c(1:length(bases))){
      currbase <- bases[[i]]
      aa <- unique(wgs.drugres.summary$AAref[wgs.drugres.summary$ID == col & wgs.drugres.summary$Ref == currbase])
      if (identical(aa, character(0))){
        aa <- wgs.drugres.summary$AAalt[wgs.drugres.summary$ID == col & wgs.drugres.summary$Alt == currbase]
      }
      if (identical(aa, character(0))){
        aa <- "X"
      }
      indexboollist <- aabarcodes.bigsnps$Sample == sample &
                               aabarcodes.bigsnps$Alt == paste0("Alt",i)
      if(!any(indexboollist)){
        if(i == 1){
          aabarcodes.bigsnps <- rbind(aabarcodes.bigsnps,data.frame("Sample" = sample,
                                                                    "AAbarcode" = "",
                                                                    "Alt" = paste0("Alt",i)))
          
        }else{
          prevbcode <- sub("[A-Z]$", "", aabarcodes.bigsnps$AAbarcode[aabarcodes.bigsnps$Sample==sample &
                                                      aabarcodes.bigsnps$Alt == "Alt1"])
          aabarcodes.bigsnps <- rbind(aabarcodes.bigsnps,data.frame("Sample" = sample,
                                                                    "AAbarcode" = prevbcode,
                                                                    "Alt" = paste0("Alt",i)))
        }
        indexboollist <- aabarcodes.bigsnps$Sample == sample &
          aabarcodes.bigsnps$Alt == paste0("Alt",i)
      }
      indextoupdate <- which(indexboollist)
      aabarcodes.bigsnps$AAbarcode[indextoupdate] <- paste0(aabarcodes.bigsnps$AAbarcode[indextoupdate],substr(aa,1,1))
      if (i == 1){
        indexboollist <- aabarcodes.bigsnps$Sample == sample &
          !aabarcodes.bigsnps$Alt %in% paste0("Alt",alts)
        indextoupdate <- which(indexboollist)
        aabarcodes.bigsnps$AAbarcode[indextoupdate] <- paste0(aabarcodes.bigsnps$AAbarcode[indextoupdate],substr(aa,1,1))
      }
    }
  }
}
aabarcodes.bigsnps <- aabarcodes.bigsnps[-1,]
row.names(aabarcodes.bigsnps) <- paste(sample.bigsnps$Internal.Sample.ID[aabarcodes.bigsnps$Sample],aabarcodes.bigsnps$Alt,
                                       sep = ".")
wgs.drugres.phenotype <- wgs.drugres.summary %>%
  group_by(ID)%>%
  summarise(Annotation = unique(Annotation),
            CodonPosition = unique(AA),
            Ref = unique(AAref))
wgs.drugres.phenotype$Ref <- substr(wgs.drugres.phenotype$Ref, 1, 1)
#The code needs to change if there are multiple amino acids affected per variant (with an INDEL for example)
wgs.drugres.phenotype <- cbind(wgs.drugres.phenotype,t(do.call("rbind",str_split(aabarcodes.bigsnps$AAbarcode,""))))
colnames(wgs.drugres.phenotype) <- c("ID", "Annotation", "CodonPosition", "Ref", row.names(aabarcodes.bigsnps))

########################

############DEAL WITH MIXED POSITIONS############
#Combine the different amino acids at identical positions in the same sample barcode (not between alts)
#Since amino acids are obtained from 3-bases codons, mutations at different codon positions can produce different amino acids
#The best in the future would be to develop a varif like tool to take into accout all variants from a gene before translating
wgs.drugres.phenotype <- wgs.drugres.phenotype %>%
  group_by(Annotation, Ref, CodonPosition)%>%
  summarise_at(c(2:(nrow(aabarcodes.bigsnps)+1)),
               function(x){
                 if (all(unique(x) == "X")){
                   return("X")
                 }else{
                   nonnullcall=x[x!="X"]
                   return(paste0(unique(nonnullcall), collapse = ""))
                 }
               })

wgs.drugres.phenotype <- merge(wgs.drugres.phenotype,drugphenotype,
                               by.x = c("Annotation", "CodonPosition"),
                               by.y = c("Gene", "AA"))
samplesnum <- nrow(aabarcodes.bigsnps)
wgs.drugres.phenotype <- wgs.drugres.phenotype[,c(1,2,(samplesnum+4):ncol(wgs.drugres.phenotype),4:(samplesnum+3))]
metacolsnum <- ncol(wgs.drugres.phenotype)-samplesnum
startsamcol <- metacolsnum+1
#Double calls, I remove all amino acids that are identical to ref, keeping the mutant
which(apply(wgs.drugres.phenotype[,c(startsamcol:ncol(wgs.drugres.phenotype))],c(1,2),nchar)>1)
for (i in 1:nrow(wgs.drugres.phenotype)){
  newwgs.drugres.phenotype <- sapply(X = wgs.drugres.phenotype[i,c(startsamcol:ncol(wgs.drugres.phenotype))],
                                     FUN = function(x){
                                       if (nchar(x) > 1){
                                         return(gsub(wgs.drugres.phenotype$AAsen[i],"",x))
                                       }else{
                                         return(x)
                                       }
                                     })
  wgs.drugres.phenotype[i,c(startsamcol:ncol(wgs.drugres.phenotype))] <- as.list(newwgs.drugres.phenotype)
}
#Check if no double calls remaining
#Here double-calls were only a mixture of one type of mutant and the ref allele
which(apply(wgs.drugres.phenotype[,c(startsamcol:ncol(wgs.drugres.phenotype))],c(1,2),nchar)>1)
wgs.drugres.phenotype$CodonPosition <- as.numeric(wgs.drugres.phenotype$CodonPosition)
wgs.drugres.phenotype <- wgs.drugres.phenotype[with(wgs.drugres.phenotype,order(Group,Annotation,CodonPosition)),]
########################

############FILTER SAMPLES AND ADD METADATA############
#Filter: keep samples only if they have one known amino acid
goodsamples <- which(apply(wgs.drugres.phenotype[,startsamcol:ncol(wgs.drugres.phenotype)],
                           2,
                           FUN = function(x){
                             length(x[x=="X"]) < nrow(wgs.drugres.phenotype)
                           }))+metacolsnum
wgs.drugres.phenotype <- wgs.drugres.phenotype[,c(1:metacolsnum,goodsamples)]

#Format samples
samplesavail <- data.frame(colnames(wgs.drugres.phenotype[,c(startsamcol:ncol(wgs.drugres.phenotype))]))
samplesavail$sample <- sub("\\.Alt\\S+$", "", samplesavail[,1])
colnames(samplesavail) <- c("complete_ID", "ID")

metadata <- samplesavail
metadata$individual <- substr(metadata$ID, 1, 8)
metadata$date <- as.Date(paste0(substr(metadata$ID, 10, 13),"01"),format = "%y%m%d")
metadata$village <- substr(metadata$ID,1,1)
metadata$compound <- substr(metadata$ID,1,4)
wetseason <- c(8,12)
cycle1 <- c(as.Date("2014-12-01"),
            as.Date("2015-07-01"))
cycle2 <- c(as.Date("2015-08-01"),
            as.Date("2016-07-01"))
cycle3 <- c(as.Date("2016-08-01"),
            as.Date("2017-07-01"))
metadata$season <- ifelse(month(metadata$date)>=wetseason[1],
                          "wet","dry")
metadata$cycle[metadata$date <= cycle1[[2]]] <- 1
metadata$cycle[metadata$date <= cycle2[[2]] & metadata$date >= cycle2[[1]] ] <- 2
metadata$cycle[metadata$date <= cycle3[[2]] & metadata$date >= cycle3[[1]] ] <- 3
metadata <- metadata[!is.na(metadata$complete_ID),]
length(unique(metadata$ID))
########################

############CREATE ONE DR BARCODE PER RES GROUP############
wgs.drugres.phenotype.summary <- wgs.drugres.phenotype[,c(6,7,3,4, c(metacolsnum:ncol(wgs.drugres.phenotype)))] %>%
  group_by(Group,Name)%>%
  summarise_at(c(1:(ncol(wgs.drugres.phenotype)-5)),
               function(x){
                 paste0(x, collapse = "")
               })
wgs.drugres.phenotype.summary <- melt(wgs.drugres.phenotype.summary,id.vars = c(1:4),value.name = "Drugbarcode",variable.name = "Sample")
wgs.drugres.phenotype.groups <- wgs.drugres.phenotype.summary%>% 
  group_by(Group,Name,AAsen,AAres) %>%
  reframe(Drugbarcode = unique(Drugbarcode))
wgs.drugres.phenotype.groups$Phenotype <- NA
wgs.drugres.phenotype.groups$Phenotype[wgs.drugres.phenotype.groups$Drugbarcode==wgs.drugres.phenotype.groups$AAsen] <- "Sensitive"
wgs.drugres.phenotype.groups$Phenotype[wgs.drugres.phenotype.groups$Drugbarcode!=wgs.drugres.phenotype.groups$AAres] <- "Resistant"
wgs.drugres.phenotype.groups$Phenotype[str_detect(wgs.drugres.phenotype.groups$Drugbarcode,"X")] <- NA
wgs.drugres.phenotype.summary <- merge(wgs.drugres.phenotype.summary,wgs.drugres.phenotype.groups)
wgs.drugres.phenotype.summary <- merge(wgs.drugres.phenotype.summary,metadata,
                                       by.x = "Sample", by.y = "complete_ID")
wgs.drugres.phenotype.summary$Type = "WGS"
wgs.drugres.phenotype.summary$Sample <- wgs.drugres.phenotype.summary$ID
wgs.drugres.phenotype.summary <- wgs.drugres.phenotype.summary[,-which(colnames(wgs.drugres.phenotype.summary)%in%c("complete_ID",'ID'))]
wgs.drugres.phenotype.summary <- wgs.drugres.phenotype.summary[!duplicated(wgs.drugres.phenotype.summary),]
########################

############SAVING FILES############
write.csv(sample.bigsnps, file = paste("out/WGS-DR_DNA.csv", sep=""),
          quote=F, row.names = F)
write.table(wgs.drugres.phenotype,"out/WGS-DR_AA.tsv", quote = F, row.names = F, sep = "\t")
####################################


