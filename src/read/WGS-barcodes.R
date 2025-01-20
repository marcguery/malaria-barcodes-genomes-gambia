##################INTRO##################
#Generates the WGS barcodes
####################################
##################READ##################
#gzip broken pipe expected because only the header is fetched
snpsfile <- "rawdata/WGS-101snps"
outfile <- "WGS-barcodes"
vcf_header <- system(str_glue('cat "{snpsfile}".vcf | grep -m1 "#CHROM"'), ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]
snps <- read.delim(paste0(snpsfile, ".vcf"), 
                   comment.char = "#", header = FALSE)
colnames(snps) <- make.names(vcf_header)
snps.ad <- read.table(paste0(snpsfile, "_allele-depth.tsv"), h = T)
allSNPs <- read.csv("rawdata/101-SNP-position.csv")
####################################

##################SHAPE##################
snps <- snps[,c(1:5)]
colnames(snps)[1] <- "CHROM"

snps.merged <- merge(snps, snps.ad, by = c("CHROM", "POS"))
snps.merged <- merge(snps.merged, allSNPs[,c("Chromosome", "Position")],
                     by.x = c("CHROM", "POS"),
                     by.y = c("Chromosome", "Position"),
                     all = T)
####################################
##################USING AD TO GET GENOTYPE##################
gbarcodes <- snps.merged
gbarcodes[,-c(1:5)] <- apply(gbarcodes[,-c(1:5)], c(1,2),
                             FUN = function(x){
                               if(is.na(x))
                                 return("0,0")
                               else
                                 return(x)
                             })
gbarcodes[,-c(1:5)] <- apply(gbarcodes[,-c(1:5)], MARGIN = c(1,2),
               FUN = getgenotype)
#Conversion
gbarcodes[,-c(1:5)] <- apply(gbarcodes[,-c(1:5)], MARGIN = 2,
                   FUN = convertbarcode, 
                   ref=as.vector(snps.merged$REF), alt=as.vector(snps.merged$ALT))

gbarcodes$Locus_ID <- paste(gbarcodes$CHROM, gbarcodes$POS, sep = ".")
gbarcodes$Locus_ID <- factor(gbarcodes$Locus_ID, 
                                levels = allSNPs$Locus_ID[order(allSNPs$Barcode.order)])
gbarcodes <- gbarcodes[order(gbarcodes$Locus_ID),c(6:(ncol(gbarcodes)-1))]
#################TRANSFORMING DATAFRAME#################
sample.gbarcodes <- data.frame(t(gbarcodes))
sample.gbarcodes$Barcode <- apply(X = sample.gbarcodes, MARGIN = 1,
                                  FUN =  function(x) {paste(x, sep="", collapse = "")})
sample.gbarcodes$Genome_ID <- row.names(sample.gbarcodes)

n.num <- apply(X = sample.gbarcodes[,-c(which(colnames(sample.gbarcodes)%in%c("Genome_ID","Barcode")))], 
               FUN = function (x) { length(which(x=="N"))}, 
               MARGIN = 1)
x.num <- apply(X = sample.gbarcodes[,-c(which(colnames(sample.gbarcodes)%in%c("Genome_ID","Barcode")))], 
               FUN = function (x) { length(which(x=="X"))}, 
               MARGIN = 1)
sample.gbarcodes$Ns <-  n.num
sample.gbarcodes$Xs <-  x.num

sample.gbarcodes <- sample.gbarcodes[,c("Genome_ID", "Barcode", "Ns", "Xs")]
sample.gbarcodes <- sample.gbarcodes[sample.gbarcodes$Barcode != paste0(rep("X", 101), collapse = ""),]
########################################################

############SAVING FILES############
write.csv(sample.gbarcodes, file = paste0("out/", wgsbarcodes.filename),
          quote = F, row.names = F)
####################################

