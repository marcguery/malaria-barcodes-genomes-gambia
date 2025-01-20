################WGS snps#################
bigbarcodes <- read.csv(paste0("../read/out/", wgssnps.filename), 
                        header = T)
bigsnpsfile <- "../read/rawdata/WGS-27577snps"
vcf_header <- system(str_glue('zcat "{bigsnpsfile}".vcf.gz | grep -m1 "#CHROM"'), ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]

bigsnps.data <- read_delim(paste0(bigsnpsfile, ".vcf.gz"), delim = "\t",
                           comment = "#", col_types = cols(.default = "c"), col_names = vcf_header)
colnames(bigsnps.data) <- make.names(vcf_header)
bigsnps.data <- bigsnps.data[,c(1,2,4,5)]
colnames(bigsnps.data) <- c("chrom", "pos", "REF", "ALT")
################################

#####FOR WGS SNPS#########

#These SNPs selected for their quality 
#will be compared with hmmIBD output of consensus barcodes
bcodes <- bigbarcodes$Barcode
names(bcodes) <- bigbarcodes$Genome_ID
bigbarcodes.loci <- as.data.frame(do.call("rbind", sapply(bcodes, FUN = function(x){
  res <- strsplit(x, split = "")
})))

hmmibdwgssnps <- t(bigbarcodes.loci)
hmmibdwgssnps <- cbind(bigsnps.data[,c(1,2)], hmmibdwgssnps)
hmmibdwgssnps$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdwgssnps$chrom)
hmmibdwgssnps$chrom <- sub(pattern = "_v3", replacement = "", hmmibdwgssnps$chrom)

hmmibdwgssnps$chrom <- as.integer(hmmibdwgssnps$chrom)
hmmibdwgssnps$pos <- as.integer(hmmibdwgssnps$pos)

hmmibdwgssnps[,-c(1:2)] <- data.frame(apply(X = hmmibdwgssnps[,-c(1:2)], MARGIN = 2,
                                  FUN = cv4ibd, ref=bigsnps.data$REF, alt=bigsnps.data$ALT))

hmmibdwgssnps <- hmmibdwgssnps[with(hmmibdwgssnps, order(chrom, pos)), ]
##############################

############WRITE DATA############
write.table(hmmibdwgssnps, file = paste0("out/", hmmibdwgssnps.filename),
            sep="\t", quote = F, col.names = T, row.names = F)
########################
