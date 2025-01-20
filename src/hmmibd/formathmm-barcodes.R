barcodes <- read.csv("../read/out/barcodes-consensus.csv")[,c(1,2)]
allsnps <- read.csv("../read/out/consensus-SNPs.csv")
colnames(allsnps) <- c("chrom", "pos", "id", "ref", "alt")
######################CONVERT FORMAT FOR HMMIBD###############


#FOR BARCODES
bcodes <- barcodes$Barcode
names(bcodes) <- barcodes$ID
bcodes.loci <- as.data.frame(do.call("rbind", sapply(bcodes, FUN = function(x){
  res <- strsplit(x, split = "")
})))

hmmibdbcode <- t(bcodes.loci)
hmmibdbcode <- cbind(allsnps[,c(1,2)], hmmibdbcode)
hmmibdbcode$chrom <- sub(pattern = "Pf3D7_(0)?", replacement = "", hmmibdbcode$chrom)
hmmibdbcode$chrom <- sub(pattern = "_v3", replacement = "", hmmibdbcode$chrom)

hmmibdbcode$chrom <- as.integer(hmmibdbcode$chrom)
hmmibdbcode$pos <- as.integer(hmmibdbcode$pos)

hmmibdbcode[,-c(1:2)] <- data.frame(apply(X = hmmibdbcode[,-c(1:2)], MARGIN = 2,
                                            FUN = cv4ibd, ref=allsnps$ref, alt=allsnps$alt))

hmmibdbcode <- hmmibdbcode[with(hmmibdbcode, order(chrom, pos)), ]
########################

############WRITE DATA############
write.table(hmmibdbcode, file = "out/hmmIBD-barcodes.tsv", sep="\t", quote = F, col.names = T, row.names = F)
########################
