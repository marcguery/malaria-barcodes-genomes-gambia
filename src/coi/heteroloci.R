
nnum.snp <- read.csv("../read/out/barcodes-consensus.csv")
barcodesnps <- unique(nchar(nnum.snp$Barcode))
nnum.snp <- nnum.snp[,c("ID", "Ns", "Xs")]
nnum.snp$prop <- nnum.snp$Ns/(barcodesnps-nnum.snp$Xs)

nnum.wgs <- read.csv("../read/out/WGS-SNPs.csv")
wgssnps <- unique(nchar(nnum.wgs$Barcode))
nnum.wgs <- nnum.wgs[,c("Genome_ID", "Ns", "Xs")]
colnames(nnum.wgs) <- c("ID", "Ns", "Xs")
nnum.wgs <- nnum.wgs[wgssnps - nnum.wgs$Xs > minsnpnum,]
nnum.wgs$prop <- nnum.wgs$Ns/(wgssnps-nnum.wgs$Xs)

nnum <- merge(nnum.snp, nnum.wgs, by = "ID",
              suffixes = c(".SNP", ".WGS"), all = T)

write.csv(nnum, "out/Ns-prop.csv", quote = F, row.names = F)





