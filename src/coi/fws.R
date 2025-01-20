######################READING######################
rawdir <- "../read/rawdata/"
vcfname <- "WGS-27577snps.vcf"
adname <- "WGS-27577snps_allele-depth.tsv"

#gzip broken pipe expected because only the header is fetched
vcf_header <- system(str_glue('zcat "{rawdir}{vcfname}".gz | grep -m1 "#CHROM"'), ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]
snps <- read.delim(pipe(paste("zcat", paste0(rawdir, vcfname, ".gz"))), 
                           comment.char = "#", header = FALSE)
colnames(snps) <- make.names(vcf_header)
snps.ad <- read.delim(pipe(paste("zcat", paste0(rawdir, adname, ".gz"))), 
                      header = TRUE)
############################################

######################FORMATTING######################
snps <- snps[,c(1:5)]
colnames(snps)[1] <- "CHROM"

snps.merged <- merge(snps, snps.ad, by = c("CHROM", "POS"))

snps.ref <- as.data.frame(snps.merged[,-c(1:5)])
snps.ref <- as.data.frame(apply(snps.ref, 2, function(x){as.integer(sub(",\\S+", "", x))}),
                          row.names=paste(snps.merged[,1], snps.merged[,2], sep=":"))

snps.alt <- as.data.frame(snps.merged[,-c(1:5)], 
                          row.names = paste(snps.merged[,1], snps.merged[,2], sep=":"))
snps.alt <- as.data.frame(apply(snps.alt, 2, function(x){as.integer(sub("\\S+,", "", x))}),
                          row.names=paste(snps.merged[,1], snps.merged[,2], sep=":"))

############################################
#################CALCULATE FWS (MANSKE 2012)#################

lowdepth <- 5
snpsnum <- nrow(snps.merged)

snps.sum <- snps.ref+snps.alt
snps.ref <- snps.ref[,which(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))})>minsnpnum)]
snps.alt <- snps.alt[,which(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))})>minsnpnum)]

snps.sum <- snps.ref+snps.alt
snps.ref[snps.sum<=lowdepth] <- 0
snps.alt[snps.sum<=lowdepth] <- 0
snps.ref[snps.sum==0] <- NA
snps.alt[snps.sum==0] <- NA

f1 <- snps.ref/(snps.ref + snps.alt) #f1 = r1/(r1 + r2)
f2 <- snps.alt/(snps.ref + snps.alt) #f2 = r2/(r1 + r2)
h_sx <- 1 - (f1^2 + f2^2)  # H_s,x = 1 - (f1^2 + f2^2)

#Calculate heterozygosity for each SNP across all samples as 1-(F1^2 + F2^2)
snps.ref.sums <- apply(snps.ref, 1, sum, na.rm=T)
snps.alt.sums <- apply(snps.alt, 1, sum, na.rm=T)

f1.pop <- snps.ref.sums/(snps.ref.sums + snps.alt.sums) #F1 = R1/(R1 + R2)
f2.pop <- snps.alt.sums/(snps.ref.sums + snps.alt.sums) #F2 = R1/(R1 + R2)
h_px.pop <- 1 - (f1.pop^2 + f2.pop^2) # H_p,x = 1 - (f_MAF^2 + (1-f_MAF)^2)

maf <- apply(data.frame(f1.pop, f2.pop), 1, min)
write.table(maf, "out/maf.tsv", row.names = T, col.names = F, quote = F, sep ="\t")

#10 bins of SNPs with MAF values from 0 to 0.5
bins <- paste0("MAF", seq(0, 45, 5))
snps.binned <- cut(maf, right = T, include.lowest = F, 
          breaks = seq(0,0.5, 0.05),
          labels = bins)

#within-sample heterozygosity h_s,f for each bin
h_sx.binned <- bins
h_sx.binned <- lapply(h_sx.binned, function(x){apply(h_sx[which(snps.binned == x),], 2, mean, na.rm=T)})
h_sx.binned <- as.data.frame(h_sx.binned)
colnames(h_sx.binned) <- bins

#within-population heterozygosity H_p,f for each bin
h_px.pop.binned <- bins
h_px.pop.binned <- lapply(h_px.pop.binned, function(x){mean(h_px.pop[which(snps.binned == x)])})
h_px.pop.binned <- unlist(h_px.pop.binned)
names(h_px.pop.binned) <- bins

#Linear regression
h_combined <- h_sx.binned[sample(c(1:199), 10, replace = F),]
h_combined$sample <- row.names(h_combined)
h_combined <- melt(h_combined, id.vars = "sample", variable.name = "bin", value.name = "MAF")
h_combined$Hpop <- sapply(h_combined$bin, function(x){h_px.pop.binned[names(h_px.pop.binned) == x]})

gg <- ggplot(h_combined)+
  geom_point(aes(x = MAF, y = Hpop, color = sample), show.legend = F)+
  geom_smooth(aes(x = MAF, y = Hpop, color = sample, group = sample), 
              method='lm', formula = y ~ 0 + x,se = F, linetype = 2, show.legend = F)
gg

#Fws = 1 - beta
# with beta the regression coefficient (fixed intercept at 0,0) in h_s,f = beta H_p,f
fws <- 1 - apply(h_sx.binned, 1, function(x){lm(x ~ 0 + h_px.pop.binned)$coefficients[1]})
names(fws) <- row.names(h_sx.binned)

fws <- data.frame(fws)
snps.sum <- snps.ref+snps.alt
snpsUsed <- data.frame(apply(snps.sum, 2, function(x){
  length(which(x>lowdepth))}))
colnames(snpsUsed) <- "SNPs"
fws <- merge(fws, snpsUsed, by="row.names")
colnames(fws) <- c("Sample", "Fws", "SNPs")
write.csv(fws, paste("out/fws", snpsnum, lowdepth, minsnpnum, "persample.csv", sep="-"), quote = F, row.names = F)

##################################
