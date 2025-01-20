############WGS SNPS##############
bigsnpsfile <- "rawdata/WGS-27577snps"
vcf_header <- system(str_glue('zcat "{bigsnpsfile}".vcf.gz | grep -m1 "#CHROM"'), ignore.stderr = T, intern = T)
vcf_header <- strsplit(vcf_header, split = "\t")[[1]]

bigsnps.data <- read_delim(paste0(bigsnpsfile, ".vcf.gz"), delim = "\t",
                   comment = "#", col_types = cols(.default = "c"), col_names = vcf_header)
colnames(bigsnps.data) <- make.names(vcf_header)

bigsnps.ad <- read_delim(paste0(bigsnpsfile, "_allele-depth.tsv.gz"),
                         delim="\t", col_types = cols(.default = "c"))

##################SHAPE##################
bigsnps.data <- bigsnps.data[,c(1:5)]
colnames(bigsnps.data)[1] <- "CHROM"

bigsnps.merged <- merge(bigsnps.data, bigsnps.ad, by = c("CHROM", "POS"), sort = F)
####################################

##################USING AD TO GET GENOTYPE##################
bigsnps <- bigsnps.merged
bigsnps[,-c(1:5)] <- apply(bigsnps.merged[,-c(1:5)], MARGIN = c(1,2),
                           FUN = getgenotype)
#Conversion
bigsnps[,-c(1:5)] <- apply(bigsnps[,-c(1:5)], MARGIN = 2,
                           FUN = convertbarcode, 
                           ref=as.vector(bigsnps.merged$REF), alt=as.vector(bigsnps.merged$ALT))
####################################

#################TRANSFORMING DATAFRAME#################
sample.bigsnps <- t(bigsnps)
colnames(sample.bigsnps) <- paste(sample.bigsnps[1,], 
                                  sub(pattern =" ", replacement = "", sample.bigsnps[2,]), 
                                  sep = ".")
sample.bigsnps <- as.data.frame(sample.bigsnps[-c(1,2,3),])
sample.bigsnps <- cbind(row.names(sample.bigsnps), sample.bigsnps)
colnames(sample.bigsnps)[1] = "Genome_ID"

sample.bigsnps$Barcode <- apply(X = sample.bigsnps[,-1], MARGIN = 1,
                                FUN =  function(x) {paste(x, sep="", collapse = "")})
sample.bigsnps$Study <- "BigWGS"
########################################################

sample.bigsnps$Genome_ID <- make.names(sample.bigsnps$Genome_ID)
dim(sample.bigsnps[,-c(which(colnames(sample.bigsnps)%in%c("Genome_ID","Study","Barcode")))])
##Add Ns and Xs
n.num <- apply(X = sample.bigsnps[,-c(which(colnames(sample.bigsnps)%in%c("Genome_ID","Study","Barcode")))], 
               FUN = function (x) { length(which(x=="N"))}, 
               MARGIN = 1)
x.num <- apply(X = sample.bigsnps[,-c(which(colnames(sample.bigsnps)%in%c("Genome_ID","Study","Barcode")))], 
               FUN = function (x) { length(which(x=="X"))}, 
               MARGIN = 1)
sample.bigsnps$Ns <-  n.num
sample.bigsnps$Xs <-  x.num
sort(nrow(bigsnps.data) - sample.bigsnps$Xs)
length(which(nrow(bigsnps.data) - sample.bigsnps$Xs - sample.bigsnps$Ns > 3750))
vqslodfiltsamples <- sample.bigsnps$Genome_ID[which(nrow(bigsnps.data) - sample.bigsnps$Xs > 3750)]
qualfiltsamples <- read.table(textConnection("J0010103_1412
J0020302_1506
J0030406_1701
J0030406_1702
J0030406_1705
J0030409_1610
J0050607_1412
J0050611_1611
J0050617_1611
J0060701_1511
J0070804_1610
J0101103_1701
J0101103_1702
J0101103_1703
J0111201_1412
J0111208_1610
J0131403_1511
J0131406_1512
J0141503_1511
J0151608_1701
J0151612_1511
J0171802_1611
J0171804_1610
J0171805_1611
J0171806_1511
J0181911_1611
J0242509_1611
J0242517_1611
J0262707_1412
J0262707_1506
J0262723_1511
J0262731_1412
J0272803_1512
J0272804_1511
J0272809_1611
J0272809_1612
J0272812_1611
J0272813_1610
J0272813_1611
J0272813_1612
J0272813_1701
J0272813_1702
J0272813_1703
J0272813_1705
J0272815_1611
J0272815_1701
J0272815_1702
J0272815_1703
J0282915_1511
J0282916_1701
J0282916_1702
J0282916_1703
J0282916_1705
K0010101_1611
K0010105_1506
K0022008_1512
K0030301_1504
K0030302_1512
K0030302_1610
K0030302_1611
K0030303_1512
K0030303_1605
K0030303_1612
K0050502_1512
K0050512_1611
K0060617_1512
K0060617_1603
K0070714_1611
K0080809_1702
K0090903_1512
K0090903_1605
K0101005_1512
K0101005_1701
K0101007_1412
K0111111_1611
K0141403_1504
K0141403_1511
K0161621_1506
K0171716_1612
K0171716_1701
K0171716_1702
K0171716_1703
K0171719_1701
K0171722_1512
K0171728_1512
K0171732_1512
K0171738_1412
K0181809_1512
K0181810_1705
K0202002_1512
K0212110_1512
K0222208_1611
K0222213_1612
K0222213_1701
K0232412_1701
K0252713_1512
K0252717_1512
K0252808_1611
K0252808_1701
K0262902_1511
K0262913_1611
K0262922_1612
K0262922_1704
K0262922_1705
K0262940_1511
K0262940_1512
K0262940_1703
K0262940_1704
K0263004_1511
K0263004_1512
K0263008_1612
K0263008_1701
K0263008_1704
K0263008_1705
K0283406_1512
K0323702_1611
K0344032_1512
K0364309_1511
K0364311_1511
K0364312_1607
K0364314_1701
K0364319_1511
K0364319_1612
K0364324_1701
K0374406_1612
K0374406_1701
K0374408_1412
K0374413_1611
K0374502_1412
K0374502_1612
K0384601_1702
K0384602_1512
K0384604_1511
K0384611_1702
K0394906_1611
K0394938_1611
K0405016_1611
K0415108_1701
K0415108_1702
K0425208_1512
K0425216_1511
N0050505_1610
N0131303_1701
N0131318_1611
N0131323_1610
N0141403_1612
P0010104_1611
P0010121_1611
P0010126_1611
P0010127_1610
P0010137_1610
P0020206_1701
P0020206_1702
P0020206_1703
P0020206_1704
P0020206_1705
P0020213_1612
P0020213_1701
P0020213_1705
P0020215_1701
P0030304_1607
P0050504_1611
P0050514_1611
P0050519_1705
P0060609_1611
P0060616_1611
P0060620_1611
P0060622_1611
P0060623_1612
P0060623_1702
P0060623_1703
P0060623_1704
P0060623_1705
P0070705_1611
P0070707_1610
P0070710_1607
P0070722_1612
P0070722_1701
P0070722_1702
P0080805_1612
P0080805_1705
P0080836_1612
P0080836_1701
P0080840_1611
P0080842_1612
P0080842_1703
P0080842_1704
P0080842_1705
P0080843_1611
P0080854_1611
P0080857_1702
P0080880_1612
P0111104_1611
P0111105_1701
P0111106_1611
P0121204_1701
P0121204_1703
P0121204_1704
P0121204_1705"))
table(vqslodfiltsamples%in%qualfiltsamples$V1)

sample.bigsnps <- sample.bigsnps[which(nrow(bigsnps.data) - sample.bigsnps$Xs>minwgssnps),]
############SAVING FILES############
write.csv(sample.bigsnps[-c(1,2),c("Genome_ID","Study","Barcode","Ns", "Xs")],
          file = paste0("out/", wgssnps.filename),
          quote=F, row.names = F)
####################################