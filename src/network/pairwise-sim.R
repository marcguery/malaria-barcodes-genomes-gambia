#####################PAIRWISE COMPARISON#####################
barcodes.dnabin <- barcodes$Barcode
names(barcodes.dnabin) <- barcodes$ID
barcodes.dnabin <- gsub("N", "X", barcodes.dnabin)
barcodes.dnabin <- as.DNAbin(sapply(tolower(barcodes.dnabin),
                       FUN = function(x){str_split(x,"")}))
pairwise.edges <- melt(dist.dna(barcodes.dnabin, model = "raw",
                          pairwise.deletion = T,
                          as.matrix = T))
pairwise.edges$Var1 <- as.character(pairwise.edges$Var1)
pairwise.edges$Var2 <- as.character(pairwise.edges$Var2)
pairwise.edges <- pairwise.edges[duplicated(paste0(pmin(pairwise.edges$Var1, pairwise.edges$Var2),
                                 pmax(pairwise.edges$Var1, pairwise.edges$Var2))),]
pairwise.edges$ID1 <- pmin(as.character(pairwise.edges$Var1), as.character(pairwise.edges$Var2))
pairwise.edges$ID2 <- pmax(as.character(pairwise.edges$Var1), as.character(pairwise.edges$Var2))
pairwise.edges <- pairwise.edges[,c(4,5,3)]
colnames(pairwise.edges) <- c("ID1", "ID2", "pairwise")
pairwise.edges <- merge(pairwise.edges,
                        barcodes.comp.num)
##########################################
#####################SAVE DATA#####################
write.table(pairwise.edges,"out/IBDcorr/barcode-pairwise.tsv",quote = F,row.names = F)
##########################################
