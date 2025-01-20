##############NUMBER OF COMPARABLE SITES##############
for (dataset in c("WGS", "barcodes")){
  if (dataset == "WGS"){
    bcodes <- sample.bigsnps$Barcode[-c(1,2)]
    ids <- sample.bigsnps$Genome_ID[-c(1,2)]
  }else{
    
    bcodes <- barcodes$Barcode
    ids <- barcodes$ID
    }
  
    
    bcodes.num <- gsub("X", "0", bcodes)
    bcodes.num <- gsub("N", "0", bcodes.num)
    bcodes.num <- gsub("(A|T|G|C)", "1", bcodes.num)
    names(bcodes.num) <- ids
    bcodes.num <- melt(bcodes.num)
    bcodes.comp.num <- expand.grid(row.names(bcodes.num), row.names(bcodes.num))
    bcodes.comp.num$ID1 <- pmin(as.character(bcodes.comp.num$Var1), as.character(bcodes.comp.num$Var2))
    bcodes.comp.num$ID2 <- pmax(as.character(bcodes.comp.num$Var1), as.character(bcodes.comp.num$Var2))
    bcodes.comp.num <- bcodes.comp.num[,c(3,4)]
    bcodes.comp.num <- bcodes.comp.num[!duplicated(bcodes.comp.num),]
    bcodes.comp.num <- bcodes.comp.num[bcodes.comp.num$ID1!=bcodes.comp.num$ID2,]
    bcodes.comp.num$N_comp_sites <- apply(bcodes.comp.num,1,
                                            FUN = function(x){
                                              bc1 <- as.numeric(do.call("rbind", str_split(bcodes.num[row.names(bcodes.num)==x[1],1], "")))
                                              bc2 <- as.numeric(do.call("rbind", str_split(bcodes.num[row.names(bcodes.num)==x[2],1], "")))
                                              return(length(which(bc1+bc2==2)))
                                            })
    write.csv(bcodes.comp.num, paste0("out/", dataset, "-", compsites.filename), 
              quote = F, row.names = F)
}
############################
