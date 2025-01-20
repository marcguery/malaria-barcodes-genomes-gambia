##################OPTIONS##################
lowdepth <- 6 #if depth is inferior strict to NUM, SNP is 'X'
if (!exists("minprop")){
  minprop <- 0.8
} #if prop is inferior strict to NUM, SNP is 'N'
#minprop is set up from the pipeline calling this script
####################################
##################FUNCTIONS##################
parsedepth <- function(ad){
  intlist <- as.integer(unlist(strsplit(ad, ",")))
  res <- sum(intlist)
  return(c(intlist[1], list(intlist[-1]), res))
}

getgenotype <- function(ad){
  depth <- parsedepth(ad)
  refdepth <- depth[1][[1]]
  altdepth <- unlist(depth[2])
  totdepth <- depth[3][[1]]
  if (max(refdepth, max(altdepth)) < lowdepth){
    return("X")
  }
  if (max(refdepth, max(altdepth))/totdepth < minprop){
    if (length(altdepth) > 1 ){
      shared <- paste(which(c(refdepth, c(altdepth))>=lowdepth),collapse = ",")
      return(paste0("N",shared))
    }else{
    return("N")}
  }
  if (refdepth >= max(altdepth)){
    return('ref')
  }else{
    if (length(altdepth) > 1 ){
      return(paste0('alt', which.max(altdepth)))
    }
    return('alt')
  }
}

convertbarcode <- function(barcode, ref, alt, maxaltn=1){
  barcode[barcode=="ref"] <- ref[barcode=="ref"]
  barcode[barcode=="alt"] <- alt[barcode=="alt"]
  if (maxaltn > 1){
    for (i in c(1:maxaltn)){
      barcode[barcode==paste0("alt", i)] <- unlist(strsplit(alt[barcode==paste0("alt", i)], ","))[i]
    }
  }
  return(barcode)
}


#Attribute a score of alignment between SNP and WGS barcodes from the same blood sample
positionqual <- function(bcode1, bcode2){
  bcode <- rbind(bcode1, bcode2)
  bcode[bcode=="A"] <- 2
  bcode[bcode=="T"] <- 4
  bcode[bcode=="C"] <- 6
  bcode[bcode=="G"] <- 8
  bcode[bcode=="N"] <- 0
  bcode[bcode=="X"] <- "NA"
  scores <- as.integer(bcode[1,])/as.integer(bcode[2,])
  scores[scores==Inf|scores==0] <- -2 #N-ATGC : e. g. 2/0 or 0/2 
  scores[scores>0 & scores!=1 & scores<Inf] <- -1 #Mismatch : e. g. 6/2, 2/6
  scores[is.nan(scores)] <- 2 #N-N : 0/0=nan
  scores[is.na(scores)] <- 0 #X-ATGCNX : e.g. 2/NA, NA/0 = NA
  return(scores)
}

statfromscore <- function(scores, side = 1){
  scores.summary <- as.data.frame(
    matrix(c(apply(scores, side,
                   FUN = function(x) {length(which(x==1))/length(which(x!=0))}),
             apply(scores, side,
                   FUN = function(x) {length(which(x==1))}),
             apply(scores, side,
                   FUN = function(x) {length(which(x==-1))}),
             apply(scores, side,
                   FUN = function(x) {length(which(x==0))}),
             apply(scores, side,
                   FUN = function(x) {length(which(x==2))}),
             apply(scores, side,
                   FUN = function(x) {length(which(x==-2))})
    ),
    ncol = 6
    ),
    row.names = row.names(scores)
  )
  colnames(scores.summary) <- c("Agreement", "same", "different", "unknown", "nn", "nl")
  return(scores.summary)
}

make_consensus <- function(molbarcode, genbarcode){
  #Replace SNP Xs by WGS call
  cbarcode <- molbarcode
  cbarcode[which(cbarcode=="X")] <- genbarcode[which(cbarcode=="X")]
  cbarcode[which(cbarcode=="N")] <- genbarcode[which(cbarcode=="N")]
  cbarcode[which(cbarcode!="X" & cbarcode!="N" & genbarcode!="X" & genbarcode!="N" & cbarcode != genbarcode)] <- "X"
  res <- as.vector(t(cbarcode))
  names(res) <- colnames(cbarcode)
  return(paste0(res, collapse = ""))
}

####################################
