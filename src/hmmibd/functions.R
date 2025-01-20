#All positions that are unknown or mixed are ignored
cv4ibd <- function(sample, ref, alt){
  sample <- as.vector(sample)
  sample[sample==ref] <- 1
  sample[sample==alt] <- 0
  sample[sample=="N"] <- -1
  sample[sample=="X"] <- -1
  return(sample)
}
