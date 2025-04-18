#!/usr/bin/env Rscript
#03/12/2024

#This script will scrutinize inbreeding levels of malaria parasite Plasmodium falciparum
#in 4 nearby villages of the eastern part of The Gambia between 2014 and 2017

#Steps of the analysis pipeline:

#1. Read and homogenize data input files of epidemiological, barcode and WGS data

#2. Format barcodes to be read by hmmIBD and calculate pairwise IBDs

#3. Build networks of samples based on IBD on barcodes improved with WGS

#4. Calculate the drug resistance marker prevalence of barcode consensus and WGS snps

#5. Estimate COI using barcode consensus and WGS SNPs

library("reshape2")
library("lubridate")
library("scales")
library("dplyr")
library("tidyr")
library("readr")
library("stringr")
library("igraph")
library("ape")
library("IRanges")
library("binom")
library("ggplot2")
library("ggnewscale")
library("ggforce")
library("ggpattern")
library("RColorBrewer")
library("cowplot")

options(stringsAsFactors = F, tz="GMT")
imgfmt <- "tiff" #image extension without the "."
setwd("./src/read")
###########################1. READ###########################
setwd("../read")
wgsbarcodes.filename <- "WGS-barcodes.csv"
wgssnps.filename <- "WGS-SNPs.csv"
consbarcodes.filename <- "barcodes-consensus.csv"
compsites.filename <- "comp-sites.csv"
minwgssnps <- 3750 #Weird cutoff for legacy purposes: 
# this includes only the 199 high quality genomes that were kept 
# with the initial filters of QUAL > 10000 (now MAF > 0.01) 
# and 4000 SNPs at least. Out of the 199 genomes, 198 had more than 4000
# SNPs while one had 3756 SNPs, hence this cutoff. When using a 'cleaner' 
# value of 3000 SNPs at least, you get 5 more genomes that have 
# between 3000 and 3750 SNPs.
minbarcodesnps <- 30
source("functions.R") #Various functions used in this step
source("WGS-barcode-mixed-cutoff.R") #Estimate WGS barcodes N cutoff
minprop <- 0.8
source("WGS-barcodes.R") #Load WGS snp for barcodes
source("merge-SNP-WGS.R") #Merge barcodes from genotyping and WGS
source("WGS-SNPs.R") #Load WGS whole good quality SNPs
#Quite long to run, comment after first run:
source("get-comp-sites.R") #Calculate the number of comparable sites between barcodes
source("figures.R")
rm(list=setdiff(ls(), c("imgfmt")))
######################################################
###########################2. HMMIBD###########################
setwd("../hmmibd")
wgssnps.filename <- "WGS-SNPs.csv"
consbarcodes.filename <- "barcodes-consensus.csv"
hmmibdwgssnps.filename <- "hmmIBD-WGS-snps.tsv"
source("functions.R") #Various functions used in this step
source("formathmm-barcodes.R") #Format barcode data for hmmIBD
source("formathmm-WGS-snps.R") #Format WGS SNP data for hmmIBD
rm(list=setdiff(ls(), c("imgfmt")))
#system("./run-hmmibd-local.sh", wait=FALSE) #uncomment to run hmmIBD

######################################################

###########################3. NETWORKS###########################
setwd("../network")
source("functions.R") #Various functions used in this step
wgs.hmmfile <- "IBD-WGS.hmm_fract.txt.zip"
mininfsites <- 100
source("WGS-network.R") #Create network from WGS SNPs
rm(list=setdiff(ls(), c("imgfmt")))
source("functions.R") #Various functions used in this step
barcode.hmmfile <- "IBD-barcodes.hmm_fract.txt.zip"
clusterIBDmin <- 0.5
minsites <- 30
mininfsites <- 10
source("barcodes-network.R") #Build network
source("network-doi-strains.R") #Estimate individual duration of infection
source("remove-continf.R") #Remove oversampling due to continuous infections
source("diversity.R")  #Get the distribution of IBD values
removecontinf <- F #Revert back to network without removal of oversampling
source("group-network.R") #Group samples by location and date
source("figures.R")
rm(list=setdiff(ls(), c("imgfmt")))
######################################################

###########################4. DRUG RES###########################
setwd("../drugres")
source("../read/functions.R") #Functions to generate a barcode from WGS SNPs
source("WGS-DR.R") #Format drug resistance markers from WGS
source("barcode-DR.R") #Format durg resistance markers from barcodes
source("DR-consensus.R") #Combine WGS and barcode data
source("figures.R") #Statistics of drug resistance prevalence
rm(list=setdiff(ls(), c("imgfmt")))
######################################################

###########################5. COI###########################
setwd("../coi")
minsnpnum <- 3750 #Weird cutoff for legacy purposes (see 1.READ section)
source("fws.R") #Calculate Fws according to Manske et al. 2012
source("heteroloci.R") #Calculate the proportion of Ns in each barcode/WGS
source("figures.R")
rm(list=setdiff(ls(), c("imgfmt")))
######################################################


