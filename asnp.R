#######################################
# asnp.R
# v1.1
# 3 nov 2015
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# extract allele frequencies from 450k data for in-silico genotyping (fingerprinting).
# returns intensities for ref and alt and betavalue-like allele frequency.

options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen=999)

suppressMessages(library(minfi))


# args
a <- commandArgs(trailingOnly = TRUE)
a1 <- a[1]
#a1 <- "/cb0801/hovestad/chaos/fingerprint/fingerprint/snp141Common_RS-CG.n436.vh20151103.bed"
a2 <- a[2]
#a2 <- "/cb0803/hovestad/data450k/minfi/6752625054/6752625054_R04C02.RData"
a3 <- a[3]
#a3 <- "120618__60AB_Bender_Jones_SampleMethylationProfile_all_noAnno.MB113"
a4 <- a[4]
#a4 <- "/cb0801/hovestad/chaos/fingerprint/test/120618__60AB_Bender_Jones_SampleMethylationProfile_all_noAnno.MB113.fp"

# load SNP positions
#message("loading ", a1)
pos <- read.table(a1, sep="\t")

# load 450k data
message("loading ", a2)
load(a2)

# function
intensitySNP <- function(rg, a, b, a_col, b_col) {
  rg.green <- getGreen(rg)[, 1]
  rg.red <- getRed(rg)[, 1]
  
  a.split <- strsplit(a, ",")
  b.split <- strsplit(b, ",")
  
  a.green <- sapply(split(rg.green[unlist(a.split)], rep(seq(length(a.split)), sapply(a.split, length))), sum)
  a.red <- sapply(split(rg.red[unlist(a.split)], rep(seq(length(a.split)), sapply(a.split, length))), sum)
  b.green <- sapply(split(rg.green[unlist(b.split)], rep(seq(length(b.split)), sapply(b.split, length))), sum)
  b.red <- sapply(split(rg.red[unlist(b.split)], rep(seq(length(b.split)), sapply(b.split, length))), sum)
  
  A <- ifelse(a_col == "Green", a.green, a.red)
  B <- ifelse(b_col == "Green", b.green, b.red)
  
  cbind(A, B)
}

# process
snp <- intensitySNP(RGset, pos$V10, pos$V11, pos$V12, pos$V13)

# write
dnp.df <- data.frame(sample=a3, chr=pos$V1, pos=pos$V3, strand="*", name=pos$V4, snp, beta=round(snp[, 1]/rowSums(snp), 3))

message("writing ", a4)
write.table(c("# asnp.R v1.1 output",
              paste0("# positions: ", a1),
              paste0("# idat: ", a2)),
            file=a4, quote=FALSE, col.names=FALSE, row.names=FALSE)

write.table(dnp.df, file=a4, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

