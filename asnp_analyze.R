#######################################
# asnp_analyze.R
# v1.1
# 6 nov 2015
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# analyze fingerprint files to identify genotype matches.
# this script is optimized for analyzing 450k data only (using the same SNP list, in the same order, no missoing values!).
# for all other purposes use bsnp_analyze.R.
#
# parameters:
# [-[-cor|t] [<double>]]    min pairwise pearson correlation to form group/cluster. default: 0.9
# [-[-out|o] <character>]   output file basename. required.
# [-[-inp|i] <character>]   input files/directories, separated by comma. required.

options(max.print=1000)
options(stringsAsFactors=FALSE)
options(scipen=999)

suppressMessages(library(parallel))
suppressMessages(library(amap))
suppressMessages(library(getopt))


###############################################
# Implementation of simple union-find algorithm
# usage: union_find(pair_matrix) where pair_matrix
# is a set of validated pairs (in rows).
# Return value: a list of group of connected elements
# 
# by Simon Papillon-Cavanagh, 11 Apr 2014
###############################################

set_find <- function(i, target_array) {
  if(target_array[i] < 0 || i < 0)  {
    return(i) 
  } else {
    return(set_find(target_array[i], target_array))
  }
}

set_union <- function(i, j, target_array) {
  if(i == j) {
    return(target_array)
  } else {
    target_array[i] <- j
    return(target_array)
  }
}

union_find <- function(pair_matrix) {
  pair_matrix <- as.matrix(pair_matrix)
  d <- dim(pair_matrix)
  dim(pair_matrix) <- NULL
  sample_list <- unique(pair_matrix)
  num_pair_matrix <- match(pair_matrix, sample_list)
  dim(num_pair_matrix) <- d
  
  network <- rep(-1, length(sample_list))
  
  for(row in 1:nrow(num_pair_matrix)) {
    network <- set_union(set_find(num_pair_matrix[row, 1], network),
                         set_find(num_pair_matrix[row, 2], network),
                         network)
  }
  
  network <- sapply(network, set_find, network)
  
  sapply(which(network == -1), function(root) {
    children <- which(network == root)
    sample_list[c(root, children)]
  }, simplify = FALSE)
}

###############################################


### arguments
spec <- matrix(c("help", "h", 0, "logical",
                 "cor", "t", 2, "double",
                 "out", "o", 1, "character",
                 "inp", "i", 1, "character"), ncol=4, byrow=TRUE)

a <- getopt(spec)
if(!is.null(a$help)) {
  message(getopt(spec, usage=TRUE, command="bsnp_analyze.R"), appendLF=FALSE)
  q()
}
if(is.null(a$out) | is.null(a$inp)) {
  message(getopt(spec, usage=TRUE, command="bsnp_analyze.R"), appendLF=FALSE)
  message("define input and output")
  q()
} else {
  a.out <- a$out
  a.in <- strsplit(a$inp, ",")[[1]]
}
a.cor <- ifelse(is.null(a$cor), 0.9, a$cor)


### load data
# list files in directories
message("list files")
aa <- lapply(a.in, function(ai) {
  ai <- sub("/$", "", ai)
  if(length(list.files(ai)) == 0) {
    return(ai)
  } else {
    ail <- paste(ai, list.files(ai, recursive = TRUE), sep="/")
    ail <- ail[grepl(".fp$", ail)]    
    message(paste(" -", length(ail), "fp files in directory", ai))
    return(ail)
  }
})

# read files
message("read files (mc.cores=15)")
aa <- unique(sort(unlist(aa)))
fp <- mclapply(aa, function(aai) {
  read.table(aai, sep="\t", comment.char="#", colClasses=c(rep("NULL", 7), "numeric"), nrows = 1000, stringsAsFactors=FALSE)[, 1]
}, mc.cores=15)
message(" - done")


### coerce
message("coerce into matrix")
d <- do.call(rbind, fp)
rownames(d) <- sub(".fp$", "", sapply(strsplit(aa, "/"), tail, 1))
colnames(d) <- read.table(aa[1], sep="\t", comment.char="#", colClasses=c(rep("NULL", 4), "character", rep("NULL", 3)), nrows = 1000, stringsAsFactors=FALSE)[, 1]

message(paste(" -", nrow(d), "samples"))
message(paste(" -", ncol(d), "positions"))


### cluster
# calc correlation
message("calculate pairwise correlation (nbproc=45)")
d.dist <- Dist(d, method="correlation", nbproc=45)
d.cor <- matrix(NA, ncol=attr(d.dist, "Size"), nrow=attr(d.dist, "Size"))
d.cor[lower.tri(d.cor)] <- d.dist <= (1-a.cor)

d.cor.w <- apply(d.cor, 1, which)
m.cor <- cbind(rep(attr(d.dist, "Labels"), sapply(d.cor.w, length)),
               attr(d.dist, "Labels")[unlist(d.cor.w)])
message(" - ", paste(nrow(m.cor), "correlated samples"))

# form clusters
message("define clusters")
if(nrow(m.cor) > 0) {
  l <- union_find(m.cor[, c(1,2)])
  l <- sapply(l, sort, decreasing=TRUE)  # newest samples front
  lv <- sapply(l, paste, collapse=" /// ")
  lv <- sort(lv, decreasing=TRUE)
} else {
  l <- lv <- NULL
}
message(" - ", paste(length(lv), "clusters found"))


### output
message("write output")
write.table(lv, sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE, file=a.out)
message(" - done")


