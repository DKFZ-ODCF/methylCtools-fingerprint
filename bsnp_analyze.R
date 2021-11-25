#######################################
# bsnp_analyze.R
# v1.1
# 4 nov 2015
#
# volker hovestadt
# german cancer research center
# v.hovestadt@dkfz.de
#
#
# analyze fingerprint files to identify genotype matches.
# returns groups of samples that share genotypes (_cluster.txt), pairwise correlation (_cor.txt)
# and pairwise number SNPs for which information is available (_pairobs.txt)
#
# parameters:
# [-[-cov|c] [<integer>]]   min coverage to consider SNP position. default: 10
# [-[-obs|p] [<integer>]]   min number of SNPs with sufficient coverage to calculate pairwise correlation. default: 15
# [-[-cor|t] [<double>]]    min pairwise pearson correlation to form group/cluster. default: 0.9
# [-[-rg|g] [<logical>]]    use read group information (e.g. split sample by sequencing lane). default: FALSE
# [-[-out|o] <character>]   output file basename. required.
# [-[-inp|i] <character>]   input files/directories, separated by comma. required.

options(max.print=1000)
options(stringsAsFactors=FALSE)
options(scipen=999)


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
library(getopt)
spec <- matrix(c("help", "h", 0, "logical",
                 "cor", "t", 2, "double",
                 "obs", "p", 2, "integer",
                 "cov", "c", 2, "integer",
                 "rg", "g", 2, "logical",
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
a.obs <- ifelse(is.null(a$obs), 15, a$obs)
a.cov <- ifelse(is.null(a$cov), 10, a$cov)
a.rg <- ifelse(is.null(a$rg), FALSE, TRUE)


### load data

# list files in directories
message("list files")
aa <- lapply(a.in, function(ai) {
  ai <- sub("/$", "", ai)
  if(length(list.files(ai)) == 0) {
    return(ai)
  } else {
    message(paste(" -", length(list.files(ai, recursive = TRUE)), "files in directory", ai))
    return(paste(ai, list.files(ai, recursive = TRUE), sep="/"))
  }
})

# read files
message("read files")
fp <- list()
for(aai in unique(sort(unlist(aa)))) {
  #message(aai)
  if(!file.exists(aai)) {
    message(paste("does not exist:", aai))
  } else if(grepl(".fp$", aai)) {
    aaid <- read.table(aai, sep="\t", comment.char="#", col.names=c("sample", "chr", "pos", "strand", "id", "A", "B", "ratio"))
    aaid <- aaid[(aaid$A + aaid$B) >= a.cov, ]
    if(nrow(aaid) < a.obs) {
      message(paste(" - too few observations:", aai))
      next
    }
    fp[[aai]] <- aaid
  } else {
    if(!grepl(".mix_gender.txt", aai) & !grepl(".bam", aai)) message(paste(" - unrecognized file format:", aai))
  } 
}
message(paste(" -", length(fp), "files read in total"))


### coerce
message("coerce into matrix")
fpcat <- do.call(rbind, fp)
dx <- unique(sort(fpcat[, 1]))
dy <- unique(sort(fpcat[, "id"]))

# filter read group data
if(!a.rg) {
  dx <- dx[!grepl(":", dx)]
  fpcat <- fpcat[is.element(fpcat[, 1], dx), ]
}

message(paste(" -", length(dx), "samples"))
message(paste(" -", length(dy), "positions"))
d <- matrix(NA, ncol=length(dx), nrow=length(dy))
colnames(d) <- dx
rownames(d) <- dy

if(length(fpcat)) {
  fplist <- split(fpcat[, c("id", "ratio")], fpcat[, "sample"])
  for(s in names(fplist)) d[fplist[[s]][, 1], s] <- fplist[[s]][, 2]
}


### cluster
message("calculate pairwise correlation")
d.cor <- suppressWarnings(cor(d, use="pair"))
d.ev <- crossprod(!is.na(d))
d.cor[d.ev < a.obs] <- NA

m <- data.frame(X=rep(dx, each=length(dx)), Y=dx, r=as.vector(d.cor))[as.vector(lower.tri(d.cor, diag=FALSE)), ]
m.cor <- m[which(m$r > a.cor), ]
message(" - ", paste(nrow(m.cor), "correlated samples"))

# form clusters
message("define clusters")
if(nrow(m.cor) > 0) {
  l <- union_find(m.cor[, c(1,2)])
  l <- l[order(sapply(l, length), decreasing=TRUE)]
  lv <- sapply(l, paste, collapse=" /// ")
  lv <- lv[sapply(l, length) > 1]
} else {
  l <- lv <- NULL
}
message(" - ", paste(length(lv), "clusters found"))


### output
message("write output")
dx.out <- lapply(dx, function(dxi) {
  #dxi <- "control_MB113"
  dxi.cor <- setdiff(unique(sort(unlist(m.cor[c(which(is.element(m.cor[, 1], dxi)), which(is.element(m.cor[, 2], dxi))), 1:2]))), dxi)
  dxi.cluster.i <- sapply(l, function(i) any(is.element(i, dxi)))
  if(any(dxi.cluster.i)) { dxi.cluster <- setdiff(l[[which(dxi.cluster.i)]], c(dxi, dxi.cor)) } else { dxi.cluster <- NULL }
  dxi.notcor <- names(which(d.cor[dxi, ] <= a.cor))
  dxi.notpo <- names(which(is.na(d.cor[dxi, ])))
  return(sapply(list(dxi, dxi.cor, dxi.cluster, dxi.notcor, dxi.notpo), paste, collapse=" /// "))
})

#write.table(do.call(rbind, dx.out), sep="\t", quote=FALSE, col.names=c("sample", "correlated", "cluster", "not_correlated", "no_observations"), row.names=FALSE, file=paste(a.out, "_sample.txt", sep=""))
write.table(lv, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste(a.out, "_cluster.txt", sep=""))
write.table(d.cor, sep="\t", quote=FALSE, file=paste(a.out, "_cor.txt", sep=""))
write.table(d.ev, sep="\t", quote=FALSE, file=paste(a.out, "_pairobs.txt", sep=""))

message(" - done")

