# Fingerprint toolset

The fingerprint toolset provides functionality for in-silico genotyping. The toolset can be used to __confirm expected genotype matches__ (e.g. tumor & normal, WGBS & 450k array), or to __detect unexpected matches__ (e.g. different DNA sources, sample swaps).

This code is related to Volker Hovestadt's [methylCtools](https://github.com/hovestadt/methylCtools).

## Usage

### Software requirements

This is legacy code. The software versions that probably work with this code are these

 * Python 2.7.9
   * pysam 0.7.7
 * R 3.4.0
   * minfi 1.24.0
   * getopt 1.20.0

### Generate fingerprint files (.fp)

 - From BAM alignment files (python script):
   ```Shell
   python fingerprint/bsnp.py \
      fingerprint/snp141Common_RS-CG.n150.vh20151103.bed \
      sample.bam \
      > sample.bam.fp
   ```
 - From 450k data files (R script, RData input file must contain a [minfi](https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html) RGset object of a single sample):
   ```Shell
   Rscript fingerprint/asnp.R \
      fingerprint/snp141Common_RS-CG.n150.vh20151103.bed \
      sample.RData \
      sample sample.fp
   ```

### Analyze fingerprint files to identify genotype matches

 - R script, input arguments can be individual files and directories, separated by comma:
   ```Shell
   Rscript fingerprint/bsnp_analyze.R \
      -i sampleset/ \
      -o sampleset
   ```

More detailed information is provided with the help messages and headers of the tools. It is recommended to __understand__ and __adapt__ the settings used in `bsnp_analyze.R`. often it also helps to look at the pairwise correlation table directly (e.g. in Excel using conditional color formatting).

The `bsnp_analyze.R` tool follows a very simple statistical approach (i.e. just calculating the pearson correlation), but provides reasonable results. Because of that it requires SNPs that are highly variable within the population (avHet >0.2). Anyone is invited to come up with a more sophisticated approach.
 
## Included SNP files

- `snp138Common.n1000.vh20140318.bed`: A set of 1000 SNPs (avHet >0.25), selected for coverage in RNA-, Exome- and ChIP-seq datasets, including 450k RS probes. **This file is recommended for comparing sequencing data**.
- `snp141Common_RS-CG.n436.vh20151103.bed`: The original set of 436 SNPs (major allele freq <0.95) that can be assessed on the 450k arrays (RS probes and informative CG type I probes).
- `snp141Common_RS-CG.n385.vh20151104.bed`: A subset of the previous (consistency with WGS data, i.e. some probes with very high/low intensities were removed).
- `snp141Common_RS-CG.n192.vh20151103.bed`: A subset of the previous (avHet >0.2). **This file is recommended for comparing 450k data, and for comparing 450k to sequencing data**.
- `snp141Common_RS-CG.n50.vh20151104.bed`: A subset of the previous (filtered for coverage in Exome- and RNA-seq data).

## Author

* Volker Hovestadt

This code is related to [methylCtools](https://github.com/hovestadt/methylCtools).

## Changelog

- 1.1 [4. November 2015]
  - New Git repository for old code.
  - Volker Hovestadt, German Cancer Research Center, 2015, v1.1
