#Expecting the sample name as a prefix for the following files: two fasta files per sample (suffixed *.1.fa and *.2. fa), and a vcf file with the quality/phase of the various genotypes (suffixed *.phased_SNPs.vcf), and nothing else in the working directory

postprocess <- function(working_dir) {#1A
#e.g. working_dir <- "/Users/alanaalexander/Dropbox/polg_mice"

samplenames <- unique(gsub("\\..*","",list.files(working_dir)))

for (i in samplenames) {
  M1_name <- paste(i,".1.fa",sep="")
  M2_name <- paste(i,".2.fa",sep="")
  VCF_name <- paste(i,".phased_SNPs.vcf",sep="")
  M1 <- strsplit(readLines(paste(working_dir,"/",M1_name,sep=""))[2],"")[[1]]
  output <- matrix(NA,ncol=3,nrow=length(M1))
  output[,1] <- M1
  output[,2] <- (strsplit(readLines(paste(working_dir,"/",M2_name,sep=""))[2],"")[[1]])


} #1B

 GT:AD:DP:GQ:HP:PL       
GT 0/1:
AD Ref - 949 Alt - 3247
DP 4196
GQ 99
HP 158-1,158-2
PL 100237,0,20504


GT 0/1:
AD Ref 363 Alt 1943
DP 2306
GQ 99
HP 158-2,158-1
62838,0,5211:35.01

Genotype(GT):
Allelic depths for the ref and alt alleles in the order listed(AD):
Approximate read depth (reads with MQ=255 or with bad mates are filtered)(DP):
Genotype Quality(GQ):
Read-backed phasing haplotype identifiers(HP):
Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification(PL):
