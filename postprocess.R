#Expecting the following files (outputs from the previous step): 
#combined_aligned.fasta - an aligned fasta file with the two mitogenomes titled *.1 and *.2.fa for each sample as well as the reference sample. Each sequence should only span one line
#For each phased sample, a vcf file prefixed with sample name showing the quality/phase of the various genotypes (suffixed *.phased_SNPs.vcf)

postprocess <- function(working_dir) {#1A
#e.g. working_dir <- "/Users/alanaalexander/Dropbox/polg_mice"

samplenames <- unique(gsub(".phased_SNPs.vcf","",list.files(working_dir, pattern=".phased_SNPs.vcf")))

input_fasta <- readLines(paste(working_dir,"/","combined_aligned.fasta",sep=""))
nochars <- length(strsplit(input_fasta[2],"")[[1]])

outputmat <- matrix(NA,nrow=(nochars+1),ncol=(length(seq(1,length(input_fasta),2))+length(samplenames)+1))
outputmat[1,1:(length(seq(1,length(input_fasta),2)))] <- input_fasta[(seq(1,length(input_fasta),2))]
  
for (i in 1:(length(seq(1,length(input_fasta),2)))) {
  outputmat[2:(dim(outputmat)[1]),(i+1)] <- strsplit(input_fasta[i*2],"")[[1]]
}
  
  
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
