#Expecting the following files (outputs from the previous step): 
#combined_aligned.fasta - an aligned fasta file with the two mitogenomes titled *.1 and *.2.fa for each sample as well as the reference sample. Each sequence should only span one line
#For each phased sample, a vcf file prefixed with sample name showing the quality/phase of the various genotypes (suffixed *.phased_SNPs.vcf)

postprocess <- function(working_dir) {#1A
#e.g. working_dir <- "/Users/alanaalexander/Dropbox/polg_mice"

#Getting the sample names from the vcf files in the folder
samplenames <- unique(gsub(".phased_SNPs.vcf","",list.files(working_dir, pattern=".phased_SNPs.vcf")))

#Grabbing the aligned fasta files and figuring out the total number of nucleotides in the alignment
input_fasta <- readLines(paste(working_dir,"/","combined_aligned.fasta",sep=""))
nochars <- length(strsplit(input_fasta[2],"")[[1]])

#creating an outputmatrix which will have a column for each sequence, the relative ref_pos, and whether the sites have been solidly phased per sample 
outputmat <- matrix(NA,nrow=(nochars+1),ncol=(length(seq(1,length(input_fasta),2))+length(samplenames)+1))
outputmat[1,1:(length(seq(1,length(input_fasta),2)))] <- input_fasta[(seq(1,length(input_fasta),2))]

#populating the rows with sequence from the fasta file  
for (i in 1:(length(seq(1,length(input_fasta),2)))) {
  outputmat[2:(dim(outputmat)[1]),i] <- strsplit(input_fasta[i*2],"")[[1]]
}

#finding out what the reference column is (not present in the vcf sample names)
ref_column <- which(!(gsub("\\..*","",gsub(">","",outputmat[1,]))[1:(length(seq(1,length(input_fasta),2)))] %in% samplenames))

#creating a column giving the position relative to the reference sequence  
outputmat[1,((length(seq(1,length(input_fasta),2)))+1)] <- "ref_pos"

#populating this column - giving indel positions "compound" sites e.g. 484, 484_1, 484_2 so they are still relative to the reference  
ref_site <- 1
for (i in 2:(dim(outputmat)[1])) {
  if(!(outputmat[i,ref_column]=="-")) {
    outputmat[i,((length(seq(1,length(input_fasta),2)))+1)] <- ref_site
    ref_site <- ref_site + 1
    site_suffix <- 1
  } else {
    outputmat[i,((length(seq(1,length(input_fasta),2)))+1)] <- paste((ref_site-1),"_",site_suffix,sep="")
    site_suffix <- site_suffix + 1    
  }
}
 
#pulling in the vcf to work out if sites can be confidently phased or not  
for (i in samplenames) {
  #Getting a matrix together for the VCF file per sample
  VCF_name <- paste(i,".phased_SNPs.vcf",sep="")
  tempVCF <- readLines(paste(working_dir,"/",VCF_name,sep=""))
  tempVCF <- tempVCF[grep("^\\#",tempVCF,invert=TRUE)]
  VCFmat <- matrix(NA,ncol=5,nrow=length(gsub("\\..*GT","",gsub("^.*?\t","",tempVCF))))
  VCFmat[,1] <- unlist(strsplit(tempVCF,"\t"))[seq(2,length(unlist(strsplit(tempVCF,"\t"))),10)]
  VCFmat[,2] <- unlist(strsplit(tempVCF,"\t"))[seq(4,length(unlist(strsplit(tempVCF,"\t"))),10)]
  VCFmat[,3] <- unlist(strsplit(tempVCF,"\t"))[seq(5,length(unlist(strsplit(tempVCF,"\t"))),10)]
  VCFmat[,4] <- unlist(strsplit(tempVCF,"\t"))[seq(9,length(unlist(strsplit(tempVCF,"\t"))),10)]
  VCFmat[,5] <- unlist(strsplit(tempVCF,"\t"))[seq(10,length(unlist(strsplit(tempVCF,"\t"))),10)]
  
  
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
