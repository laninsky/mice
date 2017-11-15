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

#creating an inputmatrix which will have a column for each sequence and the relative ref_pos
inputmat <- matrix(NA,nrow=(nochars+1),ncol=(length(seq(1,length(input_fasta),2))+1))
inputmat[1,1:(length(seq(1,length(input_fasta),2)))] <- input_fasta[(seq(1,length(input_fasta),2))]

#populating the rows with sequence from the fasta file  
for (i in 1:(length(seq(1,length(input_fasta),2)))) {
  inputmat[2:(dim(inputmat)[1]),i] <- strsplit(input_fasta[i*2],"")[[1]]
}

#finding out what the reference column is (not present in the vcf sample names)
ref_column <- which(!(gsub("\\..*","",gsub(">","",inputmat[1,]))[1:(length(seq(1,length(input_fasta),2)))] %in% samplenames))

#creating a column giving the position relative to the reference sequence  
inputmat[1,((length(seq(1,length(input_fasta),2)))+1)] <- "ref_pos"

#populating this column - giving indel positions "compound" sites e.g. 484, 484_1, 484_2 so they are still relative to the reference  
ref_site <- 1
for (i in 2:(dim(inputmat)[1])) {
  if(!(inputmat[i,ref_column]=="-")) {
    inputmat[i,((length(seq(1,length(input_fasta),2)))+1)] <- ref_site
    ref_site <- ref_site + 1
    site_suffix <- 1
  } else {
    inputmat[i,((length(seq(1,length(input_fasta),2)))+1)] <- paste((ref_site-1),"_",site_suffix,sep="")
    site_suffix <- site_suffix + 1    
  }
}

#Getting the output matrix ready. This is going to have the ref_pos in the first column, ref base in the next
outputmat <- inputmatrix[,c(((length(seq(1,length(input_fasta),2)))+1),ref_column)]
  
#pulling in the vcf to work out if sites can be confidently phased or not  
for (i in samplenames) { #2A
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
  
  #Pulling out the sequence for the sample we are interested in, so we can get the phasing blocks defined
  samplecols <- c(which(inputmatrix[1,] %in% paste(">",i,".1",sep="")),which(inputmatrix[1,] %in% paste(">",i,".2",sep="")))
  tempsamplematrix <- matrix(NA,ncol=6,nrow=dim(inputmatrix)[1])
  tempsamplematrix[,1] <- outputmat[,1]
  tempsamplematrix[,2:3] <- inputmatrix[,samplecols]
  tempsamplematrix[1,4] <- paste(i,"_phasing",sep="")
  
  #Identifying the sites where there is more than one allele
  for (j in 2:dim(inputmatrix)[1]) {
    if(!(tempsamplematrix[j,2]==tempsamplematrix[j,3])) {
      tempsamplematrix[j,4] <- "Unknown"
    }
  }
  #Propogating the sequence over for sites that do not show signs of alternate alleles
  tempsamplematrix[(which(is.na(tempsamplematrix[,4]))),5] <- tempsamplematrix[(which(is.na(tempsamplematrix[,4]))),2]
  tempsamplematrix[(which(is.na(tempsamplematrix[,4]))),6] <- tempsamplematrix[(which(is.na(tempsamplematrix[,4]))),2]
  
  #Looping through the VCF file and propogating those bases to the output seq if the genotype is certain
  for (j in 2:dim(VCFmat)[1]) {#3A
    if(nchar(VCFmat[j,4])==nchar(gsub("PQ","",VCFmat[j,4]))) {#4A
       if(nchar(VCFmat[j,2])==nchar(VCFmat[j,3])) {#5A
       tempsamplematrix[(which(tempsamplematrix[,1]==VCFmat[j,1])),5] <- VCFmat[j,2]
       tempsamplematrix[(which(tempsamplematrix[,1]==VCFmat[j,1])),6] <- VCFmat[j,3]
       tempsamplematrix[(which(tempsamplematrix[,1]==VCFmat[j,1])),4] <- paste("Unknown_",gsub(",",":",unlist(strsplit(VCFmat[j,5],":"))[2]),sep="")
       } else { #5AB
       # what to do when no PQ and it is an indel  
       }  #5B
    } else { #4AB
       # what to do when the phase is confident - need to do this for both substitutions and indels
    }#4B  
  }#3B
}#2B  
  
  
  #UP TO HERE ######
  ### The Vcf file data then needs to get propogated into the tempsamplematrix
  ### If the vcf file doesn't have a PQ, then that site has the readcount suffixed to "Unknown" (e.g. "Unknown"_readcount1:readcount2) and cols 2/3 get copied to cols 5/6 (but phase is unknown)
  ### If the first site has haplotypes, they get poropgated
  ### If it does have a PQ, the order of the haplotypes needs to get saved in a variable (so it can be checked later on)
  ### The respective bases need to get put into the appropriate columns in 5 and 6
  ### In the samplename_phasing column I think the output should be formatted the following
  ### Hapname1_readcount:Hapname2_readcount
  ### Column 5 always corresponds to Hapname1 and Column 6 always corresponds to Hapname6
  
  ### Need to deal with underscore rows
  ### After that is done, need to find all places with "Unknown" - both counting upwards and downwards from that point, the "Unknown" spreads
  ### After this, can push through the matrix copying the phased state above
  
  
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
