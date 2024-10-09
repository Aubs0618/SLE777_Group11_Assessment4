# SLE777_Group11_Assessment4

---
title: "Part 2- Gene and Protein Comparison "
output: html_document
date: "2024-10-08"
---

# Qestion 01 - 


# to download the packge before download the sequence file from ensemble to opetn relavant files, "seqinr" and "R.utils" package installation is needed
# to install "seqinr" package
```{r}
install.packages("seqinr")
```
# to install "R.utils" package
```{r}
install.packages("R.utils")
```


# to download E coli (GCA_003780985) to the direcotory,move to the ensemble,search the bacteria for FASTA file and  select the folder destination and pase the URL
# the file then unzipped

```{r,E coli  FASTA file download}
library("R.utils")
URL=("http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz")
download.file(URL,destfile="ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
```
# to check the downloaded file in the folder
```{r}
list.files()
```
# to download Campylobacter coli (GCA_003780985) to the direcotory,move to the ensemble,search the bacteria for FASTA file and  select the folder destination and pase the URL
# the file then unzipped
```{r,campylobacter FASTA file download}
URL=("https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_46_collection/campylobacter_coli_gca_003780985/cds/Campylobacter_coli_gca_003780985.PDT000395653.1.cds.all.fa.gz")
download.file(URL,destfile="Campylobacter_cds.fa.gz")
gunzip("Campylobacter_cds.fa.gz")
```

#to check the downloaded file
```{r}
list.files()
```

#Load the seqinr library to provides tools for biological sequence retrieval and analysis.
```{r}
library("seqinr")
```


#reads the uncompressed FASTA file ecoli_cds.fa containing E coli coding sequences into R using the read.fasta function from the seqinr package. 
```{r,Read Ecoli fasta }
ecoli <- seqinr::read.fasta("ecoli_cds.fa")
str(head(ecoli))

```

#reads the uncompressed FASTA file campyloi_cds.fa containing campylobacter coding sequences into R using the read.fasta function from the seqinr package. 
```{r,Read campylo fasta}
campylo <- seqinr::read.fasta("Campylobacter_cds.fa")
str(head(campylo))
```

# to find the number of coding sequences in each bacteria 
```{r,E-coli and Campylobacter gene coding sequence}
length(ecoli )
length(campylo )
```





#Question 02 -**



#To determine total coding DNA in E-coli and campylobacter
```{r}

```





#Question 03 -
# to calculate the length of all coding sequences in these two organisms.
```{r,total coding DNA ecoli }
length_ecoli <- as.numeric(summary(ecoli)[,1])
sum(length_ecoli)
```

```{r,total coding DNA Campylobacter}
length_campylo <- as.numeric(summary(campylo)[,1])
sum(length_campylo)
```


# to make a boxplot of coding sequence length in the organisms.
```{r}
boxplot(length_ecoli,length_campylo,xlab="Bacteria Species",ylab = "sequence length (bp)",main = "Comparison of Sequence length",names = c("Ecoli","Campylobacter"))
```
# to find lenth mean of Ecoli and Campylobacter
```{r}
mean(length_ecoli)
mean(length_campylo)
```
# to find lenth median of Ecoli and Campylobacter
```{r}
median(length_ecoli)
median(length_campylo)
```



#Question 04-

# to Calculate the frequency of DNA bases in the total coding sequences for organisms.
# to analyze the nucleotide composition of a  the DNA sequence, first need to unlist the cds, store it in the dna variable.
# then calculate the frequency of each nucleotide using dna_composition <- count(..., 1), which stores these frequencies in dna_composition.
# to count each nucleotide Srt() code is used.

# for ecoli
````{r}
dna_ecoli <- unlist(ecoli)
dna_composition_ecoli <- count(dna_ecoli ,1)
str(dna_composition_ecoli)
```
# for Campylobacter

```{r}
dna_campylo <- unlist(campylo)
dna_composition_campylo <- count(dna_campylo,1)
str(dna_composition_campylo)
```
# then visualize this data with barplot where the x-axis represents the nucleotides (A, C, G, T), the y-axis shows their frequencies, and the plot title indicates that it displays the nucleotide composition of coding sequences.


```{r}
dna_composition_ecoli <- count(dna_ecoli ,1)
barplot(dna_composition_ecoli ,xlab="Nucleotides",ylab="Frequency", main="E coli CDS composition")

dna_composition_campylo <- count(dna_campylo ,1)
barplot(dna_composition_campylo ,xlab="Nucleotides",ylab="Frequency", main="Campylobacter CDS composition")

```
#Question 04 -Part 02
#Calculate the frequency of aminoacid in the total protein sequence. 
#to translate the DNA sequences for protein sequence use a for loop  to translate all the cds sequences in the genomes
#lapply applies the translate function to each element of the cds list, resulting in prot, a list of translated protein sequences.


# for ecoli
```{r,Ecoli Protein sequence}
prot_ecoli <- lapply(ecoli, translate)
```


```{r}
unlist(prot_ecoli)
protein_sequence <- unlist(ecoli)
length(protein_sequence)
aa_ecoli <- unique(prot_ecoli)
aa_ecoli <- aa_ecoli[aa_ecoli != "*"]
 <-count(protein_sequence,wordsize=1,alphabet=aa_ecoli)

```



#for campylobacter
```{r,campylobacter Protein sequence}
prot_campylo <- lapply(campylo, translate)
str(prot_campylo)
```






