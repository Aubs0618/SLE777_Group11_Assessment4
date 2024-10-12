# SLE777_Group11_Assessment4
## Part 1: Gene Expression and Growth Data Analysis

# Step 1: Download and read in the gene_expression.tsv file

```{r,}
download.file("https://github.com/ghazkha/Assessment4/raw/main/gene_expression.tsv", destfile = "gene_expression.tsv")
gene_data <- read.table("gene_expression.tsv", header = TRUE, sep = "\t", row.names = 1)
```

#Show the first six genes
```{r,top 6}
head(gene_data, 6)
```

# Step 2: Add a new column for the mean expression of each gene
```{r,mean}
gene_data$Mean_Expression <- rowMeans(gene_data)
```

#Show the first six genes with the mean expression
```{r,}
head(gene_data, 6)
```

# Step 3: Top 10 Genes with Highest Mean Expression

#List the top 10 genes with the highest mean expression
```{r,}
top_genes <- gene_data[order(-gene_data$Mean_Expression), ]
top10_genes <- head(top_genes, 10)
top10_genes
```

# Step 4: number of genes with a mean <10
```{r}
meanless10<- subset(gene_data,`Mean_Expression` < 10)
nrow(meanless10)
```

# step 5:histogram plot of the mean values
#Create a histogram of the mean expression values

```{r}
hist(top10_genes$Mean_Expression,xlab="Top 10 genes",main =" Top gene 10 mean expression")
```

# Step 6: Import and Analyze Growth Data
#Download and read in the growth_data.csv file
```{r}
download.file("https://github.com/ghazkha/Assessment4/raw/main/growth_data.csv", destfile = "growth_data.csv")
growth_data <- read.csv("growth_data.csv")
```

#Show the column names
```{r}
colnames(growth_data)
```

# Step 7: subset southwest site
#Reshape the data from wide to long format
```{r}
southwest <- subset(growth_data,Site == "southwest", c("Circumf_2005_cm", "Circumf_2020_cm"))
head(southwest)
```

# step 8:southwest mean over the years
```{r}
colMeans(southwest)
```

# step 9 : circumferenc standard deviation each year
```{r}
sd(southwest[,1])
```
```{r}
sd(southwest[,2])
```

# step 10:  subset northeast data site

```{r}
northeast <- subset(growth_data,Site == "northeast", c("Circumf_2005_cm", "Circumf_2020_cm"))
head(northeast)
```

# step 11: northeast mean of the years
```{r}
colMeans(northeast)
```

# step 12: circumferenc standard deviation  northeeast of the years
```{r}
sd(northeast[,1])
```

```{r}
sd(northeast[,2])
```

# Step 13 Box Plot of Growth Data
#Create a boxplot of tree circumference by year
```{r}
boxplot(northeast [,1],northeast [,2],southwest[,1],southwest[,2],xlab="number of cylinders",ylab = "circumfrence",main="Tree circumference",names = c("northeast 2005","northeast 2020","southwest 2005","southwest 2020"))
```
# step 14 - subset data 2020 and 2010 to take means value of diffrences over the 10 year
#find the mean values of the different of the growth between two years
#for northeast
```{r}
northeast_grwdif<- subset(growth_data,Site == "northeast", c("Circumf_2010_cm","Circumf_2020_cm"))
northeast_grwdif$growth_diffe <-northeast_grwdif$Circumf_2020_cm -northeast_grwdif$Circumf_2010_cm
colMeans(northeast_grwdif)
```
#for soutwest
```{r}
southwest_grwdif <- subset(growth_data,Site == "southwest", c("Circumf_2010_cm","Circumf_2020_cm"))
southwest_grwdif$growth_diffe <-southwest_grwdif$Circumf_2020_cm -southwest_grwdif$Circumf_2010_cm
colMeans(southwest_grwdif)
```




title: "Part 2- Gene and Protein Comparison of E-coli and Campylobacter coli "
output: html_document
date: "2024-10-08"
---



# Qestion 01 - 


#before downloading the sequence file from ensemble to opetn relavant files, it is needed to downnload "seqinr" and "R.utils" package installation.
#to install "seqinr" package 

```{r, installation package seq}
install.packages("seqinr")
```
#to install "R.utils" package

```{r installation package R utilities}
install.packages("R.utils")
```

#to download E coli (GCA_003780985) to the direcotory, from the ensemble, FASTA file was searched and selected the folder destination and pase the URL.
#the file then unzipped

```{r,E coli download FASTA}
library("R.utils")
URL=("http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz")
download.file(URL,destfile="ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
```
#to check the downloaded file whether the file is downloaded properly.
```{r}
list.files()
```
#to download Campylobacter coli (GCA_003780985) to the direcotory, from the ensemble, FASTA file was searched and selected the folder destination and pase the URL.
#the file then unzipped
```{r,campylobacter FASTA file download}
URL=("https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_46_collection/campylobacter_coli_gca_003780985/cds/Campylobacter_coli_gca_003780985.PDT000395653.1.cds.all.fa.gz")
download.file(URL,destfile="Campylobacter_cds.fa.gz")
gunzip("Campylobacter_cds.fa.gz")
```

#to check the downloaded file whether the file is downloaded properly.
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

#to find the number of coding sequences in each bacteria 
```{r,E-coli and Campylobacter gene coding sequence}
length(ecoli )
length(campylo )
```

#to compare the result,data frame is constructed
#variable is made for the species names as "Bacterial_strains"
#variable is made for the number of coding sequence as "number of coding sequence"
#make data frame by using data.frame code.
#run str() to check that the resulting dataframe has the intended structure. 
#include stringsAsFactors=FALSE to protect character strings being converted to factors.

```{r}
Bacterial_strains<- c("E coli","Campylobacter Coli")
Number_of_coding_sequence <- c(4239,1976)
species_codseq <- data.frame(Bacterial_strains,Number_of_coding_sequence,stringsAsFactors = FALSE)
str(species_codseq)
species_codseq
```



# Question 02 -

#To determine total coding DNA in E-coli and campylobacter ,
#unlisted the csv and saved as a variable tot_cod_dna_ecoli 
#using"count" command , the DNA nucleotides are counted
#To determine total coding DNA in E-coli 
```{r,total coding DNA ecoli}
tot_cod_dna_ecoli <- unlist(ecoli)
count(tot_cod_dna_ecoli,1)
```
#To determine total coding DNA in campylobacter 
```{r,total coding DNA campylo}
tot_cod_dna_campylo <- unlist(campylo)
count(tot_cod_dna_campylo,1)
```

#to present the result,data frame is constructed
#variables are made for neucleotides for each species named as "Number_of_coding_DNA_Ecoli " and "Number_of_coding_DNA_Campylobacter" 
#variable is made for the number of Nucleotides as "Nucleiotieds"
#make data frame by using data.frame code.
#included the "stringsAsFactors=FALSE" to protect character strings being converted to factors.

```{r,dataframe codDNA}
Number_of_coding_DNA_Ecoli <- c(955768,977594,1088501,956665 )
Number_of_coding_DNA_Campylobacter <- c(623683,229717,318984,554434 )
Nucleiotieds <-c("A","C","G","T")
cod_DNA<- data.frame(Nucleiotieds,Number_of_coding_DNA_eco,Number_of_coding_DNA_camp,stringsAsFactors = FALSE)
cod_DNA
```

#to calculte the total colsum code is used.
```{r}
colSums(cod_DNA[,2:3])
```
#to construct a dataframe
#variablea are made for each species names as "Bacterial_strains " 
#variable is made for the coding DNA "Bacterial_strains"
#make data frame by using data.frame code.
#include stringsAsFactors=FALSE to protect character strings being converted to factors.
```{r}
Bacterial_strains<- c("E coli","Campylobacter Coli")
Toal_Coding_DNA <-c(3978528,1726818 )
df_tot_DNA <- data.frame(Bacterial_strains,Toal_Coding_DNA,stringsAsFactors = FALSE)
df_tot_DNA 
```


# Question 03 -

#to calculate the length of all coding sequences in these two organisms,
#variable is made the strings to convert numetric values  and sum is used to calculate the total lengths.

```{r,total coding DNA ecoli }
length_ecoli <- as.numeric(summary(ecoli)[,1])
sum(length_ecoli)
```
```{r,total coding DNA Campylobacter}
length_campylo <- as.numeric(summary(campylo)[,1])
sum(length_campylo)
```

#to make a boxplot of coding sequence length in the organisms.
```{r,Boxplot seqences}
boxplot(length_ecoli,length_campylo,xlab="Bacteria Species",ylab = "sequence length (bp)",main = "Comparison of Sequence length",names = c("Ecoli","Campylobacter"))
```

#to find lenth mean of Ecoli and Campylobacter mean funtion is used.
```{r}
mean(length_ecoli)
mean(length_campylo)
```

#to find lenth mean of Ecoli and Campylobacter median funtion is used. 
```{r}
median(length_ecoli)
median(length_campylo)
```




# Question 04-

#part 01 - for DNA

#to Calculate the frequency of DNA bases in the total coding sequences for organisms.
#to analyze the nucleotide composition of a  the DNA sequence, first need to unlist the cds, store it in the dna variable.
#then calculate the frequency of each nucleotide using dna_composition <- count(..., 1), which stores these frequencies in dna_composition.
#to count each nucleotide Srt() code is used.

#for ecoli
````{r,ecoi frequency}
dna_ecoli <- unlist(ecoli)
dna_composition_ecoli <- count(dna_ecoli ,1)
str(dna_composition_ecoli)
```
# for Campylobacter
```{r,campylo frequency}
dna_campylo <- unlist(campylo)
dna_composition_campylo <- count(dna_campylo,1)
str(dna_composition_campylo)
```

# then visualize this data with barplot where the x-axis represents the nucleotides (A, C, G, T), the y-axis shows their frequencies, and the plot title indicates that it displays the nucleotide composition of coding sequences.
```{r,barplot for species}
dna_composition_ecoli <- count(dna_ecoli ,1)
barplot(dna_composition_ecoli ,xlab="Nucleotides",ylab="Frequency", main="E coli CDS composition")

dna_composition_campylo <- count(dna_campylo ,1)
barplot(dna_composition_campylo ,xlab="Nucleotides",ylab="Frequency", main="Campylobacter CDS composition")
```

# -Part 02 - for protein

#Calculate the frequency of amino acid in the total protein sequence. 
#to translate the DNA sequences for protein sequence use a for loop  to translate all the cds sequences in the genomes
#lapply applies the translate function to each element of the cds list, resulting in prot, a list of translated protein sequences.
# to convert the list to string ,unlits the protein sequence using "unlist" command

# for ecoli ,to detect total  each amino acid frequencey in protein sequence 
```{r,Ecoli Protein sequence}
prot_ecoli <- lapply(ecoli, translate)
aa <- unique(prot_ecoli[[2]])
aa <- aa[aa != "*"]
Proseq_ecoli <- unlist(prot_ecoli)
count(Proseq_ecoli,wordsize=1, alphabet =aa)
```

# for campylobacter ,# for ecoli ,to detect total  each amino acid frequencey in protein sequence 
```{r,campylobacter Protein sequence}
prot_campylo <- lapply(campylo, translate)
aa <- unique(prot_campylo[[2]])
aa <- aa[aa != "*"]
Proseq_campylo <- unlist(prot_campylo)
count(Proseq_campylo,wordsize=1, alphabet =aa)
```

# to determine the frequency of each amino acid in all genes,this counts the occurrences of each amino acid in the sequence of the prot list, using aa as the set of possible amino acids. 
# unlist the sequence using "unlist" and count the sequence and alphabet used to reperesent the each amonoacid,word size give to repersent as one alphabet

```{r,protein composition charts}
prot_composition_ecoli <- count(Proseq_campylo,wordsize=1, alphabet =aa)
barplot(prot_composition_ecoli ,xlab="Amino acids",ylab="Frequency", main="E coli aminoacid composition")

prot_composition_campylo <- count(Proseq_ecoli,wordsize=1, alphabet =aa)
barplot(prot_composition_campylo ,xlab="Amino acids",ylab="Frequency", main="Campylobacter amino acid composition")
```


# Question 05 -

# to Create a codon usage table and quantify the codon usage bias among all coding sequences.
# computes the Relative Synonymous Codon Usage (RSCU) values, which is a measure as to how preferred a codon of a particular amino acid is compared to the other codons which also encode it,
#for that use "rscu" fucntion
# unlist teh protein sequecne and then define newvariales relatively to calcutlte the sequences.
#to get these values in a data frame ,use the option as.data.frame=True
#uco command is used to determine the calculates some codon usage information.

# for Ecoli ,
```{r,codon usage ecoli}
prosec_ecoli <- unlist(prot_ecoli)
uco(prosec_ecoli,index="rscu",as.data.frame=TRUE)

```
#for campylobacter
```{r,codon usage campylobacter}
prosec_campylo <- unlist(prot_campylo)
uco(prosec_campylo,index="rscu",as.data.frame=TRUE)

```

# providing evidece as chart for frequencies of genes in each species
# define the output of variable of codons for create a dataframe

```{r, ecoli RASCU chart}
dataframe_ecoli <- uco(prosec_ecoli,index="rscu",as.data.frame=TRUE)
hist(dataframe_ecoli$RSCU,x=,xlab="RSCU Value", main="Relative synonymous codon usage  for Ecoli codons")
```
```{r,campylobacter RSCU chart}
dataframe_campylo <- uco(prosec_campylo,index="rscu",as.data.frame=TRUE)
hist(dataframe_campylo$RSCU,x=,xlab="RSCU Value", main="Relative synonymous codon usage for Campylobacter codons")
```

# Question 06-


#to identify particular sequences that are over- or under-represented in an organism,K-mer profiling is done for each species.

# unlist the sequence using "unlist" and define it as new variable "ecoli_kcounts" 
# for output acmino acid, nmber of alphabets with 3 leter is given for represent the codon.


# for ecoli
```{r,kmer Ecoli}
prots_ecoli <- unlist(prot_ecoli)
ecoli_kcounts <- count(prots_ecoli ,wordsize=3,alphabet=aa)
head(ecoli_kcounts)
```


# for campylobacter
```{r,Kmer campylobacter}
prots_campylo <- unlist(prot_campylo)
campylo_kcounts <- count(prots_campylo ,wordsize=3,alphabet=aa)
head(campylo_kcounts)
```



```{r}
ecofreq <- count(prots_ecoli,wordsize=3,alphabet=aa,freq=TRUE)
datafre_eco <- as.data.frame(ecofreq)
sorted_eco<-datafre_eco[order(datafre_eco $`Freq`),]
sorted_eco[c(1:10),c(1,ncol(sorted_eco))]
```


# to compare the sequnec,barplots are created using "barplot" and x = protein sequence and y=frequencies
```{r}

ecolifreq <- count(prots_ecoli,wordsize=1,alphabet=aa,freq=TRUE)
barplot(ecolifreq ,xlab="Proteins sequence",ylab="Frequency", main="ecoli  protein composition")


campylofreq <- count(prots_campylo,wordsize=1,alphabet=aa,freq=TRUE)
barplot(campylofreq ,xlab="Proteins sequence",ylab="Frequency", main="Campylobacter  protein composition")

```

