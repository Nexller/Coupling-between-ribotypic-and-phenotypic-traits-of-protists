############################################################################################################################################################################################################
############################################################################################################################################################################################################
# Colpoda steinii and Colpoda inflata


# Removing primers using cutadapt
for sample in $(cat samples)
do 
	cutadapt -g CCAGCASCYGCGGTAATTCC -a TYRATCAAGAACGAAAGT -G ACTTTCGTTCTTGATYRA -A GGAATTACCGCRGSTGCTGG -m 180 
	-o ${sample}_R1_trimmed.fq.gz -p ${sample}_R2_trimmed.fq.gz ${sample}.raw_1.fq.gz ${sample}.raw_2.fq.gz &
done


# Processing with DADA2 in R
#!/usr/bin/Rscript
# Setting up R working environment
setwd("/database2/Zou_Colpoda_polymorphisms/DADA2_V2/soil")
library(dada2)
library(DECIPHER)
load("SILVA_SSU_r138_2019.RData")

samples <- scan("samples", what="character")

for (i in samples){
  
  print (paste0("Processing sample: ", i))
  
  i.forward_reads <- paste0(i, "_R1_trimmed.fq.gz")
  i.reverse_reads <- paste0(i, "_R2_trimmed.fq.gz")
  
  
  i.filtered_forward_reads <- paste0(i, "_R1_filtered.fq.gz")
  i.filtered_reverse_reads <- paste0(i, "_R2_filtered.fq.gz")
  

  # Quality trimming/filtering
  i.filtered_out <- filterAndTrim(i.forward_reads, i.filtered_forward_reads,
                                  i.reverse_reads, i.filtered_reverse_reads, 
                                  minLen=180,
                                  # truncLen=c(290,200),
                                  # trimLeft = c(5,10),
                                  rm.phix = TRUE,
                                  maxN = 0,
                                  truncQ = 2,
                                  maxEE=c(2,2),
                                  compress = T, 
                                  verbose = TRUE, 
                                  multithread = TRUE)
  
  
  # Generating an error model
  i.err_forward_reads <- learnErrors(i.filtered_forward_reads, multithread=TRUE)
  i.err_reverse_reads <- learnErrors(i.filtered_reverse_reads, multithread=TRUE)
  
  # Dereplication
  i.derep_forward <- derepFastq(i.filtered_forward_reads, verbose=TRUE)
  names(i.derep_forward) <- i 
  
  i.derep_reverse <- derepFastq(i.filtered_reverse_reads, verbose=TRUE)
  names(i.derep_reverse) <- i
  
  
  # Inferring ASVs using core command dada
  i.dada_forward <- dada(i.derep_forward, err=i.err_forward_reads, pool="pseudo", multithread=TRUE)
  i.dada_reverse <- dada(i.derep_reverse, err=i.err_reverse_reads, pool="pseudo", multithread=TRUE)
  
  
  # Merging forward and reverse reads
  i.merged_amplicons <- mergePairs(i.dada_forward, i.derep_forward, 
                                   i.dada_reverse, i.derep_reverse, 
                                   # trimOverhang=TRUE, 
                                   # minOverlap=170,
                                   verbose = TRUE)
  

  # Generating a count table
  i.seqtab <- makeSequenceTable(i.merged_amplicons)

  
  # Chimera identification
  i.seqtab.nochim <- removeBimeraDenovo(i.seqtab, verbose=T) 
  

  # set a little function
  getN <- function(x) sum(getUniques(x))


  # making a data table for paired-end reads
  i.summary_tab <- data.frame(row.names=i, dada2_input=i.filtered_out[,1],
                          filtered=i.filtered_out[,2], dada_f=sapply(i.dada_forward, getN),
                          dada_r=sapply(i.dada_reverse, getN), merged=sapply(i.merged_amplicons, getN),
                          nonchim=rowSums(i.seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(i.seqtab.nochim)/i.filtered_out[,1]*100, 1))

  # making a data table for single reads
  # i.track <- cbind(i, getN(i.filtered_out), getN(i.dada_forward), getN(i.merged_amplicons), rowSums(i.seqtab), rowSums(i.seqtab.nochim))


  # save table into a txt file
  write.table(i.summary_tab, file = paste0(i, "_data_summary.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)
  

  # Extracting the standard goods from DADA2
  i.asv_seqs <- colnames(i.seqtab.nochim)
  i.asv_headers <- vector(dim(i.seqtab.nochim)[2], mode="character")
  
  
  for (j in 1:dim(i.seqtab.nochim)[2]) {
    i.asv_headers[j] <- paste(paste0(">", i), j, sep="_")
  }
  

  # save ASVs sequences
  i.asv_fasta <- c(rbind(i.asv_headers, i.asv_seqs))
  head(i.asv_fasta)
  write(i.asv_fasta, paste0(i, ".fasta"))
  

  # save ASVs count table
  i.asv_tab <- t(i.seqtab.nochim)
  head(i.asv_tab)

  row.names(i.asv_tab) <- sub(">", "", i.asv_headers)
  write.table(i.asv_tab, paste0(i, "_counts.tsv"), sep="\t", quote=F, col.names=NA)


  # Assigning taxonomy using IdTaxa command in DECIPHER packages 
  i.dna <- DNAStringSet(getSequences(i.seqtab.nochim))
  i.tax_info <- IdTaxa(test=i.dna, trainingSet=trainingSet, strand="both", processors=NULL)
  
  ranks <- c("rootrank", "kingdom", "division" , "phylum", "class", "order", "family", "genus", "species")
  i.asv_tax <- t(sapply(i.tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(i.asv_tax) <- ranks
  
  rownames(i.asv_tax) <- gsub(pattern=">", replacement="", x=i.asv_headers)
  
  # Save taxonomy results
  write.table(i.asv_tax, paste0(i, "_Arch_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)

  print (paste0(i, " of Colpoda have done!"))
  print ("--------------------------------------------------------------------------")

}



############################################################################################################################################################################################################
############################################################################################################################################################################################################
# Euplotes vannus and Strombidium sulcatum

# Removing primers using cutadapt
for sample in $(cat samples)
do
	echo "正在处理的样品: $sample"
	cutadapt -a ^GAADCTGYGAAYGGCTC...GGAGGGCAAGTCTGGT \
	-m 50 -o ${sample}_trimmed.fq.gz ${sample}.fastq >> cutadapt_primer_trimming_stats.txt 2>&1
done


# Processing with DADA2 in R
#!/usr/bin/Rscript
# Setting up R working environment
setwd("/database2/Zou_Colpoda_polymorphisms/DADA2_V2/marine")

library(dada2)
library(DECIPHER)
load("SILVA_SSU_r138_2019.RData")

samples <- scan("samples", what="character")


for (i in samples){
  
  print (paste0("Processing sample: ", i))
  
  i.forward_reads <- paste0(i, "_trimmed.fq.gz")
  
  
  i.filtered_forward_reads <- paste0(i, "_filtered.fq.gz")
  
  # Quality trimming/filtering
  i.filtered_out <- filterAndTrim(i.forward_reads, i.filtered_forward_reads,
                              minLen=180,
                              # truncLen=c(290,200),
                              # trimLeft = c(5,10),
                              rm.phix = TRUE,
                              maxN = 0,
                              truncQ = 2,
                              maxEE=2,
                              compress = T, 
                              verbose = TRUE, 
                              multithread = TRUE)
  
  # Generating an error model
  i.err_forward_reads <- learnErrors(i.filtered_forward_reads, multithread=TRUE)
  
  
  # Dereplication
  i.derep_forward <- derepFastq(i.filtered_forward_reads, verbose=TRUE)
  names(i.derep_forward) <- i 
  

  # Inferring ASVs
  i.dada_forward <- dada(i.derep_forward, err=i.err_forward_reads, pool="pseudo", multithread=TRUE)
  
  
  
  # i.merged_amplicons <- mergePairs(i.dada_forward, i.derep_forward, 
  #                                  i.dada_reverse, i.derep_reverse, 
  #                                  # trimOverhang=TRUE, 
  #                                  # minOverlap=170,
  #                                  verbose = TRUE)
  
  # Generating a count table
  i.seqtab <- makeSequenceTable(i.dada_forward)


  # Chimera identification
  i.seqtab.nochim <- removeBimeraDenovo(i.seqtab, verbose=T) 
  

  # set a little function
  getN <- function(x) sum(getUniques(x))


  # making a data table for paired-end reads
  # i.summary_tab <- data.frame(row.names=i, dada2_input=i.filtered_out[,1],
  #                         filtered=i.filtered_out[,2], dada_f=sapply(i.dada_forward, getN),
  #                         dada_r=sapply(i.dada_reverse, getN), merged=sapply(i.merged_amplicons, getN),
  #                         nonchim=rowSums(i.seqtab.nochim),
  #                         final_perc_reads_retained=round(rowSums(i.seqtab.nochim)/i.filtered_out[,1]*100, 1))


  # making a data table for single reads
  i.track <- cbind(i, getN(i.filtered_out), getN(i.dada_forward), rowSums(i.seqtab), rowSums(i.seqtab.nochim))


  # save summary table into a txt file
  write.table(i.track, file = paste0(i, "_data_summary.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)
  

  # Extracting the standard goods from DADA2
  i.asv_seqs <- colnames(i.seqtab.nochim)
  i.asv_headers <- vector(dim(i.seqtab.nochim)[2], mode="character")
  
  
  for (j in 1:dim(i.seqtab.nochim)[2]) {
    i.asv_headers[j] <- paste(paste0(">", i), j, sep="_")
  }
  
  
  # save ASVs sequences
  i.asv_fasta <- c(rbind(i.asv_headers, i.asv_seqs))
  head(i.asv_fasta)
  write(i.asv_fasta, paste0(i, ".fasta"))
  
  
  # save ASVs count table
  i.asv_tab <- t(i.seqtab.nochim)
  
  row.names(i.asv_tab) <- sub(">", "", i.asv_headers)

  write.table(i.asv_tab, paste0(i, "_counts.tsv"), sep="\t", quote=F, col.names=NA)

  
  # Assigning taxonomy using IdTaxa command in DECIPHER packages
  i.dna <- DNAStringSet(getSequences(i.seqtab.nochim))
  
  i.tax_info <- IdTaxa(test=i.dna, trainingSet=trainingSet, strand="both", processors=NULL)
  
  ranks <- c("rootrank", "kingdom", "division" , "phylum", "class", "order", "family", "genus", "species")

  i.asv_tax <- t(sapply(i.tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(i.asv_tax) <- ranks
  
  rownames(i.asv_tax) <- gsub(pattern=">", replacement="", x=i.asv_headers)
  
  # save taxonomy results
  write.table(i.asv_tax, paste0(i, "_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)

  print (paste0(i, "of marine have done!"))
  print ("--------------------------------------------------------------------------")
  
}


############################################################################################################################################################################################################
############################################################################################################################################################################################################
# ASVs to OTUs
# To extract target ASV sequences for each species
# For two colpoda species
grep -E "Colpodea" Soil_taxonomy.tsv | cut -f 1 | sort -V > Soil.list

cat *_counts.tsv > Soil_counts.tsv

for i in $(cat Soil.list); do grep -w $i Soil_counts.tsv >> Soil_counts2.tsv & done


# For two marine speices
grep -E "Strombidiida" SS_taxonomy.tsv | cut -f 1 | sort -V > SS.list

grep -E "Euplotia" EV_taxonomy.tsv | cut -f 1 | sort -V > EV.list

cat EV.list SS.list > Marine.list

cat *_counts.tsv > Marine_counts.tsv

for i in $(cat Marine.list); do grep -w $i Marine_counts.tsv >> Marine_counts2.tsv & done


# Combine four species together
cat Soil.list Marine.list > All.list
cat Soil_counts2.tsv Marine_counts2.tsv > All_counts.tsv
sed -i '1i Representseq\tTotal' All_counts.tsv

# Concatenate all fasta together
cat Soil.fasta Marine.fasta > All_asv.fasta

# Extract sequences in regions contained in file 'All.list'
seqtk subseq All_asv.fasta All.list > All_sub.fasta


# Convert unique sequence into redundant sequence based on count table
# the deunique.seqs command is the reverse of the unique.seqs command
mothur
deunique.seqs(fasta=All_sub.fasta, count=All_counts.tsv)
#All_sub.redundant.fasta


# The sequence Variants were clustered using UCLUST at various thresholds ranging from 90%-100%.
for i in $(seq 0.9 0.01 1.0);do pick_otus.py -i All_sub.redundant.fasta -m uclust -s $i -o $i.All; done


# To generate the OTU table
for i in $(seq 0.9 0.01 1.0);do make_otu_table.py -i $i.All/All_sub.redundant_otus.txt -o All_sub.$i.biom & done


# To summarize the OTU table.
biom summarize-table -i All_sub.1.00.biom -o All_sub.100.biom.txt

# SCR28.6: 6.000
# SS16RN1: 13.000
# SCD28.7: 18.000
# SLD28.6: 60.000
# ICD18.5: 1,641.000
# SCR18.5: 2,110.000
# ICD28.1: 2,364.000
# SCD18.1: 2,893.000
# ICD28.3: 2,925.000
# SCD18.2: 3,011.000
# ICD28.4: 3,114.000
# ICD28.2: 7,040.000
# SCD28.4: 9,058.000
# SCR28.4: 11,574.000
# SCD18.5: 11,994.000
# ICD28.5: 14,252.000
# ICR18.4: 17,361.000
# SS21DN1: 21,094.000
# EV21RN2: 21,956.000
# ICD18.1: 24,276.000
# EV16RN1: 25,783.000
# SS21RN2: 26,925.000
# ICR18.5: 28,172.000
# EV21DN1: 29,027.000
# EV2516DN2: 29,162.000
# SS16DN1: 29,805.000
# SS16DN2: 29,823.000
# SS25RN2: 29,843.000
# SS2516DN2: 30,127.000
# SS25DN2: 31,890.000
# SS25RN1: 32,325.000
# ICD18.4: 32,604.000
# SS2516DN1: 33,348.000
# EV21DN2: 33,869.000
# EV2516D: 34,096.000
# SS2516RN2: 34,122.000
# SS25DN1: 34,467.000
# SS21DN2: 34,663.000
# SS2516RN1: 35,038.000
# EV2516RN1: 35,909.000
# EV25DN2: 36,096.000
# SS2516D: 36,181.000
# EV16DN2: 36,393.000
# EV16DN1: 36,706.000
# SS16RN2: 36,975.000
# SCR18.3: 37,213.000
# SS21D: 37,338.000
# EV25RN2: 37,375.000
# EV2516DN1: 37,791.000
# SS21RN1: 39,195.000
# EV2516RN2: 39,639.000
# EV25RN1: 39,758.000
# EV16RN2: 40,709.000
# EV2516R: 41,296.000
# SCR18.2: 41,346.000
# EV25DN1: 41,471.000
# EV16D: 41,603.000
# EV21RN1: 42,749.000
# EV21D: 42,814.000
# SS21R: 43,018.000
# SS25D: 43,190.000
# ICD18.2: 44,837.000
# EV25D: 45,593.000
# SS2516R: 46,151.000
# EV21R: 49,915.000
# SS16R: 51,788.000
# EV16R: 52,723.000
# ICR18.1: 54,220.000
# SS25R: 54,224.000
# ICR18.3: 54,537.000
# EV25R: 56,241.000
# SCD28.1: 57,875.000
# SLD28.2: 58,998.000
# ICR28.1: 61,701.000
# ICR28.2: 62,700.000
# ILD18.2: 62,964.000
# SCR28.5: 63,788.000
# ICD18.3: 64,193.000
# ICR28.5: 65,936.000
# ILD28.3: 67,307.000
# ILD18.3: 68,749.000
# ILD28.4: 69,257.000
# ILD28.2: 69,752.000
# SLR28.3: 71,087.000
# SLD18.2: 71,205.000
# SLR28.5: 71,544.000
# SLR18.2: 71,744.000
# ILR18.3: 71,832.000
# ICR18.2: 73,255.000
# ILD28.5: 73,328.000
# SLR18.1: 73,348.000
# SLD28.3: 73,362.000
# SCR28.3: 73,692.000
# SLD18.1: 74,167.000
# SLR18.3: 74,462.000
# SCR18.1: 74,970.000
# ICR28.4: 76,207.000
# ICR28.3: 76,970.000
# ILD18.4: 77,090.000
# SLR28.4: 77,370.000
# ILD28.1: 77,543.000
# ILR28.5: 78,272.000
# SLR18.4: 78,276.000
# SCD28.5: 78,417.000
# SLD28.5: 78,712.000
# SLD18.5: 78,757.000
# SLR18.5: 79,787.000
# SLR28.2: 79,839.000
# SLD18.3: 80,229.000
# ILR18.4: 80,301.000
# SCR18.4: 80,402.000
# ILR28.1: 82,284.000
# SLD18.4: 82,345.000
# ILR18.5: 83,506.000
# SCR28.2: 83,738.000
# SLD28.1: 84,215.000
# ILR28.3: 84,288.000
# ILR18.2: 84,421.000
# ILR28.4: 84,973.000
# SLR28.1: 85,086.000
# ILD18.1: 85,310.000
# ILR28.2: 87,222.000
# ILR18.1: 87,584.000


# Rarefied based on the lowest number for each species.
for i in $(seq 0.9 0.01 1.0);do multiple_rarefactions.py -i All_sub.$i.biom -m 100 -x 1600 -s 100 -n 20 -o $i.rarefied/ & done


# Calculating the Alpha-diversity estimator (observed_otus).
for i in $(seq 0.9 0.01 1.0);do alpha_diversity.py -i $i.rarefied/ -m observed_otus -o $i.alpha/ & done


# Put them all together.
for i in $(seq 0.9 0.01 1.0);do collate_alpha.py -i $i.alpha/ -o $i.alpha_collated/ & done


# Copy and rename the results files.
for i in $(seq 0.9 0.01 1.0);do mv $i.alpha_collated/observed_otus.txt $i.txt ;done


############################################################################################################################################################################################################
############################################################################################################################################################################################################
# To visualize

library(readxl)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(Rsongbao)

setwd("G:/[项目] 2种Colpoda拷贝数实验/Submission/Colpoda_DADA2/OTU_R1600-0612")

SE_summary <- function (data, varname, groupnames) 
{
  require(dplyr)
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = sprintf("%0.2f",mean(x[[col]], na.rm = TRUE)), 
      se = sprintf("%0.2f",sd(x[[col]]/sqrt(length(x[[col]])), na.rm = TRUE)))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func, 
                    varname)
  data_sum <- rename(data_sum, c(mean = varname))
  return(data_sum)
}



# Mean ---------------------------------------------------------------------

all_mean <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>%
  select(-Species) %>% 
  column_to_rownames("Sample") %>%
  group_by(Group) %>%
  summarise_all(mean)

write.csv(all_mean,"mean_for_heatmap.csv")



# Mean±SE ------------------------------------------------------------------

X90 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X90",digtal = 2) %>% 
  unite(col = "X90",
        X90,se,
        sep = " ± ")


X91 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X91",digtal = 1) %>% 
  unite(col = "X91",
        X91,se,
        sep = " ± ")


X92 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X92",digtal = 1) %>% 
  unite(col = "X92",
        X92,se,
        sep = " ± ")


X93 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X93",digtal = 1) %>% 
  unite(col = "X93",
        X93,se,
        sep = " ± ")


X94 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X94",digtal = 1) %>% 
  unite(col = "X94",
        X94,se,
        sep = " ± ")


X95 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X95",digtal = 1) %>% 
  unite(col = "X95",
        X95,se,
        sep = " ± ")

X96 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X96",digtal = 1) %>% 
  unite(col = "X96",
        X96,se,
        sep = " ± ")


X97 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X97",digtal = 1) %>% 
  unite(col = "X97",
        X97,se,
        sep = " ± ")

X98 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X98",digtal = 1) %>% 
  unite(col = "X98",
        X98,se,
        sep = " ± ")


X99 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X99",digtal = 1) %>% 
  unite(col = "X99",
        X99,se,
        sep = " ± ")


X100 <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(-Sample) %>% 
  SE_summary(groupnames = "Group",varname = "X100",digtal = 2) %>% 
  unite(col = "X100",
        X100,se,
        sep = " ± ")


otu_all <- cbind(X100,X99,X98,X97,X96,X95,X94,X93,X92,X91,X90)

write.csv(otu_all,"for_heatmap_label.csv")


# OTU heatmap  ----------------------------------------------------------------------

library(pheatmap)
library(RColorBrewer)

Color<-c(brewer.pal(n=10,name = "RdYlBu")[11:7],
         brewer.pal(n=11,name = "RdYlBu")[5:1])


read_excel("R1600.xlsx",sheet = "For_heatmap") %>%
  column_to_rownames("Group") %>%
  pheatmap(display_numbers = F,
           number_format = "%.1f",
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           cellwidth =45,cellheight = 15,
           color = colorRampPalette(Color)(100))



# ASV_Boxplot -------------------------------------------------------------

Hap_CI <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(Species,Group,X99) %>%
  filter(Species=="CI") %>%
  SE_summary(varname = "X99", groupnames = c("Species","Group"),digtal = 2) %>% 
  ggplot(aes(x=Group, y = X99, shape= Group))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= X99-se, ymax = X99+se),width = 0.1)+
  scale_shape_manual(values = c(18,18,18,18,0,0,0,0))+
  scale_x_discrete(limits=c("ILD18", "ILD28", "ICD18", "ICD28", "ILR18", "ILR28", "ICR18", "ICR28"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))


Hap_CS <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(Species,Group,X99) %>%
  filter(Species=="CS") %>%
  SE_summary(varname = "X99", groupnames = c("Species","Group"),digtal = 2) %>% 
  ggplot(aes(x=Group, y = X99, shape= Group))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= X99-se, ymax = X99+se),width = 0.1)+
  scale_shape_manual(values = c(18,18,18,18,0,0,0,0))+
  scale_x_discrete(limits=c("SLD18", "SLD28","SCD18", "SCD28", "SLR18", "SLR28","SCR18", "SCR28"))+
  theme_base()+
  ylab("Haplotype")+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))


Hap_EV <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(Species,Group,X100) %>%
  filter(Species=="EV") %>%
  SE_summary(varname = "X100", groupnames = c("Species","Group"),digtal = 2) %>% 
  ggplot(aes(x=Group, y = X100, shape= Group))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= X100-se, ymax = X100+se),width = 0.1)+
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0))+
  scale_x_discrete(limits=c("EV16D", "EV21D","EV25D", "EV2516D", "EV16R", "EV21R","EV25R", "EV2516R"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))


Hap_SS <- read_excel("R1600.xlsx",sheet = "ASV90_100") %>% 
  select(Species,Group,X100) %>%
  filter(Species=="SS") %>%
  SE_summary(varname = "X100", groupnames = c("Species","Group"),digtal = 2) %>% 
  ggplot(aes(x=Group, y = X100, shape= Group))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= X100-se, ymax = X100+se),width = 0.1)+
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0))+
  scale_x_discrete(limits=c("SS16D", "SS21D","SS25D", "SS2516D", "SS16R", "SS21R","SS25R", "SS2516R"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))

ggarrange(Hap_CS, Hap_EV, Hap_CI, Hap_SS, nrow = 1, ncol = 4)


