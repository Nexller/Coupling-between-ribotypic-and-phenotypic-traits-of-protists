# PART1. QIIME pipeline for single-clone high-throughput sequencing data.

# Unzip the raw data files.
unzip 18S211375_1.zip
unzip 18S211375_2.zip


# Prepare the primer sequences files.
vi Colpoda.oligos
CCAGCASCYGCGGTAATTCC
ACTTTCGTTCTTGATYRA

vi Marine.oligos
GAADCTGYGAAYGGCTC
ACCAGACTTGCCCTCC


# list all fastq files, cut the sample name for soil and marine species seperately, the two sample lists will used in the following looping.
ls C*.raw.fastq | cut -d . -f1 > sample_soil
ls E*.raw.fastq && ls S*.raw.fastq | cut -d . -f1 > sample_marine


# Quality control with Prinseq-lite
for i in $(cat sample_soil); do perl prinseq-lite.pl -verbose -fastq $i.raw.fastq -out_good $i.fastq -min_len 380 -ns_max_n 0 -min_qual_mean 20; done 
for i in $(cat sample_marine); do perl prinseq-lite.pl -verbose -fastq $i.raw.fastq -out_good $i.fastq -min_len 400 -ns_max_n 0 -min_qual_mean 20; done


# Rename the filtered data
for i in $(cat sample_soil); do rename.sh in=$i.fastq out=$i.fasta prefix=$i; done
for i in $(cat sample_marine); do rename.sh in=$i.fastq out=$i.fasta prefix=$i; done


# Concatenate fasta files for each species.
cat CI*.fasta > CI.fasta
cat CS*.fasta > CS.fasta
cat EV*.fasta > EV.fasta
cat SS*.fasta > SS.fasta


# Creat a new list contain four species names.
ls *.fasta | cut -d . -f1 > species



# Vsearch was used to detected the chimera sequences based on silva database.
# Taxonomic assignments of Non-chimeras sequences using uclust.
for i in $(cat species)
do
  vsearch --threads 20 --uchime_ref $i.fasta --db silva_132_97_18S.fna --nonchimeras $i.nonchimeras.fasta
  assign_taxonomy.py -i $i.nonchimeras.fasta -r silva_132_97_18S.fna -t consensus_taxonomy_all_levels.txt \
  -m uclust --min_consensus_fraction 0.8 -o $i.silva_132/
done



# The reads assigned to each species were extracted by mothur.
mothur > get.lineage(fasta=CI.nonchimeras.fasta,taxonomy=CI.nonchimeras_tax_assignments.txt,taxon=D_6__Colpodea)
mothur > get.lineage(fasta=CS.nonchimeras.fasta,taxonomy=CS.nonchimeras_tax_assignments.txt,taxon=D_6__Colpodea)
mothur > get.lineage(fasta=EV.nonchimeras.fasta,taxonomy=EV.nonchimeras_tax_assignments.txt,taxon=D_6__Euplotia)
mothur > get.lineage(fasta=SS.nonchimeras.fasta,taxonomy=SS.nonchimeras_tax_assignments.txt,taxon=D_6__Oligotrichia)


# To trimming the primer sequences.
# 
mothur > pcr.seqs(fasta=CI.nonchimeras.pick.fasta, oligos=Colpoda.oligos, processors=46, pdiffs=1)
mothur > pcr.seqs(fasta=CS.nonchimeras.pick.fasta, oligos=Colpoda.oligos, processors=46, pdiffs=1)
mothur > pcr.seqs(fasta=EV.nonchimeras.pick.fasta, oligos=Marine.oligos, processors=46, pdiffs=1)
mothur > pcr.seqs(fasta=SS.nonchimeras.pick.fasta, oligos=Marine.oligos, processors=46, pdiffs=1)


# Homopolymers >= 6 were discard.
screen.seqs(fasta=CI.nonchimeras.pick.pcr.fasta, maxhomop=6, processors=46)
screen.seqs(fasta=CS.nonchimeras.pick.pcr.fasta, maxhomop=6, processors=46)
screen.seqs(fasta=EV.nonchimeras.pick.pcr.fasta, maxhomop=6, processors=46)
screen.seqs(fasta=SS.nonchimeras.pick.pcr.fasta, maxhomop=6, processors=46)


# Singletons were detected at 100% cutoff, and removed manually using the following four scripts.
for i in $(cat species); do pick_otus.py -i $i.nonchimeras.pick.pcr.good.fasta -m uclust -s 1.00 -o $i.singleton_100 --threads 48& done
awk NF==2 SS.nonchimeras.pick.pcr.good_otus.txt | cut -f2 > SS.Singleton.ID.txt
for i in $(cat species); do awk NF==2 $i.nonchimeras.pick.pcr.good_otus.txt | cut -f2 > $i.Singleton.ID.txt; done
for i in $(cat species); do filter_fasta.py -f $i.nonchimeras.pick.pcr.good.fasta -o $i.nonchimeras_nosingleton.fasta -s $i.Singleton.ID.txt -n; done



# The high-quality sequences were clustered using UCLUST at various thresholds ranging from 90%-100%.
for i in $(seq 0.9 0.01 1.0)
do
  pick_otus.py -i CI.nonchimeras_nosingleton.fasta -m uclust -s $i -o $i.CI
  pick_otus.py -i CS.nonchimeras_nosingleton.fasta -m uclust -s $i -o $i.CS
  pick_otus.py -i EV.nonchimeras_nosingleton.fasta -m uclust -s $i -o $i.EV
  pick_otus.py -i SS.nonchimeras_nosingleton.fasta -m uclust -s $i -o $i.SS
done



# To generate the OTU table
for i in $(seq 0.9 0.01 1.0)
do
  make_otu_table.py -i $i.CI/CI.nonchimeras_nosingleton_otus.txt -o $i.CI/CI.$i.biom
  make_otu_table.py -i $i.CS/CS.nonchimeras_nosingleton_otus.txt -o $i.CS/CS.$i.biom
  make_otu_table.py -i $i.EV/EV.nonchimeras_nosingleton_otus.txt -o $i.EV/EV.$i.biom
  make_otu_table.py -i $i.SS/SS.nonchimeras_nosingleton_otus.txt -o $i.SS/SS.$i.biom
done


# Counts/sample detail:
# CI3: 15,805.000
# CI1: 17,871.000
# CI2: 19,897.000
# CI4: 23,751.000
# CI5: 26,060.000

# EV2: 22,076.000
# EV3: 23,332.000
# EV4: 24,146.000
# EV1: 28,141.000
# EV5: 33,554.000

# CS5: 19,926.000
# CS2: 30,564.000
# CS4: 30,934.000
# CS1: 33,690.000
# CS3: 35,845.000

# SS4: 10.000
# SS3: 74.000
# SS1: 98.000
# SS2: 98.000
# SS5: 113.000

# Rarefied based on the lowest number for each species.
for i in $(seq 0.9 0.01 1.0)
do
  multiple_rarefactions.py -i $i.CI/CI.$i.biom -m 1000 -x 15000 -s 1000 -n 1 -o $i.CI/$i.rarefied/
  multiple_rarefactions.py -i $i.CS/CS.$i.biom -m 1000 -x 19000 -s 1000 -n 1 -o $i.CS/$i.rarefied/
  multiple_rarefactions.py -i $i.EV/EV.$i.biom -m 1000 -x 22000 -s 1000 -n 1 -o $i.EV/$i.rarefied/
  multiple_rarefactions.py -i $i.SS/SS.$i.biom -m 1 -x 10 -s 1 -n 1 -o $i.SS/$i.rarefied/
done



# Calculating the Alpha-diversity estimator (observed_otus).
for i in $(seq 0.9 0.01 1.0)
do
  alpha_diversity.py -i $i.CI/$i.rarefied/ -m observed_otus -o $i.CI/alpha/
  alpha_diversity.py -i $i.CS/$i.rarefied/ -m observed_otus -o $i.CS/alpha/
  alpha_diversity.py -i $i.EV/$i.rarefied/ -m observed_otus -o $i.EV/alpha/
  alpha_diversity.py -i $i.SS/$i.rarefied/ -m observed_otus -o $i.SS/alpha/
done



# Concatenate all.
for i in $(seq 0.9 0.01 1.0)
do 
  collate_alpha.py -i $i.CI/alpha/ -o $i.CI/alpha_collated/
  collate_alpha.py -i $i.CS/alpha/ -o $i.CS/alpha_collated/
  collate_alpha.py -i $i.EV/alpha/ -o $i.EV/alpha_collated/
  collate_alpha.py -i $i.SS/alpha/ -o $i.SS/alpha_collated/
done



# Copy and rename the results files.
for i in $(seq 0.9 0.01 1.0)
do
  mv $i.CI/alpha_collated/observed_otus.txt OTU_Clustering/$i.CI.txt
  mv $i.CS/alpha_collated/observed_otus.txt OTU_Clustering/$i.CS.txt
  mv $i.EV/alpha_collated/observed_otus.txt OTU_Clustering/$i.EV.txt
  mv $i.SS/alpha_collated/observed_otus.txt OTU_Clustering/$i.SS.txt
done


# Extract the observed_otus results.
for i in $(ls *.txt); do awk 'NR == 1; END{print}' $i | cut -f 4-8 > $i.tsv; done



###########################################################################################################################################################################
# PART2 DADA2 pipeline for single-clone high-throughput sequencing data.

# Creating a new file and prepare the single-clone raw reads and sample list file. 
# The following scripts were running in linux platform.
mkdir DADA2

cp *.raw.fastq DADA2/
cp sample_marine DADA2/
cp sample_soil DADA2/
cat sample_marine cp sample_soil > samples
 

# Trimming primer sequences using cutadapt v1.18.
# For two soil species
for sample in $(cat sample_soil)
do
  echo "Processing sample: $sample"
  cutadapt -a ^CCAGCASCYGCGGTAATTCC...ACTTTCGTTCTTGATYRA$ -m 50 -o ${sample}_raw_trimmed.fq ${sample}.raw.fastq \
  >> cutadapt_soil_stats.txt 2>&1
done


# For two marine species
for sample in $(cat sample_marine)
do
  echo "Processing sample: $sample"
  cutadapt -a ^GAADCTGYGAAYGGCTC...ACCAGACTTGCCCTCC$ \
  -m 50 -o ${sample}_raw_trimmed.fq ${sample}.raw.fastq \
  >> cutadapt_marine_stats.txt 2>&1
done


###############################################################################
# The following analysis will switch R and start using DADA2 in Rstudio-sever.
setwd("/database2/Zou_Colpoda_polymorphisms/Copoda_clone/DADA2")


## Downloading DECIPHER-formatted SILVA v138 reference 
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")


# Install and load R packages.
# BiocManager::install('dada2')
library(dada2) 
library(DECIPHER)
# packageVersion("DECIPHER") #2.14.0


# load the reference database.
load("SILVA_SSU_r138_2019.RData")


# samples contains all of filenames 
samples <- scan("samples", what="character")


# Looping for each sample using DADA2 pipeline;
# filterAndTrim, was used to quality trim/filter;
# learnErrors, was used to generate an error model of our data;
# derepFastq, was used to dereplicate sequences;
# dada, the core algorithm, to infer ASVs on both forward and reverse reads independently;
# mergePairs, to merge forward and reverse reads to further refine ASVs;
# makeSequenceTable, to generate a count table;
# removeBimeraDenovo, to screen for and remove chimeras;
# IdTaxa, to assign taxonomy.

for (i in samples){
  
  print (i)
  
  i.forward_reads <- paste0(i, "_raw_trimmed.fq")
  # i.reverse_reads <- paste0(i, "_R2_trimmed.fq.gz")
  
  
  i.filtered_forward_reads <- paste0(i, "_filtered.fq.gz")
  # i.filtered_reverse_reads <- paste0(i, "_R2_filtered.fq.gz")
  
  
  filtered_out <- filterAndTrim(i.forward_reads, i.filtered_forward_reads,
                                # i.reverse_reads, i.filtered_reverse_reads, 
                                minLen=180,
                                rm.phix = TRUE,
                                maxN = 0,
                                truncQ = 2,
                                maxEE = 2,
                                # maxEE = c(2,2),
                                compress = T, 
                                verbose = TRUE, 
                                multithread = TRUE)
  
  
  i.err_forward_reads <- learnErrors(i.filtered_forward_reads, multithread=TRUE)
  # i.err_reverse_reads <- learnErrors(i.filtered_reverse_reads, multithread=TRUE)
  
  
  
  i.derep_forward <- derepFastq(i.filtered_forward_reads, verbose=TRUE)
  names(i.derep_forward) <- i 
  
  # i.derep_reverse <- derepFastq(i.filtered_reverse_reads, verbose=TRUE)
  # names(i.derep_reverse) <- i
  
  
  
  i.dada_forward <- dada(i.derep_forward, err=i.err_forward_reads, pool="pseudo", multithread=TRUE)
  # i.dada_reverse <- dada(i.derep_reverse, err=i.err_reverse_reads, pool="pseudo", multithread=TRUE)
  
  
  
  # i.merged_amplicons <- mergePairs(i.dada_forward, i.derep_forward, 
  #                                  i.dada_reverse, i.derep_reverse, 
  #                                  # trimOverhang=TRUE, 
  #                                  # minOverlap=170,
  #                                  verbose = TRUE)
  
  i.seqtab <- makeSequenceTable(i.dada_forward)
  class(i.seqtab) # matrix
  dim(i.seqtab) 
  
  
  i.seqtab.nochim <- removeBimeraDenovo(i.seqtab, verbose=T) 
  
  
  i.dna <- DNAStringSet(getSequences(i.seqtab.nochim))
  
  
  i.tax_info <- IdTaxa(test=i.dna, trainingSet=trainingSet, strand="both", processors=NULL)
  
  
  i.asv_seqs <- colnames(i.seqtab.nochim)
  i.asv_headers <- vector(dim(i.seqtab.nochim)[2], mode="character")
  
  
  # Ñ­»·¸ÄÃû
  for (j in 1:dim(i.seqtab.nochim)[2]) {
    i.asv_headers[j] <- paste(paste0(">", i), j, sep="_")
  }
  
  
  i.asv_fasta <- c(rbind(i.asv_headers, i.asv_seqs))
  head(i.asv_fasta)
  write(i.asv_fasta, paste0(i, ".fasta"))
  
  
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  i.asv_tax <- t(sapply(i.tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(i.asv_tax) <- ranks
  
  #¼ÓÉÏÐÐµÄÃû³Æ
  rownames(i.asv_tax) <- gsub(pattern=">", replacement="", x=i.asv_headers)
  
  write.table(i.asv_tax, paste0(i, "_taxonomy.tsv"), sep = "\t", quote=F, col.names=NA)
  
  
}



# Concatenate the dada2 result files.
cat *_taxonomy.tsv > All_taxonomy.tsv
cat *fasta > All.fasta


# # The reads assigned to each species were extracted using seqtk. 
grep -E "Colpodea|Euplotes|Oligotrichia" All_taxonomy.tsv | cut -f 1 | sort -V > Taxo.list
cut -f 1 All2_taxonomy.tsv | sort -V > Taxo.list
seqtk subseq All.fasta Taxo.list > All_sub.fasta


# The sequence Variants were clustered using UCLUST at various thresholds ranging from 90%-100%.
for i in $(seq 0.9 0.01 1.0);do pick_otus.py -i All_sub.fasta -m uclust -s $i -o $i.All_clone; done


# To generate the OTU table
for i in $(seq 0.9 0.01 1.0);do make_otu_table.py -i $i.All_clone/All_sub_otus.txt -o All_sub.$i.biom & done


# To summarize the OTU table.
biom summarize-table -i All_sub.1.00.biom -o All_sub.100.biom.txt

# Counts/sample detail:
# EV3: 1.000
# SS4: 1.000
# EV4: 1.000
# SS5: 1.000
# SS3: 1.000
# EV2: 1.000
# EV1: 1.000
# EV5: 1.000
# CS2: 12.000
# CS5: 15.000
# CI2: 16.000
# CS1: 19.000
# CS3: 20.000
# CI3: 24.000
# CI1: 25.000
# CS4: 25.000
# CI4: 26.000
# CI5: 28.000


# Rarefied based on the lowest number for each species.
for i in $(seq 0.9 0.01 1.0);do multiple_rarefactions.py -i All_sub.$i.biom -m 1 -x 12 -s 1 -n 1 -o $i.rarefied/ & done


# Calculating the Alpha-diversity estimator (observed_otus).
for i in $(seq 0.9 0.01 1.0);do alpha_diversity.py -i $i.rarefied/ -m observed_otus -o $i.alpha/ & done


# Put them all together.
for i in $(seq 0.9 0.01 1.0);do collate_alpha.py -i $i.alpha/ -o $i.alpha_collated/ & done


# Copy and rename the results files.
for i in $(seq 0.9 0.01 1.0);do mv $i.alpha_collated/observed_otus.txt $i.txt ;done



################################################################################ SNP

# Build mapping index using BWA and samtools
for i in $(cat species)
do 
	bwa index -a is $i.fasta
	samtools faidx $i.fasta
done


# to calculate the SNP, ASVs of each sample was aligned against the dominant sequence for each species.
# Colpoda inflata
for i in $(cat CI.list)
do
	bwa mem CI_18Seqs.fasta $i.fasta > $i.sam
	samtools view -bhS $i.sam -o $i.bam
	samtools sort $i.bam $i.sort
	samtools index $i.sort.bam
	samtools mpileup -Q 20 -q 20 -s -f CI_18Seqs.fasta $i.sort.bam > $i.realigned.mpileupBwa
	perl printBasesBwa2.pl $i.realigned.mpileupBwa > $i.printBasesBwa
done

# Colpoda steinii
for i in $(cat CS.list)
do
	bwa mem CS_18Seqs.fasta $i.fasta > $i.sam
	samtools view -bhS $i.sam -o $i.bam
	samtools sort $i.bam $i.sort
	samtools index $i.sort.bam
	samtools mpileup -Q 20 -q 20 -s -f CS_18Seqs.fasta $i.sort.bam > $i.realigned.mpileupBwa
	perl printBasesBwa2.pl $i.realigned.mpileupBwa > $i.printBasesBwa
done


# Euplotes vannus
for i in $(cat EV.list)
do
	bwa mem EV_18Seqs.fasta $i.fasta > $i.sam
	samtools view -bhS $i.sam -o $i.bam
	samtools sort $i.bam $i.sort
	samtools index $i.sort.bam
	samtools mpileup -Q 20 -q 20 -s -f EV_18Seqs.fasta $i.sort.bam > $i.realigned.mpileupBwa
	perl printBasesBwa2.pl $i.realigned.mpileupBwa > $i.printBasesBwa
done


# Strombidium sulcatum
for i in $(cat SS.list)
do
	bwa mem SS_18Seqs.fasta $i.fasta > $i.sam
	samtools view -bhS $i.sam -o $i.bam
	samtools sort $i.bam $i.sort
	samtools index $i.sort.bam
	samtools mpileup -Q 20 -q 20 -s -f SS_18Seqs.fasta $i.sort.bam > $i.realigned.mpileupBwa
	perl printBasesBwa2.pl $i.realigned.mpileupBwa > $i.printBasesBwa
done


#######################################################
# To visualize the SNP results, analysis was done in R.
library(tidyverse)
library(readxl)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(patchwork)
setwd("C:/Users/songbaozou/Desktop/Colpoda_DADA2")

filelist = list.files(pattern="*printBasesBwa")
results <- list()
for(i in filelist) {
  BWA <- read.table(i) 
  colnames(BWA) <- c("Ref","Pos","Base","Coverage","A","C","G","T","X","Y")  

  list <- select(BWA,3,5:8) %>% 
    group_by(Base) %>% 
    summarise_all(sum) %>% 
    gather(-Base, key = "type",value = "value")

  results[[i]] <- list 
  openxlsx::write.xlsx(results,"SNP_results.xlsx")
}



# Boxplot -----------------------------------------------------------------

SE_summary <- function(data, varname, groupnames, digtal){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = round(mean(x[[col]], na.rm=TRUE),digits = digtal),
      se = round(sd(x[[col]]/sqrt(length(x[[col]])), na.rm=TRUE),digits = digtal))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



Mutation_CI <- read_excel("SNP_all.xlsx",sheet = "Mutation_rate",range = "A1:E124") %>%
  filter(Species=="CI") %>%
  SE_summary(varname = "Mutate_Rate100", 
             groupnames = c("Species","Treatments"), digtal = 2) %>% 
  ggplot(aes(x=Treatments, y = Mutate_Rate100, shape= Treatments))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= Mutate_Rate100-se, ymax = Mutate_Rate100+se),width = 0.1)+
  scale_shape_manual(values = c(18,18,18,18,0,0,0,0))+
  scale_x_discrete(limits=c("ILD18", "ILD28", "ICD18", "ICD28", "ILR18", "ILR28", "ICR18", "ICR28"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))



Mutation_CS <- read_excel("SNP_all.xlsx",sheet = "Mutation_rate",range = "A1:E124") %>%
  filter(Species=="CS") %>%
  SE_summary(varname = "Mutate_Rate100", 
             groupnames = c("Species","Treatments"), digtal = 2) %>% 
  ggplot(aes(x=Treatments, y = Mutate_Rate100, shape= Treatments))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= Mutate_Rate100-se, ymax = Mutate_Rate100+se),width = 0.1)+
  scale_shape_manual(values = c(18,18,18,18,0,0,0,0))+
  scale_x_discrete(limits=c("SLD18", "SLD28","SCD18", "SCD28", "SLR18", "SLR28","SCR18", "SCR28"))+
  theme_base()+
  ylab("Mutationlotype")+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))



Mutation_EV <- read_excel("SNP_all.xlsx",sheet = "Mutation_rate",range = "A1:E124") %>%
  filter(Species=="EV") %>%
  SE_summary(varname = "Mutate_Rate100", 
             groupnames = c("Species","Treatments"), digtal = 2) %>%
  ggplot(aes(x=Treatments, y = Mutate_Rate100, shape= Treatments))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= Mutate_Rate100-se, ymax = Mutate_Rate100+se),width = 0.1)+
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0))+
  scale_x_discrete(limits=c("EV16D", "EV21D","EV25D", "EV2516D", "EV16R", "EV21R","EV25R", "EV2516R"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))
  


Mutation_SS <- read_excel("SNP_all.xlsx",sheet = "Mutation_rate",range = "A1:E124") %>%
  filter(Species=="SS") %>%
  SE_summary(varname = "Mutate_Rate100", 
             groupnames = c("Species","Treatments"), digtal = 2) %>%
  ggplot(aes(x=Treatments, y = Mutate_Rate100, shape= Treatments))+
  geom_point(size = 5.0)+
  geom_errorbar(aes(ymin= Mutate_Rate100-se, ymax = Mutate_Rate100+se),width = 0.1)+
  scale_shape_manual(values = c(0,0,0,0,0,0,0,0))+
  scale_x_discrete(limits=c("SS16D", "SS21D","SS25D", "SS2516D", "SS16R", "SS21R","SS25R", "SS2516R"))+
  theme_base()+
  geom_vline(xintercept = 4.5, linetype="dashed", size= 1)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",size = 11))


ggarrange(Mutation_CS, Mutation_EV, Mutation_CI, Mutation_SS,
          nrow = 2, ncol = 2)



# SNP spectrum ------------------------------------------------------------

Ts_CT <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Ts(C/T)", groupnames = "Group", digtal = 1)


Ts_AG <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Ts(A/G)", groupnames = "Group", digtal = 1)


Tv_TA <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Tv(T/A)", groupnames = "Group", digtal = 1)


Tv_GT <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Tv(G/T)", groupnames = "Group", digtal = 1)


Tv_AC <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Tv(A/C)", groupnames = "Group", digtal = 1)


Tv_GC <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Tv(G/C)", groupnames = "Group", digtal = 1)


Ts_Tv <- read_excel("SNP_all.xlsx",sheet = "SNP_spectrum",range = "A1:I124") %>%
  select(-Samples) %>% 
  SE_summary(varname = "Ts : Tv", groupnames = "Group", digtal = 1)


lst <- bind_cols(Ts_CT, Ts_AG, Tv_TA, Tv_GT, Tv_AC, Tv_GC, Ts_Tv)

openxlsx::write.xlsx(lst, "SNP_spectrum.xlsx")


# heatmap ----------------------------------------------------------------------
library(pheatmap)
library(RColorBrewer)

Color<-c(brewer.pal(n=10,name = "RdYlBu")[11:10], 
         brewer.pal(n=11,name = "RdYlBu")[5:1])

read_excel("SNP_all.xlsx",sheet = "SNP_spectrum1",range = "A1:H33") %>%
  column_to_rownames(var = "Group") %>% 
  pheatmap(display_numbers = F,
           number_format = "%.1f",
           fontsize = 12,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           cellwidth =45,
           cellheight = 15,
           color = colorRampPalette(Color)(100))


# Density -----------------------------------------------------------------

CS_RNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Colpoda_steinii") %>%
  ggplot(aes(x=Pos, y = RNA_102)) +
  geom_col(color = "blue")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())


CI_RNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Colpoda_inflata") %>%
  ggplot(aes(x=Pos, y = RNA_102)) +
  geom_col(color = "blue")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())

EV_RNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Euplotes") %>%
  ggplot(aes(x=Pos, y = RNA_102)) +
  geom_col(color = "blue")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())


SS_RNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Strombidium") %>%
  ggplot(aes(x=Pos,y = RNA_102)) +
  geom_col(color = "blue")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  xlab("Position")+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black",size = 12))



CS_DNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Colpoda_steinii") %>%
  ggplot(aes(x=Pos, y = DNA_102)) +
  geom_col(color = "brown")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())



CI_DNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Colpoda_inflata") %>%
  ggplot(aes(x=Pos, y = DNA_102)) +
  geom_col(color = "brown")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())

EV_DNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Euplotes") %>%
  ggplot(aes(x=Pos, y = DNA_102)) +
  geom_col(color = "brown")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title = element_blank())


SS_DNA <- read_excel("SNP_all.xlsx",sheet = "SNP_density") %>%
  filter(Species == "Strombidium") %>%
  ggplot(aes(x=Pos,y = DNA_102)) +
  geom_col(color = "brown")+
  geom_hline(aes(yintercept = 0),size = 1.2)+
  theme_minimal()+
  xlab("Position")+
  theme(axis.text = element_text(color = "black",size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black",size = 12))


CS_DNA + CS_RNA + plot_layout(ncol = 1, nrow = 2)
CI_DNA + CI_RNA + plot_layout(ncol = 1, nrow = 2)
EV_DNA + EV_RNA + plot_layout(ncol = 1, nrow = 2)
SS_DNA + SS_RNA + plot_layout(ncol = 1, nrow = 2)

 
################################################################################



