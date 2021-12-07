
######################################################################################################


# 0 # set working directory & load libs

#install.packages("vcd")
library(vcd)

#install.packages("summarytools")
library(summarytools)

library(readr)
library(splitstackshape)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)

source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/" # local 


###########################################################################################


# 1 # prepares metadata files to feed guppy_epca.sh:
# one group = one time point : 1 breed : max 2 days diff in bdays

# load metadata 

mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

####
# date formatting: 
from = c("2017-01-30",
         "2017-01-31","2017-02-01",
         "2017-02-03",
         "2017-02-06","2017-02-07","2017-02-08",
         "2017-02-10",
         "2017-02-14",
         "2017-02-16","2017-02-17",
         "2017-02-21",
         "2017-02-24", 
         "2017-02-28",
         "2017-03-03",
         "2017-03-06","2017-03-07","2017-03-08","2017-03-09","2017-03-10",
         "2017-08-14", #mock community
         "2018-01-24",  #probiotics - pos controls
         NA) # neg controls

to = c("tM",
       "t0","t0",
       "t1",
       "t2","t2","t2",
       "t3",
       "t4",
       "t5","t5",
       "t6",
       "t7", 
       "t8",
       "t9",
       "t10","t10","t10","t10","t10",
       "tNONE", #mock community
       "tNONE", #probiotics - pos controls
       "tNONE") #neg controls

# replace collection dates (date format) with groups of collection dates (character format)
mdat$collection_date <- plyr::mapvalues(as.character(mdat$collection_date), from, to)
unique(mdat$collection_date)
####

# load breed and bday data 
details <- read_excel(paste0(source_data, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'isolation_source'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$isolation_source <- gsub("G","",details$isolation_source)
details$isolation_source <- gsub("T","",details$isolation_source)

details <- details %>%
  dplyr::select(isolation_source,BIRTH_DAY,breed)

##########################################################################################
###########################################################################################


# POSITIVE CONTROLS 

# filter to keep only pos controls
mdat_sel <- mdat %>% 
  filter(Cohort=="MockCommunity"|Cohort=="PosControl_D-scour"|Cohort=="PosControl_ColiGuard") %>% 
  dplyr::select(isolation_source,DNA_plate,DNA_well)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

mdat_sel$isolation_source
# we have 9 Mock , 8 Protexin, 8 ColiGuard = 25 in total 

mdat_sel <- as.character(mdat_sel$ID)

writeLines(unlist(mdat_sel), paste0(out_dir,"guppy_input/pos_controls_sel.txt"), sep = " ")
# contains 25 sample IDs in total 


###########################################################################################
###########################################################################################


# ALL PIGGIES 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

df1 <- setDT(mdat_sel)[, .(Freq = .N), by = .(Cohort)]
df1[order(df1$Cohort)]
# we have 142 ctrl, 133 d-scour, 131 coliguard, 150 neo, 139 neo+d-scour, 130 neo+coliguard

mdat_sel <- as.character(mdat_sel$ID)

NROW(mdat_sel)
# contains 825 sample IDs in total 

writeLines(unlist(mdat_sel), paste0(out_dir,"guppy_input/piggies_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


# ALL PIGGIES ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


# SUBSET PIGGIES (two breeds, 3 days bday diff within each breed); groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-scour"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# merge breed and bday details to metadata
mdat_sel <- inner_join(mdat_sel,details)
mdat_sel$BIRTH_DAY <- as.character(mdat_sel$BIRTH_DAY)

##################
##################

# check distribution of breeds and bdays across cohorts to choose a breed and bdays to form groups 

distribution_check <- mdat_sel %>%
  dplyr::select(breed,BIRTH_DAY,Cohort,collection_date)

# save distribution table 
distribution_check$grouping <- paste0(distribution_check$breed,"_",distribution_check$BIRTH_DAY)
x <- ctable(x = distribution_check$Cohort, y = distribution_check$grouping, prop = "r")

view(x, file = paste0(out_dir,"distribution_check.html"))

# print distribution table on console 
distribution_check <- xtabs(~ breed+BIRTH_DAY+Cohort, data=distribution_check)
ftable(distribution_check)

##################
##################

# decision : 
# group_A: breed: Duroc x Landrace, BIRTH_DAY: 2017-01-09 , 2017-01-10, 2017-01-11
# group_B: breed: Duroc x Large white, BIRTH_DAY: 2017-01-09 , 2017-01-10, 2017-01-11

##################
##################

# prepare metadata selections for groups A and B: 

group_A <- mdat_sel %>%
  dplyr::filter(breed=="Duroc x Landrace") %>%
  dplyr::filter(BIRTH_DAY=="2017-01-09"|BIRTH_DAY=="2017-01-10"|BIRTH_DAY=="2017-01-11")

group_B <- mdat_sel %>%
  dplyr::filter(breed=="Duroc x Large white") %>%
  dplyr::filter(BIRTH_DAY=="2017-01-09"|BIRTH_DAY=="2017-01-10"|BIRTH_DAY=="2017-01-11")


# dates selection for group A: 

group_A_t0 <- group_A %>% 
  dplyr::filter(collection_date == "t0")

group_A_t2 <- group_A %>% 
  dplyr::filter(collection_date == "t2")

group_A_t4 <- group_A %>% 
  dplyr::filter(collection_date == "t4")

group_A_t6 <- group_A %>% 
  dplyr::filter(collection_date == "t6")

group_A_t8 <- group_A %>% 
  dplyr::filter(collection_date == "t8")

group_A_t9 <- group_A %>% 
  dplyr::filter(collection_date == "t9")

A_a <- as.character(group_A_t0$ID)
A_b <- as.character(group_A_t2$ID)
A_c <- as.character(group_A_t4$ID)
A_d <- as.character(group_A_t6$ID)
A_e <- as.character(group_A_t8$ID)
A_f <- as.character(group_A_t9$ID)

writeLines(unlist(A_a), paste0(out_dir,"guppy_input/piggies_group_A_t0_sel.txt"), sep = " ")
writeLines(unlist(A_b), paste0(out_dir,"guppy_input/piggies_group_A_t2_sel.txt"), sep = " ")
writeLines(unlist(A_c), paste0(out_dir,"guppy_input/piggies_group_A_t4_sel.txt"), sep = " ")
writeLines(unlist(A_d), paste0(out_dir,"guppy_input/piggies_group_A_t6_sel.txt"), sep = " ")
writeLines(unlist(A_e), paste0(out_dir,"guppy_input/piggies_group_A_t8_sel.txt"), sep = " ")
writeLines(unlist(A_f), paste0(out_dir,"guppy_input/piggies_group_A_t9_sel.txt"), sep = " ")

# dates selection for group A: 

group_B_t0 <- group_B %>% 
  dplyr::filter(collection_date == "t0")

group_B_t2 <- group_B %>% 
  dplyr::filter(collection_date == "t2")

group_B_t4 <- group_B %>% 
  dplyr::filter(collection_date == "t4")

group_B_t6 <- group_B %>% 
  dplyr::filter(collection_date == "t6")

group_B_t8 <- group_B %>% 
  dplyr::filter(collection_date == "t8")

group_B_t9 <- group_B %>% 
  dplyr::filter(collection_date == "t9")


B_a <- as.character(group_B_t0$ID)
B_b <- as.character(group_B_t2$ID)
B_c <- as.character(group_B_t4$ID)
B_d <- as.character(group_B_t6$ID)
B_e <- as.character(group_B_t8$ID)
B_f <- as.character(group_B_t9$ID)


writeLines(unlist(B_a), paste0(out_dir,"guppy_input/piggies_group_B_t0_sel.txt"), sep = " ")
writeLines(unlist(B_b), paste0(out_dir,"guppy_input/piggies_group_B_t2_sel.txt"), sep = " ")
writeLines(unlist(B_c), paste0(out_dir,"guppy_input/piggies_group_B_t4_sel.txt"), sep = " ")
writeLines(unlist(B_d), paste0(out_dir,"guppy_input/piggies_group_B_t6_sel.txt"), sep = " ")
writeLines(unlist(B_e), paste0(out_dir,"guppy_input/piggies_group_B_t8_sel.txt"), sep = " ")
writeLines(unlist(B_f), paste0(out_dir,"guppy_input/piggies_group_B_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and NEO ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="Neomycin") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+D ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_NEONEOD_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_NEONEOD_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_NEONEOD_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_NEONEOD_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_NEONEOD_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_NEONEOD_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+C ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_NEONEOC_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_NEONEOC_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_NEONEOC_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_NEONEOC_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_NEONEOC_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_NEONEOC_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Ds ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLDs_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLDs_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLDs_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLDs_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLDs_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLDs_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Co ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLC_t0_sel.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLC_t2_sel.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLC_t4_sel.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLC_t6_sel.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLC_t8_sel.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLC_t9_sel.txt"), sep = " ")


###########################################################################################
###########################################################################################


# these are going to guppy: 

# 1 piggies_sel.txt
# 1 pos_controls_sel.txt
# 6 piggies_*.txt
# 6 piggies_group_A_*.txt
# 6 piggies_group_B_*.txt
# 6 piggies_CTRLNEO_*.txt
# 6 piggies_NEONEOD_*.txt
# 6 piggies_NEONEOC_*.txt
# 6 piggies_CTRLDs_*.txt
# 6 piggies_CTRLC_*.txt

###########################################################################################
###########################################################################################


# NEXT: 

# move all the *_sel.txt files on to HPC here /shared/homes/pig_microbiome/phy_10M_202109/guppy_groups/

# necessary input files in the sample directory: 
# all the .jplace files are 960 :
# (base) phylosift_metapigs_20200225/$ ls *.jplace | wc -l
# 960


# script to run: 

# nano run_guppy.sh

# ```
# #!/bin/bash
# #PBS -l ncpus=10
# #PBS -l walltime=120:00:00
# #PBS -l mem=70g
# #PBS -N run_guppy.sh
# #PBS -M daniela.gaio@uts.edu.au
# 
# cd /shared/homes/s1/pig_microbiome/phy_10M/PS_temp
# 
# # compute per-sample alpha diversity with various diversity metrics
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fpd plate_7_*1.fastq.gz/* > all.alphadiv
# 
# # cluster the samples: performs squash clustering - NOT RUN
# # /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy squash plate_*.fastq.gz/*jplace
# 
# # edge PCA to explore variation in community composition among samples: performs edge principal components
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_all plate_*.fastq.gz/*jplace
# 
# # run guppy fat: makes trees with edges fattened in proportion to the number of reads
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_ plate_*.fastq.gz/*jplace
# ```
# works! 



# output files are: (extensions .xml, .edgediff, .trans, .jplace)
# move all the *_sel.txt.jplace and all the *_sel.txt.xml files to local machine ~/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output
# .xml and .jplace files will be processed  with guppy_XML_process.R 

###########################################################################################
###########################################################################################


# ADDITIONAL ANALYSIS: COHORT/TREATMENT GROUPS EXCLUDING PIGS > 3 DAYS BDAY APART

details <- details %>%
  dplyr::select(isolation_source,BIRTH_DAY)
unique(details$BIRTH_DAY)

keep <- details %>%
  dplyr::filter(BIRTH_DAY > "2017-01-08 00:00:00") 
unique(keep$BIRTH_DAY)
NROW(keep)
keep %>%
  group_by(BIRTH_DAY) %>%
  tally()
keep_ID <- unique(keep$isolation_source)


mdat <- mdat[mdat$isolation_source %in% keep_ID,]

###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and NEO ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="Neomycin") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t0_sel_bdays.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t2_sel_bdays.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t4_sel_bdays.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t6_sel_bdays.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t8_sel_bdays.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLNEO_t9_sel_bdays.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+D ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_NEONEOD_t0_sel_bdays.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_NEONEOD_t2_sel_bdays.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_NEONEOD_t4_sel_bdays.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_NEONEOD_t6_sel_bdays.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_NEONEOD_t8_sel_bdays.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_NEONEOD_t9_sel_bdays.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : NEO and NEO+C ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Neomycin"|Cohort=="Neomycin+ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_NEONEOC_t0_sel_bdays.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_NEONEOC_t2_sel_bdays.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_NEONEOC_t4_sel_bdays.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_NEONEOC_t6_sel_bdays.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_NEONEOC_t8_sel_bdays.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_NEONEOC_t9_sel_bdays.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Ds ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="D-scour") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLDs_t0_sel_bdays.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLDs_t2_sel_bdays.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLDs_t4_sel_bdays.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLDs_t6_sel_bdays.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLDs_t8_sel_bdays.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLDs_t9_sel_bdays.txt"), sep = " ")


###########################################################################################
###########################################################################################


#  PIGGIES : CTRL and Co ; groups by collection date 

# filter to keep piggies samples only 
mdat_sel <- mdat %>% 
  filter(Cohort=="Control"|Cohort=="ColiGuard") %>%
  dplyr::select(isolation_source,Cohort,DNA_plate,DNA_well,collection_date)

mdat_sel$ID <- paste0(mdat_sel$DNA_plate,"_",mdat_sel$DNA_well,"_*")

# groups by date 

mdat_t0 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t0")

mdat_t2 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t2")

mdat_t4 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t4")

mdat_t6 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t6")

mdat_t8 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t8")

mdat_t9 <- mdat_sel %>% 
  dplyr::filter(collection_date == "t9")

a <- as.character(mdat_t0$ID)
b <- as.character(mdat_t2$ID)
c <- as.character(mdat_t4$ID)
d <- as.character(mdat_t6$ID)
e <- as.character(mdat_t8$ID)
f <- as.character(mdat_t9$ID)

writeLines(unlist(a), paste0(out_dir,"guppy_input/piggies_CTRLC_t0_sel_bdays.txt"), sep = " ")
writeLines(unlist(b), paste0(out_dir,"guppy_input/piggies_CTRLC_t2_sel_bdays.txt"), sep = " ")
writeLines(unlist(c), paste0(out_dir,"guppy_input/piggies_CTRLC_t4_sel_bdays.txt"), sep = " ")
writeLines(unlist(d), paste0(out_dir,"guppy_input/piggies_CTRLC_t6_sel_bdays.txt"), sep = " ")
writeLines(unlist(e), paste0(out_dir,"guppy_input/piggies_CTRLC_t8_sel_bdays.txt"), sep = " ")
writeLines(unlist(f), paste0(out_dir,"guppy_input/piggies_CTRLC_t9_sel_bdays.txt"), sep = " ")


###########################################################################################
###########################################################################################


# these are going to guppy: 

# 1 piggies_sel.txt
# 1 pos_controls_sel.txt
# 6 piggies_*.txt
# 6 piggies_group_A_*.txt
# 6 piggies_group_B_*.txt
# 12 piggies_CTRLNEO_*.txt  --> 6 including all, 6 excluding bdays > 3 days apart 
# 12 piggies_NEONEOD_*.txt  --> 6 including all, 6 excluding bdays > 3 days apart 
# 12 piggies_NEONEOC_*.txt  --> 6 including all, 6 excluding bdays > 3 days apart 
# 12 piggies_CTRLDs_*.txt   --> 6 including all, 6 excluding bdays > 3 days apart 
# 12 piggies_CTRLC_*.txt    --> 6 including all, 6 excluding bdays > 3 days apart 

###########################################################################################
###########################################################################################


# NEXT: 

# move all the *_sel.txt files on to HPC here /shared/homes/s1/pig_microbiome/phy_10M_202109/guppy_groups/
# and also to /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/

# necessary input files in the sample directory: 
# all the .jplace files are 960 :
# (base) phylosift_metapigs_20200225/$ ls *.jplace | wc -l
# 960


# script to run: 

# nano run_guppy.sh

# ```
# #!/bin/bash
# #PBS -l ncpus=10
# #PBS -l walltime=120:00:00
# #PBS -l mem=70g
# #PBS -N run_guppy.sh
# #PBS -M daniela.gaio@uts.edu.au
# 
# cd /shared/homes/s1/pig_microbiome/phy_10M/PS_temp
# 
# # compute per-sample alpha diversity with various diversity metrics
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fpd plate_*.fastq.gz/* > all.alphadiv
# 
# # cluster the samples: performs squash clustering - NOT RUN
# # /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy squash plate_*.fastq.gz/*jplace
# 
# # edge PCA to explore variation in community composition among samples: performs edge principal components
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_all plate_*.fastq.gz/*jplace
# 
# # run guppy fat: makes trees with edges fattened in proportion to the number of reads
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_ plate_*.fastq.gz/*jplace
# ```
# works! 

# guppy by group: 
# $ cat run_guppy_bygroup.sh
# #!/bin/bash
# #PBS -l ncpus=20
# #PBS -l walltime=120:00:00
# #PBS -l mem=100g
# #PBS -N run_guppy_bygroup.sh
# #PBS -M daniela.gaio@student.uts.edu.au
# 
# 
# cd /shared/homes/s1/pig_microbiome/phy_10M/PS_temp
# 
# 
# # edge PCA to explore variation in community composition among samples - by group
# for f in /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/*.txt # or: /shared/homes/12705859/phylosift_metapigs_20200225/*.txt
# do N=$(basename $f)
# #cat ../guppy_groups/$N
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_$N `cat ../guppy_groups/$N`
# done
# 
# 
# # run guppy fat - by group
# for f in /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/*.txt # or: /shared/homes/12705859/phylosift_metapigs_20200225/*.txt
# do N=$(basename $f)
# #cat ../guppy_groups/$N
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_$N `cat ../guppy_groups/$N`
# done


###########################################################################################
###########################################################################################


# NEW ANALYSIS BY TREATMENT/COHORT: EXCLUDES PIGS BORN > 3 DAYS APART (PERFORMED ON 20211207)

# cd /shared/homes/s1/pig_microbiome/phy_10M/PS_temp
# 
# 
# # edge PCA to explore variation in community composition among samples - by group
# for f in /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/piggies_*_t*_sel_bdays.txt 
# do N=$(basename $f)
# #cat ../guppy_groups/$N
# /shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_$N `cat ../guppy_groups/$N`
# done
# 

###########################################################################################
###########################################################################################


# output files are: (extensions .xml, .edgediff, .trans, .jplace)
# move all the *_sel.txt.jplace and all the *_sel.txt.xml files to local machine ~/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output
# NEW as of 20211207: move all the *_sel_bdays.txt.jplace and all the *_sel_bdays.txt.xml files to local machine ~/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output
# .xml files will be processed  with guppy_XML_process.R 

###########################################################################################
###########################################################################################

