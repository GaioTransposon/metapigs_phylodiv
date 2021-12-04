
# CONTENT WARNING

# This is the longest chunk of code I ever wrote (also my first script)
# please please do not ask me to refactor it 
# Anticipated apologies to those who will attempt reading it 

######################################################
#                                                    #
#       "Abandon all hope, ye who enter here."       #
#                                                    #
######################################################

#################################
#  R version 3.6.1 (2019-07-05)

# 0   # loading input data
# 1   # batch effect (plot,removal and plot again)
# 1.1 # merge alpha div data with metadata -> boggo
# 1.2 # merge beta div data with metadata -> coggo
# 2   # boggo and coggo: formatting & averaging duplicates 
# 2.1 # plot samples distribution (plate vs time; cohorts vs time)
# 3   # DELTAS
# 4   # plot ALPHA diversity (all timepoints)
# 5   # plot BETA diversity (all timepoints)
# 6   # plot ALPHA diversity (at pig trial start) 
# 7   # plot BETA diversity (at pig trial start) 
# 8   # plot distribution of cross_breeds and bdays among cohorts + cohorts distr plates 
# 9   # p-values
# 10 # plot p-values (starting/individual factors)
# 11  # prepare input files for guppy (select by time interval) and plot output 

##### 0. LOAD_LIBRARIES_fast (not used)
##### 0. LOAD_LIBRARIES ###################################################################

#install.packages("BiocManager")
#BiocManager::install("sva")
library(sva)

#install.packages("tiff")
library(tiff) # non-zero status

#install.packages("rstatix")
library(rstatix)

#install.packages("devtools")
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot) # ggbiplot_0.55 

#install.packages("ggpubr")
library(ggpubr) # ggpubr_0.2.4 

#install.packages("tidyverse")
library(tidyverse) # tidyverse_1.3.0 

library(broom) # broom_0.5.2
library(cowplot) # cowplot_1.0.0
library(data.table) # data.table_1.12.8  

#install.packages("dunn.test")
library(dunn.test) # dunn.test_1.3.5 

library(plyr) # plyr_1.8.5 
library(dplyr) # dplyr_0.8.3 
library(forcats) # forcats_0.4.0 
library(ggplot2) # ggplot2_3.2.1
library(gridExtra) # gridExtra_2.3     

#install.packages("plotrix")
library(plotrix) # plotrix_3.7-7

library(readr) # readr_1.3.1
library(readxl) # readxl_1.3.1
library(tidyr) # tidyr_1.0.0

#install.packages("varhandle")
library(varhandle) # varhandle_2.0.4

library(tibble) # tibble_2.1.3 
library(purrr) # purrr_0.3.3
library(openxlsx)

library(genefilter)

#install.packages("compareGroups")
library(compareGroups)

#install.packages("splitstackshape")
library(splitstackshape)

#install.packages("pheatmap")
library(pheatmap) # used in pos_controls_reads.R

library(magrittr)
#install.packages("EnvStats")
library(EnvStats)

#install.packages("magick")
# install.packages("~/Downloads/magick_2.7.2.tar.gz", repos = NULL, type = "source")
library(magick)



source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/out/" # local 

###### 0. LOAD input data and parse ##################################################################

# 0   # loading input data

# tiffs (timelines)
timeline <- image_read(paste0(out_dir,"timeline.tiff"))

# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$Cohort <- gsub("Mothers","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

mdat <- mdat %>%
  dplyr::select(isolation_source,collection_date,Cohort,DNA_plate,DNA_well,PigPen)

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


# load details (cross_breed, line, bday, Sows)
details <- read_excel(paste0(source_data, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse_sow'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'maternal_sow'
colnames(details)[colnames(details) == '...8'] <- 'cross_breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)

details <- details %>%
  dplyr::select(pig,BIRTH_DAY,LINE,cross_breed,maternal_sow,nurse_sow)
unique(details$cross_breed)

# load alpha pd
fpddat<-read.table(paste0(middle_dir,"all.alphadiv"),header=T,stringsAsFactors=F) # old file: fpdalpha_div.tsv
# load beta div
pcadat<-read_csv(paste0(middle_dir,"pca_piggies_sel.txt.proj"),col_names = FALSE) # olf file: piggies_sel.txt.proj

pcadat <- cSplit(pcadat, "X1","_")

pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
pcadat$DNA_well <- pcadat$X1_3
colnames(pcadat)[1:5] <- c("pc1","pc2","pc3","pc4","pc5")
pcadat <- pcadat %>%
  dplyr::select(DNA_plate,DNA_well,pc1,pc2,pc3,pc4,pc5)

pcadat$DNA_plate <- gsub("plate_","P",pcadat$DNA_plate)


##### 1. BATCH_EFFECT_DETECTION ##################################################################


# 1   # batch effect

# Plots the batch effect (both alpha and beta div)

fpddat <- cSplit(fpddat, "placerun","_")
fpddat$DNA_plate <- paste0(fpddat$placerun_1,"_",fpddat$placerun_2)
fpddat$DNA_plate <- gsub("plate_","P",fpddat$DNA_plate)
fpddat$DNA_well <- fpddat$placerun_3
fpddat <- fpddat %>%
  dplyr::select(DNA_plate,DNA_well, phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd)

fpddat$sid <- NULL

# re-order plates
fpddat$DNA_plate <- factor(fpddat$DNA_plate, 
                           levels=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))
pcadat$DNA_plate <- factor(pcadat$DNA_plate, 
                           levels=c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))

palette <- c("black","red","green3","blue","cyan","magenta","yellow","gray","orange","brown")

pdf(paste0(out_dir,"batch_pca.pdf"))
par(mfrow=c(3,2), mai = c(0.3, 0.3, 0.3, 0.3))
plot(pcadat$pc1,pcadat$pc2,main="PC1 PC2",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc2,pcadat$pc3,main="PC2 PC3",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc3[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc1,pcadat$pc3,main="PC1 PC3",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc3[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc1,pcadat$pc4,main="PC1 PC4",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc1[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc4[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot(pcadat$pc2,pcadat$pc4,main="PC2 PC4",xlab="",ylab="",pch=NA,type="p")
DNA_plates=unique(sort(pcadat$DNA_plate))
for(plate in 1:length(DNA_plates)){
  points(pcadat$pc2[pcadat$DNA_plate==DNA_plates[plate]],pcadat$pc4[pcadat$DNA_plate==DNA_plates[plate]],col=palette,pch=1)
}
plot.new()
legend("center", legend=DNA_plates, title="DNA extraction plate", fill=palette, cex=1.2, ncol=3)
dev.off()


# plot the phylogenetic diversity based on DNA plate (batch)

####### get sample size within each dna plate
cw_summary <- fpddat %>% 
  group_by(DNA_plate) %>% 
  tally()

#font size for pvalues 
your_font_size <- 3

# other fonts
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9))

batch_unroo <- ggboxplot(fpddat, x = "DNA_plate", y = "unrooted_pd",
                         color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 
batch_bw <- ggboxplot(fpddat, x = "DNA_plate", y = "bwpd",
                      color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 

####### get sample size within each dna plate
cw_summary <- pcadat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_pc1 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc1",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-2.5, size = your_font_size) 

batch_pc2 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc2",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=0, size = your_font_size) 

batch_pc3 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc3",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.8, size = your_font_size) 

batch_pc4 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc4",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.5, size = your_font_size) 

batch_pc5 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=4, size = your_font_size)+
  ylim(4,7)


# Extract the legend. Returns a gtable
for_legend_only <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                             color = "DNA_plate", palette = "jco")+
  guides(color=guide_legend(ncol=3)) +
  theme(legend.position=c(0.5,0.5),  
        plot.margin=unit(c(1,1,7,1),"lines")) +
  labs(fill="") 

leg <- get_legend(for_legend_only)


figure <- grid.arrange(
  batch_unroo, batch_bw, batch_pc1, batch_pc2, batch_pc3, batch_pc4, batch_pc5, leg, nrow = 4, ncol = 2
)

pdf(paste0(out_dir,"batch_alpha_beta.pdf"))
annotate_figure(figure,
                top = text_grob("Batch effect by alpha and beta diversity", color = "black", size = 14)
)
dev.off()


# adjusted p-values 

aov.out1 = aov(unrooted_pd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out1)
aov.out1 <- as.data.frame(res$DNA_plate)
aov.out1$type="unrooted_pd"

aov.out2 = aov(bwpd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out2)
aov.out2 <- as.data.frame(res$DNA_plate)
aov.out2$type="bwpd"

aov.out3 = aov(pc1 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out3)
aov.out3 <- as.data.frame(res$DNA_plate)
aov.out3$type="PC1"

aov.out4 = aov(pc2 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out4)
aov.out4 <- as.data.frame(res$DNA_plate)
aov.out4$type="PC2"

aov.out5 = aov(pc3 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out5)
aov.out5 <- as.data.frame(res$DNA_plate)
aov.out5$type="PC3"

aov.out6 = aov(pc4 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out6)
aov.out6 <- as.data.frame(res$DNA_plate)
aov.out6$type="PC4"

aov.out7 = aov(pc5 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out7)
aov.out7 <- as.data.frame(res$DNA_plate)
aov.out7$type="PC5"

aov.out1$type <- "unrooted_pd"
aov.out1$comparison <- rownames(aov.out1)
aov.out2$type <- "bwpd"
aov.out2$comparison <- rownames(aov.out2)
aov.out3$type <- "pc1"
aov.out3$comparison <- rownames(aov.out3)
aov.out4$type <- "pc2"
aov.out4$comparison <- rownames(aov.out4)
aov.out5$type <- "pc3"
aov.out5$comparison <- rownames(aov.out5)
aov.out6$type <- "pc4"
aov.out6$comparison <- rownames(aov.out6)
aov.out7$type <- "pc5"
aov.out7$comparison <- rownames(aov.out7)

all <- rbind(aov.out1,
             aov.out2,
             aov.out3,
             aov.out4,
             aov.out5,
             aov.out6,
             aov.out7)
all$padj_method <- "TukeyHSD"

fwrite(x=all, file=paste0(stats_dir,"batch_pre_process.csv"))

##### 1. BATCH_EFFECT_REMOVAL ###############################################################################

# Removes batch effect: 

PCA_well <- pcadat$DNA_well
PCA_batch <- pcadat$DNA_plate
pcadat<- data.matrix(pcadat[,3:7], rownames.force = NA)
pcadat_clean<-ComBat(dat=t(as.matrix(pcadat)),PCA_batch,mod=NULL)

PD_well <- fpddat$DNA_well
PD_batch <- fpddat$DNA_plate
fpddat<- data.matrix(fpddat[,3:7], rownames.force = NA)
fpddat_clean<-ComBat(dat=t(as.matrix(fpddat)),PD_batch,mod=NULL)

pcadat_clean <- t(pcadat_clean)
fpddat_clean <- t(fpddat_clean)

pcadat_clean <- as.data.frame(pcadat_clean)
fpddat_clean <- as.data.frame(fpddat_clean)

pcadat_clean <- unfactor(pcadat_clean[])
fpddat_clean <- unfactor(fpddat_clean[])

pcadat <- cbind(PCA_batch,PCA_well,pcadat_clean)
fpddat <- cbind(PD_batch,PD_well,fpddat_clean)

# rename cols to dna plate 
colnames(pcadat)[colnames(pcadat)=="PCA_batch"] <- "DNA_plate"
colnames(pcadat)[colnames(pcadat)=="PCA_well"] <- "DNA_well"
colnames(fpddat)[colnames(fpddat)=="PD_batch"] <- "DNA_plate"
colnames(fpddat)[colnames(fpddat)=="PD_well"] <- "DNA_well"


# save new un-batched data 
fwrite(x = pcadat, file = paste0(middle_dir,"pcadat_clean"))
fwrite(x = fpddat, file = paste0(middle_dir,"fpddat_clean"))


######### Plot batch effect after removal of batch effect: 

####### get sample size within each dna plate

#font size for pvalues 
your_font_size <- 3

# other fonts
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9))

####### get sample size within each dna plate
cw_summary <- fpddat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_unroo <- ggboxplot(fpddat, x = "DNA_plate", y = "unrooted_pd",
                         color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 

batch_bw <- ggboxplot(fpddat, x = "DNA_plate", y = "bwpd",
                      color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=1, size = your_font_size) 


cw_summary <- pcadat %>% 
  group_by(DNA_plate) %>% 
  tally()

batch_pc1 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc1",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-2.5, size = your_font_size) 

batch_pc2 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc2",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=0, size = your_font_size) 

batch_pc3 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc3",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.8, size = your_font_size) 

batch_pc4 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc4",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=-1.5, size = your_font_size) 

batch_pc5 <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                       color = "DNA_plate", palette = "jco")+
  My_Theme+
  theme(axis.text.x=element_blank(), legend.position="")+
  geom_text(data = cw_summary,
            aes(DNA_plate, Inf, label = n), vjust="inward", size = your_font_size) +
  stat_compare_means(method = "anova", label.x=2, label.y=4, size = your_font_size)+
  ylim(4,7)


# Extract the legend. Returns a gtable
for_legend_only <- ggboxplot(pcadat, x = "DNA_plate", y = "pc5",
                             color = "DNA_plate", palette = "jco")+
  guides(color=guide_legend(ncol=3)) +
  theme(legend.position=c(0.5,0.5),  
        plot.margin=unit(c(1,1,7,1),"lines")) +
  labs(fill="") 

leg <- get_legend(for_legend_only)

figure <- grid.arrange(
  batch_unroo, batch_bw, batch_pc1, batch_pc2, batch_pc3, batch_pc4, batch_pc5, leg, nrow = 4, ncol = 2
)


# adjusted p-values 

aov.out1 = aov(unrooted_pd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out1)
aov.out1 <- as.data.frame(res$DNA_plate)
aov.out1$type="unrooted_pd"

aov.out2 = aov(bwpd ~ DNA_plate, data=fpddat)
res <- TukeyHSD(aov.out2)
aov.out2 <- as.data.frame(res$DNA_plate)
aov.out2$type="bwpd"

aov.out3 = aov(pc1 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out3)
aov.out3 <- as.data.frame(res$DNA_plate)
aov.out3$type="PC1"

aov.out4 = aov(pc2 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out4)
aov.out4 <- as.data.frame(res$DNA_plate)
aov.out4$type="PC2"

aov.out5 = aov(pc3 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out5)
aov.out5 <- as.data.frame(res$DNA_plate)
aov.out5$type="PC3"

aov.out6 = aov(pc4 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out6)
aov.out6 <- as.data.frame(res$DNA_plate)
aov.out6$type="PC4"

aov.out7 = aov(pc5 ~ DNA_plate, data=pcadat)
res <- TukeyHSD(aov.out7)
aov.out7 <- as.data.frame(res$DNA_plate)
aov.out7$type="PC5"

aov.out1$type <- "unrooted_pd"
aov.out1$comparison <- rownames(aov.out1)
aov.out2$type <- "bwpd"
aov.out2$comparison <- rownames(aov.out2)
aov.out3$type <- "pc1"
aov.out3$comparison <- rownames(aov.out3)
aov.out4$type <- "pc2"
aov.out4$comparison <- rownames(aov.out4)
aov.out5$type <- "pc3"
aov.out5$comparison <- rownames(aov.out5)
aov.out6$type <- "pc4"
aov.out6$comparison <- rownames(aov.out6)
aov.out7$type <- "pc5"
aov.out7$comparison <- rownames(aov.out7)

all <- rbind(aov.out1,
             aov.out2,
             aov.out3,
             aov.out4,
             aov.out5,
             aov.out6,
             aov.out7)
all$padj_method <- "TukeyHSD"

fwrite(x=all, file=paste0(stats_dir,"batch_post_process.csv"))

pdf(paste0(out_dir,"NO_batch_alpha_beta.pdf"))
annotate_figure(figure,
                top = text_grob("Batch effect after batch effect removal", color = "black", size = 14)
)
dev.off()

##### 1. MERGE_ALPHA_BETA ##############################################################################


# 1.1 # merge alpha div data with metadata -> boggo
# 1.2 # merge beta div data with metadata -> coggo

# Merges with metadata and plots alpha diversity

# alpha and beta div df formatting 
fpddat$DNA_plate <- gsub("P","plate_", fpddat$DNA_plate)
pcadat$DNA_plate <- gsub("P","plate_", pcadat$DNA_plate)

# Merge alpha and div to metadata one at a time

# merge metadata with alpha div 
boggo <- inner_join(fpddat,mdat)
NROW(boggo)
boggo$Cohort <- gsub("D-scour","D-Scour", boggo$Cohort)

sum(boggo$Cohort == "PosControl_ColiGuard")
sum(boggo$Cohort == "MockCommunity")
sum(boggo$Cohort == "PosControl_D-Scour")
# NB:
# here cohorts "PosControl_ColiGuard" "PosControl_D-Scour"   "MockCommunity" are present

# merge metadata with beta div 
coggo <- inner_join(pcadat,mdat)
NROW(coggo)

# NB: 
# here cohorts "PosControl_ColiGuard" "PosControl_D-Scour"   "MockCommunity" are missing
# and another 12 samples are missing (911-874=37 which is 8+9+8 pos controls and another 12 samples) 
# These are the non-matching samples: 
nomatch <- anti_join(boggo,coggo)
NROW(nomatch)
unique(nomatch$Cohort)
head(nomatch,37)

###### 2 # boggo and coggo: formatting & averaging duplicates############################################

# use alpha div df (which includes the pos controls) to see how the samples are distributed in the plates 
# how are the samples distributed over the plates ? 
# some cohorts or dates over-represented in a plate? 

# frequency of date by Cohort

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,Cohort)]

df1[order(df1$collection_date)]

p1 <- ggplot(df1, aes(fill=Cohort, y=Freq, x=collection_date)) + 
  geom_bar(position="dodge", stat="identity")+
  labs(x = "sample collection date",
       y = "number of samples",
       fill = "Cohort")+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())

# frequency of DNA_plate by date

df1 <- setDT(boggo)[, .(Freq = .N), by = .(collection_date,DNA_plate)]

df1[order(df1$collection_date)]

p2 <- ggplot(df1, aes(fill=DNA_plate, y=Freq, x=collection_date)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=palette)+
  labs(x = "sample collection date",
       y = "number of samples",
       fill = "DNA extraction plate") +
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())


pdf(paste0(out_dir,"distribution_cohorts_DNA_plate_time.pdf"))
ggarrange(
  p1,p2,nrow=2, labels=c("A","B")
)
dev.off()


###

# distribution of subjects across cohorts, birth days, and breeds
mdat_small <- mdat %>%
  dplyr::select(isolation_source,Cohort) %>%
  distinct()
colnames(details) <- c("isolation_source","BIRTH_DAY","LINE","cross_breed","maternal_sow","nurse_sow")
df <- left_join(mdat_small, details)

df1 <- df %>% 
  dplyr::select(isolation_source,Cohort,BIRTH_DAY,cross_breed) 
df2 <- as.data.frame(na.omit(df1)) # removing Sows 

df3 <- df2 %>% group_by(Cohort,cross_breed,BIRTH_DAY) %>% tally()
df3 <- as.data.frame(df3)

pdf(paste0(out_dir,"distribution_subjects.pdf"))
ggplot(df3, aes(x=Cohort,y=n,fill=BIRTH_DAY)) +
  geom_bar(stat="identity")+
  facet_wrap(~cross_breed)+
  theme_bw()+
  labs(y="number of subjects") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7))
dev.off()


####### 2.1 # plot samples distribution across plates (time and cohorts)#################

# Now to get the stats and plot we need to remove (by averaging) 
# samples with identical collection date and identical isolation_source
# these are technical replicates

# aggregating dups for alpha
# select the necessary columns
boggo <- boggo %>%
  dplyr::select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source,
                collection_date, Cohort)
# necessary to remove and add later the pos and neg controls 
# otherwise when aggregating below we would end up with only 
# one replicate per pos/neg control
controls <- boggo %>%
  filter(Cohort == "MockCommunity"|
           Cohort == "PosControl_D-Scour"|
           Cohort == "PosControl_ColiGuard"|
           Cohort == "NegativeControl")

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
boggo <- setDT(boggo)[, lapply(.SD, mean), by=c(names(boggo)[6:8]), .SDcols=cols]
# unite controls and samples 
boggo <- rbind(boggo, controls)
unique(boggo$Cohort)
boggo$Cohort <- factor(boggo$Cohort, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard",
                                "Sows",
                                "MockCommunity",
                                "PosControl_D-Scour",
                                "PosControl_ColiGuard"))

# aggregating dups for beta
# select the necessary columns
coggo <- coggo %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source,
                collection_date, Cohort)

# for beta we don't need to filter out then rbind the pos controls as we don't have them

# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
coggo <- setDT(coggo)[, lapply(.SD, mean), by=c(names(coggo)[6:8]), .SDcols=cols]
coggo$Cohort <- factor(coggo$Cohort, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard",
                                "Sows",
                                "NegativeControl"))





##### 3. DELTAS #############################################################################################


# filtering out piglets that had dysentery
boggo1 <- boggo %>%
  dplyr::filter(!isolation_source=="29665"|isolation_source=="29865"|isolation_source=="29702")

pigs_1 <- boggo1 %>%
  dplyr::filter(collection_date == "t0") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_1)

#hist(pigs_1$unrooted_pd,breaks=100)
#hist(pigs_1$bwpd,breaks=100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_1 <- pigs_1 %>%
  dplyr::filter(!unrooted_pd < 50) %>%
  dplyr::filter(!bwpd >2.6) %>%
  dplyr::filter(!bwpd <1.6)

####

pigs_2 <- boggo1 %>%
  dplyr::filter(collection_date == "t2") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_2)

#hist(pigs_2$unrooted_pd,breaks=100)
#hist(pigs_2$bwpd,breaks=100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_2 <- pigs_2 %>%
  dplyr::filter(!unrooted_pd < 50) %>%
  dplyr::filter(!bwpd >2.5) 

####

pigs_3 <- boggo1 %>%
  dplyr::filter(collection_date == "t4") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_3)

#hist(pigs_3$unrooted_pd)
#hist(pigs_3$bwpd)

####

pigs_4 <- boggo1 %>%
  dplyr::filter(collection_date == "t6") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_4)

#hist(pigs_4$unrooted_pd,breaks=100)
#hist(pigs_4$bwpd,100)

# unrooted_pd has low values indicative of spurious sample. 
# remove those rows
pigs_4 <- pigs_4 %>%
  dplyr::filter(!unrooted_pd < 50) %>%
  dplyr::filter(!bwpd < 1.6)

####

pigs_5 <- boggo1 %>%
  dplyr::filter(collection_date == "t8") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_5)

#hist(pigs_5$unrooted_pd)
#hist(pigs_5$bwpd)

####

pigs_6 <- boggo1 %>%
  dplyr::filter(collection_date == "t9") %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd)
NROW(pigs_6)

#hist(pigs_6$unrooted_pd)
#hist(pigs_6$bwpd)

####

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

theme_4diffs = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 9)) 

####

# Ja31 vs Fe7 - t0 vs t2

df1 <- merge(pigs_1,pigs_2, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  dplyr::arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))



# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t2"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t2"
res2$type <- "bwpd"
A_B <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t2"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t2"
res2$type <- "bwpd"
A_B_adj <- rbind(res1,res2)

####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_B_delta <- both
A_B_delta$time_delta <- "A_B"

####

# Fe7 vs Fe14 - t2 vs t4

df1 <- merge(pigs_2,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t2_vs_t4"
res1$type <- "unrooted_pd"
res2$time_delta <- "t2_vs_t4"
res2$type <- "bwpd"
B_C <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t2_vs_t4"
res1$type <- "unrooted_pd"
res2$time_delta <- "t2_vs_t4"
res2$type <- "bwpd"
B_C_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
B_C_delta <- both
B_C_delta$time_delta <- "B_C"


####


# Fe14 vs Fe21 - t4 vs t6

df1 <- merge(pigs_3,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t6"
res2$type <- "bwpd"
C_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t6"
res2$type <- "bwpd"
C_D_adj <- rbind(res1,res2)

####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_D_delta <- both
C_D_delta$time_delta <- "C_D"

####

# Fe21 vs Fe28 - t6 vs t8

df1 <- merge(pigs_4,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t6 vs t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t6 vs t8"
res2$type <- "bwpd"
D_E <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t6 vs t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t6 vs t8"
res2$type <- "bwpd"
D_E_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
D_E_delta <- both
D_E_delta$time_delta <- "D_E"

####

# Fe28 vs Ma3 - t8 vs t9

df1 <- merge(pigs_5,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t8_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t8_vs_t9"
res2$type <- "bwpd"
E_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t8_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t8_vs_t9"
res2$type <- "bwpd"
E_F_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
E_F_delta <- both
E_F_delta$time_delta <- "E_F"

####

# Ja31 vs Fe14 - t0 vs t4

df1 <- merge(pigs_1,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))

# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t4"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t4"
res2$type <- "bwpd"
A_C <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t4"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t4"
res2$type <- "bwpd"
A_C_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_C_delta <- both
A_C_delta$time_delta <- "A_C"


####

# Fe7 vs Fe21 - t2 vs t6

df1 <- merge(pigs_2,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t2_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t2_vs_t6"
res2$type <- "bwpd"
B_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t2_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t2_vs_t6"
res2$type <- "bwpd"
B_D_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
B_D_delta <- both
B_D_delta$time_delta <- "B_D"

####

# Fe14 vs Fe28 - t4 vs t8

df1 <- merge(pigs_3,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t8"
res2$type <- "bwpd"
C_E <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t8"
res2$type <- "bwpd"
C_E_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_E_delta <- both
C_E_delta$time_delta <- "C_E"


####

# Ja31 vs Fe28 - t0 vs t8

df1 <- merge(pigs_4,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t8"
res2$type <- "bwpd"
D_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t8"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t8"
res2$type <- "bwpd"
D_F_adj <- rbind(res1,res2)


####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
D_F_delta <- both
D_F_delta$time_delta <- "D_F"


####

# Ja31 vs Fe21 - t0 vs t6

df1 <- merge(pigs_1,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t6"
res2$type <- "bwpd"
A_D <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t6"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t6"
res2$type <- "bwpd"
A_D_adj <- rbind(res1,res2)


####


# Fe14 vs Ma3 - t4 vs t9

df1 <- merge(pigs_3,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t9"
res2$type <- "bwpd"
C_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t4_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t4_vs_t9"
res2$type <- "bwpd"
C_F_adj <- rbind(res1,res2)

####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
C_F_delta <- both
C_F_delta$time_delta <- "C_F"


####

# Ja31 vs Ma3 - t0 vs t9

df1 <- merge(pigs_1,pigs_6, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

# does it look right? proof 1: 
test <- setDT(df1)[, .(Freq = .N), by = .(Cohort.x)]
test

# does it look right? proof 2:
test <- setDT(df1)[, .(Freq = .N), by = .(isolation_source)]
test
which(test$Freq!=1)

head(df1)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

NROW(unique(df1$isolation_source))

# reorder
df1$Cohort.x <- factor(df1$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))


# stats without p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "none")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "none")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t9"
res2$type <- "bwpd"
A_F <- rbind(res1,res2)


# stats with p adjustment:

res1 <- pairwise.t.test(df1$diff_unroo, df1$Cohort.x, p.adj = "BH")
res2 <- pairwise.t.test(df1$diff_bw, df1$Cohort.x, p.adj = "BH")
res1 <- res1$p.value
res2 <- res2$p.value
# formatting matrix format (wide) to long format (NA automatically form, omit)
res1 <- na.omit(as.data.frame(matrix(res1, dimnames=list(
  t(outer(colnames(res1), rownames(res1), FUN=paste)), NULL))))
res2 <- na.omit(as.data.frame(matrix(res2, dimnames=list(
  t(outer(colnames(res2), rownames(res2), FUN=paste)), NULL))))
res1$time_delta <- "t0_vs_t9"
res1$type <- "unrooted_pd"
res2$time_delta <- "t0_vs_t9"
res2$type <- "bwpd"
A_F_adj <- rbind(res1,res2)

####

new_df <- merge(df1,details)
new_df$BIRTH_DAY <- as.character(new_df$BIRTH_DAY)
aov.out1 = aov(diff_unroo ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
aov.out2 = aov(diff_bw ~ Cohort.x + cross_breed + BIRTH_DAY, data=new_df)   
res1 <- TukeyHSD(aov.out1)
res2 <- TukeyHSD(aov.out2)
x <- as.data.frame(res1$Cohort)
x$type <- "unrooted_pd"
x <- setDT(x, keep.rownames = TRUE)[]
y <- as.data.frame(res2$Cohort)
y$type <- "bwpd"
y <- setDT(y, keep.rownames = TRUE)[]
both <- rbind(x,y)
both <- both[order(both$rn, both$type), , drop = FALSE]
A_F_delta <- both
A_F_delta$time_delta <- "A_F"

####


# this is for extracting the legend 
for_legend_only <- ggboxplot(df1, x = "Cohort.x", y = "diff_unroo", color = "Cohort.x", 
                             legend = "right")+
  scale_color_manual(labels = c("Control", 
                                "D-Scour",
                                "ColiGuard",
                                "Neo",
                                "Neo+D",
                                "Neo+C"), 
                     values = c("#F8766D", 
                                "#B79F00",
                                "#00BA38",
                                "#00BFC4",
                                "#619CFF",
                                "#F564E3")) +
  guides(color=guide_legend("Cohort")) +
  My_Theme
leg <- get_legend(for_legend_only)



####

# DELTAS p values  WITHOUT p-adjustment

# gather stats of the deltas: 
deltas_p <- rbind(A_B,B_C,C_D,D_E,E_F,A_C,B_D,C_E,D_F,A_D,C_F,A_F)
# convert rownames to first column
deltas_p <- setDT(deltas_p, keep.rownames = TRUE)[]
deltas_p$test <- "pairwise t-test"

####

# DELTAS p values WITH p-adjustment

# gather stats of the deltas: 
deltas_padj <- rbind(A_B_adj,B_C_adj,C_D_adj,D_E_adj,E_F_adj,
                     A_C_adj,B_D_adj,C_E_adj,D_F_adj,
                     A_D_adj,C_F_adj,A_F_adj)
# convert rownames to first column
deltas_padj <- setDT(deltas_padj, keep.rownames = TRUE)[]
deltas_padj$test <- "pairwise t-test"
deltas_padj$method <- "BH"

####

colnames(deltas_p)[colnames(deltas_p) == 'V1'] <- 'pvalues'
colnames(deltas_padj)[colnames(deltas_padj) == 'V1'] <- 'pvalues_adj'

deltas_stats <- merge(deltas_p,deltas_padj, by=c("rn","time_delta","type","test"))
deltas_stats <- cSplit(deltas_stats, "rn"," ")


deltas_stats <- deltas_stats %>%
  # join the two groups for which the p-value has been computed
  dplyr::mutate(comparison=paste0(rn_1,"_vs_",rn_2)) %>%
  dplyr::select(test,type,time_delta,comparison,pvalues,pvalues_adj,method)

# remove useless digits in strings "comparison"
deltas_stats$comparison <- deltas_stats$comparison %<>%
  gsub('[0-9]+', '', .)  # removes digits

# filtering to keep only meaningful comparisons 
# to be kept: 
meaningfulcomparisons <- c("Control_vs_ColiGuard", "ColiGuard_vs_Control",
                           "Control_vs_D-Scour", "D-Scour_vs_Control",
                           "Control_vs_Neomycin", "Neomycin_vs_Control",
                           "Neomycin_vs_Neomycin+D-Scour", "Neomycin+D-Scour_vs_Neomycin",
                           "Neomycin_vs_Neomycin+ColiGuard", "Neomycin+ColiGuard_vs_Neomycin")

# eliminate useless comparisons
deltas_stats <- deltas_stats[deltas_stats$comparison %in% meaningfulcomparisons,]


# add data to workbook 
# addWorksheet(wb, "deltas_stats")
# writeData(wb, sheet = "deltas_stats", deltas_stats, rowNames = FALSE)


####
####


# DELTAS p values WITH Tukey p-adjustment

# gather stats of the deltas: 
deltas_padj_w_model <- rbind(A_B_delta,B_C_delta,C_D_delta,D_E_delta,E_F_delta,
                             A_C_delta,B_D_delta,C_E_delta,D_F_delta,
                             C_F_delta,A_F_delta)

deltas_padj_w_model$test <- "pairwise t-test"
deltas_padj_w_model$padj_method <- "Tukey"
deltas_padj_w_model$model <- "value~Cohort+cross_breed+BIRTH_DAY"
# add data to workbook 
# addWorksheet(wb, "deltas_padj_w_model")
# writeData(wb, sheet = "deltas_padj_w_model", deltas_padj_w_model, rowNames = FALSE)

####


# DELTAS - plots to visualize the proportion 
# of increase-decrease of alpha diversity by cohort
# within time frame 


####

# settings for plots: 

#font size for pvalues 
your_font_size <- 2 

My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

theme_4diffs = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 9),
  axis.ticks.x = element_blank())

####

# Ja31 vs Fe7 - t0 vs t2

df1 <- merge(pigs_1,pigs_2, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_a_b <- df1
df_a_b$interval <- "t0-t2"


####

# Fe7 vs Fe14 - t2 vs t4

df1 <- merge(pigs_2,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_b_c <- df1
df_b_c$interval <- "t2-t4"


####

# Fe14 vs Fe21 - t4 vs t6

df1 <- merge(pigs_3,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_c_d <- df1
df_c_d$interval <- "t4-t6"

####

# Fe21 vs Fe28 - t6 vs t8

df1 <- merge(pigs_4,pigs_5, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_d_e <- df1
df_d_e$interval <- "t6-t8"

####

# Ja31 vs Fe14 - t0 vs t4

df1 <- merge(pigs_1,pigs_3, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_a_c <- df1
df_a_c$interval <- "t0-t4"

####


# Fe7 vs Fe21 - t2 vs t6

df1 <- merge(pigs_2,pigs_4, by=c("isolation_source"))
NROW(df1)

# pivot long
df1 <- df1 %>%
  dplyr::select(isolation_source,Cohort.x,unrooted_pd.x,unrooted_pd.y,bwpd.x,bwpd.y) %>% 
  group_by(isolation_source) %>% slice(1) %>%
  arrange(Cohort.x, isolation_source)

df1$diff_unroo = ((df1$unrooted_pd.y-df1$unrooted_pd.x)/df1$unrooted_pd.y)*100
df1$diff_bw = ((df1$bwpd.y-df1$bwpd.x)/df1$bwpd.y)*100

df1 <- df1 %>%
  dplyr::select(isolation_source,diff_unroo,diff_bw,Cohort.x)
df_b_d <- df1
df_b_d$interval <- "t2-t6"


####

all <- rbind(df_a_b,df_b_c, df_c_d, df_d_e,
             df_a_c, df_b_d)


all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t0-t2", 
  replacement = "A-B", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t2-t4", 
  replacement = "B-C", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t4-t6", 
  replacement = "C-D", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t6-t8", 
  replacement = "D-E", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t0-t4", 
  replacement = "A-C", 
  fixed = TRUE)
all[5] <- lapply(
  all[5], 
  gsub, 
  pattern = "t2-t6", 
  replacement = "B-D", 
  fixed = TRUE)

all$interval <- factor(all$interval,levels=c("A-B",
                                             "B-C",
                                             "C-D",
                                             "D-E",
                                             "A-C",
                                             "B-D"))


####


# Group: Control, D-Scour, ColiGuard, Neo

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  dplyr::mutate(.r = row_number()) 

palette <- c("#F8766D","#B79F00","#00BA38","#00BFC4")

new_df = new_df %>% 
  dplyr::filter(Cohort.x=="Control"|Cohort.x=="D-Scour"|Cohort.x=="ColiGuard"|Cohort.x=="Neomycin")

CTRL_unroo_deltas_facets <- ggplot(data = new_df,
                                   mapping = aes(x = .r, y = diff_unroo, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-100,50)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("unrooted PD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette)

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_bw) %>%
  dplyr::mutate(.r = row_number()) # Add a row number variable

new_df = new_df %>% 
  dplyr::filter(Cohort.x=="Control"|Cohort.x=="D-Scour"|Cohort.x=="ColiGuard"|Cohort.x=="Neomycin")

CTRL_bw_deltas_facets <- ggplot(data = new_df,
                                mapping = aes(x = .r, y = diff_bw, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-30,20)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("BWPD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette)

####

# Group: Neo, Neo+D-Scour, Neo+ColiGuard

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  dplyr::mutate(.r = row_number()) 

palette_neo <- c("#00BFC4","#619CFF","#F564E3")

new_df = new_df %>% 
  dplyr::filter(Cohort.x=="Neomycin"|Cohort.x=="Neomycin+D-Scour"|Cohort.x=="Neomycin+ColiGuard")

NEO_unroo_deltas_facets <- ggplot(data = new_df,
                                  mapping = aes(x = .r, y = diff_unroo, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-100,50)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("unrooted PD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette_neo)

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_bw) %>%
  dplyr::mutate(.r = row_number()) # Add a row number variable

new_df = new_df %>% 
  dplyr::filter(Cohort.x=="Neomycin"|Cohort.x=="Neomycin+D-Scour"|Cohort.x=="Neomycin+ColiGuard")

NEO_bw_deltas_facets <- ggplot(data = new_df,
                               mapping = aes(x = .r, y = diff_bw, fill = Cohort.x))+
  geom_col(width=1) + 
  ylim(-30,20)+
  facet_grid(~interval, 
             scales = "free")+
  theme_bw()+
  theme_4diffs +
  ylab("BWPD - change (%)") +
  guides(fill = FALSE)+
  scale_fill_manual(values=palette_neo)

####

both_bw <- plot_grid(CTRL_bw_deltas_facets,
                     NEO_bw_deltas_facets,
                     nrow=2,
                     rel_heights=c(0.5,0.5))

empty_space <- plot_grid(NULL, NULL, NULL, NULL,
                         ncol=4)

all_plots <- plot_grid(empty_space,
                       both_bw,
                       nrow=2,
                       rel_heights=c(0.25,0.5))

pdf(paste0(out_dir,"cohorts_deltas_bwpd.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(all_plots)
dev.off()


####


both_unroo <- plot_grid(CTRL_unroo_deltas_facets,
                        NEO_unroo_deltas_facets,
                        nrow=2,
                        rel_heights=c(0.5,0.5))

all_plots <- plot_grid(empty_space,
                       both_unroo,
                       nrow=2,
                       rel_heights=c(0.25,0.5))

pdf(paste0(out_dir,"cohorts_deltas_unroo.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(all_plots)
dev.off()

####

##### percentage increase/decrease ##################################################################

# means, standard deviations, percentage incr and decr: 

new_df <- all
new_df = new_df %>% 
  ungroup() %>%
  arrange(interval,Cohort.x,diff_unroo) %>%
  dplyr::mutate(.r = row_number()) 


numbers_bw <- new_df %>% 
  group_by(interval,Cohort.x) %>% 
  dplyr::summarize(mean = mean(diff_bw),
                   sum = sum(diff_bw),
                   sd = sd(diff_bw),
                   n = n(),
                   n_pos = sum(diff_bw > 0),
                   perc_pos = n_pos / n,
                   n_neg = sum(diff_bw < 0),
                   perc_neg = n_neg / n)

numbers_unroo <- new_df %>% 
  group_by(interval,Cohort.x) %>% 
  dplyr::summarize(mean = mean(diff_unroo),
                   sum = sum(diff_unroo),
                   sd = sd(diff_unroo),
                   n = n(),
                   n_pos = sum(diff_unroo > 0),
                   perc_pos = n_pos / n,
                   n_neg = sum(diff_unroo < 0),
                   perc_neg = n_neg / n)

numbers_all_bw <- new_df %>% 
  group_by(interval) %>% 
  dplyr::summarize(mean = mean(diff_bw),
                   sum = sum(diff_bw),
                   sd = sd(diff_bw),
                   n = n(),
                   n_pos = sum(diff_bw > 0),
                   perc_pos = n_pos / n,
                   n_neg = sum(diff_bw < 0),
                   perc_neg = n_neg / n)

numbers_all_unroo <- new_df %>% 
  group_by(interval) %>% 
  dplyr::summarize(mean = mean(diff_unroo),
                   sum = sum(diff_unroo),
                   sd = sd(diff_unroo),
                   n = n(),
                   n_pos = sum(diff_unroo > 0),
                   perc_pos = n_pos / n,
                   n_neg = sum(diff_unroo < 0),
                   perc_neg = n_neg / n)


numbers_unroo$type = "unrooted_pd"
numbers_bw$type = "bwpd"
numbers_all_unroo$type = "unrooted_pd"
numbers_all_bw$type = "bwpd"
numbers_all_unroo$type = "unrooted_pd"
numbers_all_unroo$Cohort.x = "all"
numbers_all_bw$Cohort.x = "all"

numbers_unroo <- as.data.frame(numbers_unroo)
numbers_all_unroo <- as.data.frame(numbers_all_unroo)
numbers_bw <- as.data.frame(numbers_bw)
numbers_all_bw <- as.data.frame(numbers_all_bw)


# into a single dataframe : 

numbers <- rbind(numbers_all_unroo,numbers_unroo,numbers_all_bw,numbers_bw)

numbers$dates <- numbers$interval

numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t0-t2", 
  pattern = "A-B", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t2-t4", 
  pattern = "B-C", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t4-t6", 
  pattern = "C-D", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t6-t8", 
  pattern = "D-E", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t0-t4", 
  pattern = "A-C", 
  fixed = TRUE)
numbers[12] <- lapply(
  numbers[12], 
  gsub, 
  replacement = "t2-t6", 
  pattern = "B-D", 
  fixed = TRUE)


deltas_percent_change <- numbers 

# add data to workbook 
# addWorksheet(wb, "deltas_percent_change")
# writeData(wb, sheet = "deltas_percent_change", deltas_percent_change, rowNames = FALSE)

fwrite(x=deltas_percent_change, file=paste0(stats_dir,"deltas_percent_change.csv"))


######PLOT_ALPHA####################################################################

# 3   # plot ALPHA diversity (all timepoints)

# ALPHA diversity overall (includes pos controls):

boggo1 <- boggo %>%
  dplyr::filter(!isolation_source == "NegativeControl")

# reordering
boggo1$Cohort <- factor(boggo1$Cohort, 
                        levels=c("Control", 
                                 "D-Scour", 
                                 "ColiGuard",
                                 "Neomycin",
                                 "Neomycin+D-Scour",
                                 "Neomycin+ColiGuard",
                                 "Sows",
                                 "MockCommunity",
                                 "PosControl_D-Scour",
                                 "PosControl_ColiGuard"))

boggo2 <- boggo1

boggo2 <- na.omit(boggo2)

# boxplots again for alpha diversity, to be plotted in the same pdf
p1 <- ggplot(boggo2, aes(x=Cohort, y=phylo_entropy)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p2 <- ggplot(boggo2, aes(x=Cohort, y=bwpd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()
p3 <- ggplot(boggo2, aes(x=Cohort, y=unrooted_pd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=1, notch=FALSE) +
  xlab(NULL) +
  coord_flip()


pdf(paste0(out_dir,"alpha_BWPD_unrooted.pdf"))
plot_grid(p2,p3,nrow=2, labels = c("A","B"))
dev.off()

pdf(paste0(out_dir,"alpha_all.pdf"))
plot_grid(p1,p2,p3,nrow=3)
dev.off()

#######boggo means##############################




a <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-Scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Sows") %>% 
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,n = n()
                   ,sd = sd(unrooted_pd)
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75)) 

b <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-Scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Control") %>% 
  filter(!Cohort=="D-Scour") %>% 
  filter(!Cohort=="ColiGuard") %>% 
  filter(!Cohort=="Neomycin") %>% 
  filter(!Cohort=="Neomycin+D-Scour") %>% 
  filter(!Cohort=="Neomycin+ColiGuard") %>% 
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,n = n()
                   ,sd = sd(unrooted_pd)
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75)) 

c <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-Scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Sows") %>% 
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,n = n()
                   ,sd = sd(bwpd)
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

d <- boggo %>% 
  filter(!Cohort=="MockCommunity") %>% 
  filter(!Cohort=="PosControl_D-Scour") %>% 
  filter(!Cohort=="PosControl_ColiGuard") %>% 
  filter(!Cohort=="NegativeControl") %>% 
  filter(!Cohort=="Control") %>% 
  filter(!Cohort=="D-Scour") %>% 
  filter(!Cohort=="ColiGuard") %>% 
  filter(!Cohort=="Neomycin") %>% 
  filter(!Cohort=="Neomycin+D-Scour") %>% 
  filter(!Cohort=="Neomycin+ColiGuard") %>% 
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,n = n()
                   ,sd = sd(bwpd)
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 


e <- boggo %>% 
  na.omit(.) %>%
  group_by(Cohort) %>% 
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,n = n()
                   ,sd = sd(unrooted_pd)
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75))

f <- boggo %>% 
  na.omit(.) %>%
  group_by(Cohort) %>% 
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,n = n()
                   ,sd = sd(bwpd)
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

a$Cohort="piglets"
b$Cohort="Sows"
c$Cohort="piglets"
d$Cohort="Sows"

a$type="unrooted"
b$type="unrooted"
c$type="bwpd"
d$type="bwpd"
e$type="unrooted"
f$type="bwpd"

means <- rbind(a,b,c,d,e,f)
means$collection_date = "all"
# these are added later below into a df 


#######time_alpha############################

# alpha_timeseries_all_cohorts multiple time points (6)

doggo <- boggo
doggo$collection_date <- as.character(doggo$collection_date)

doggo <- doggo %>% dplyr::filter(collection_date == "t0" |
                            collection_date == "t2" |
                            collection_date == "t4" |
                            collection_date == "t6" |
                            collection_date == "t8" |
                            collection_date == "t9") 

# filter out heavy outliers

doggo <- doggo %>% 
  dplyr::filter(!unrooted_pd < 25)

doggo <- doggo %>% 
  dplyr::filter(!bwpd > 2.7)

# check
#hist(doggo$unrooted_pd)
#hist(doggo$bwpd)

# (no need to filter out pos, neg controls and Sows 
# as filtering by date is already doing it )

# reordering
doggo$Cohort <- factor(doggo$Cohort, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))

# reordering
doggo$collection_date <- factor(doggo$collection_date,
                                levels=c("t0" ,
                                         "t2",
                                         "t4",
                                         "t6",
                                         "t8",
                                         "t9"))

my_comparisons = list( c("t0", "t2"), 
                       c("t2", "t4"), 
                       c("t4", "t6"),
                       c("t2", "t6"),
                       c("t6", "t8"),
                       c("t4", "t8"),
                       c("t8", "t9"))

# general time change - unrooted
summs_unroo <- doggo %>% group_by(collection_date,Cohort) %>% 
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,n = n()
                   ,sd = sd(unrooted_pd)
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75)) 

gen_unrooted <- ggplot(summs_unroo, aes(x=collection_date, y=mean, group=Cohort, color=Cohort)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, size = 2)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "unrooted PD - mean")


# general time change - unrooted
summs_bw <- doggo %>% group_by(collection_date,Cohort) %>% 
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,n = n()
                   ,sd = sd(bwpd)
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

gen_bwpd <- ggplot(summs_bw, aes(x=collection_date, y=mean, group=Cohort, color=Cohort)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(size=0.8)+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, size=2)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "collection date",
       y = "BWPD - mean")



tosave <- ggarrange(gen_unrooted, gen_bwpd, ncol = 2, labels=c("A","B"), common.legend = TRUE)

pdf(paste0(out_dir,"time_alpha.pdf"),width = 9,height=6)
tosave
dev.off()




#######alpha_means########################################################

# stats 

unroo <- doggo %>%
  group_by(collection_date,Cohort) %>%
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,sd = sd(unrooted_pd)
                   ,n = n()
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75)) 

unroo_all <- doggo %>%
  group_by(collection_date) %>%
  dplyr::summarise(min = min(unrooted_pd)
                   ,max = max(unrooted_pd)
                   ,mean = mean(unrooted_pd)
                   ,sd = sd(unrooted_pd)
                   ,n = n()
                   ,q25 = quantile(unrooted_pd, .25)
                   ,q75 = quantile(unrooted_pd, .75)) 

bw <- doggo %>%
  group_by(collection_date,Cohort) %>%
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,sd = sd(bwpd)
                   ,n = n()
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

bw_all <- doggo %>%
  group_by(collection_date) %>%
  dplyr::summarise(min = min(bwpd)
                   ,max = max(bwpd)
                   ,mean = mean(bwpd)
                   ,sd = sd(bwpd)
                   ,n = n()
                   ,q25 = quantile(bwpd, .25)
                   ,q75 = quantile(bwpd, .75)) 

unroo_all$Cohort <- "all"
bw_all$Cohort <- "all"

# add data to workbook 

bw$type="bwpd"
bw_all$type="bwpd"
unroo$type="unrooted_pd"
unroo_all$type="unrooted_pd"
bw <- as.data.frame(bw)
bw_all <- as.data.frame(bw_all)
unroo <- as.data.frame(unroo)
unroo_all <- as.data.frame(unroo_all)
both <- rbind(bw,unroo,bw_all,unroo_all,means)
# addWorksheet(wb, "alpha_means")
# writeData(wb, sheet = "alpha_means", both, rowNames = FALSE)
fwrite(x=both, file=paste0(stats_dir,"alpha_means.csv"))


#######time_alpha_cohorts_unroo_BW####################################################################################


# alpha diversity change during time in cohorts, with unrooted PD & BWPD in the same plot:

x <- doggo %>%
  dplyr::select(isolation_source,Cohort,collection_date,unrooted_pd,bwpd) %>%
  droplevels()

summs <- x %>% group_by(Cohort,collection_date) %>% 
  dplyr::summarise(`mean BWPD` = mean(bwpd),
                   sdUnroo = sd(unrooted_pd),
                   meanUnroo = mean(unrooted_pd)) %>%
  droplevels()

summs <- as.data.frame(summs)

pdf(paste0(out_dir,"time_alpha_cohorts_unroo_BW.pdf"))
ggplot(summs, aes(x=collection_date, y=meanUnroo,group=Cohort,color=Cohort)) +
  geom_errorbar(aes(ymin=meanUnroo-sdUnroo, ymax=meanUnroo+sdUnroo), 
                width=0.1,
                size=0.2,
                position=position_dodge(0.5)) +
  geom_line() + geom_point(aes(size=`mean BWPD`)) +
  facet_wrap(~Cohort)+
  labs(y="mean unrooted PD")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
dev.off()


#####################################################

# powerpoint "Desktop/metapgs_base/phylosift/PD_change_cohorts_time.pptx"
# 
# xx <- summs %>%
#   dplyr::mutate(rankedUnroo = rank(meanUnroo))
# 
# # number of groups large enough to appreciate diffs between timepoints 
# xx$groupUnroo <- as.numeric(cut_number(xx$rankedUnroo,4))
# 
# Relabel = c(4,6,8,10)
# xx$groupUnroo = Relabel[xx$groupUnroo]
# 
# # BWPD into 4 categories
# 
# # number of groups large enough to appreciate diffs between timepoints 
# xx$groupBW <- as.numeric(cut_number(xx$meanBW,4))
# 
# Relabel = c(0.25,0.50,0.75,1.00)
# xx$groupBW = Relabel[xx$groupBW]
# xx

#######time_alpha_cohorts########################################################################


# alpha_timeseries_all_cohorts 4 time points 

doggo <- boggo
doggo$collection_date <- as.character(doggo$collection_date)

doggo <- doggo %>% dplyr::filter(collection_date == "t0" |
                            collection_date == "t2" |
                            collection_date == "t4" |
                            collection_date == "t6" ) 


doggo[,3] <- lapply(
  doggo[,3], 
  gsub, 
  pattern = "Neomycin+D-Scour", 
  replacement = "Neo+D-Scour", 
  fixed = TRUE)

doggo[,3] <- lapply(
  doggo[,3], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "Neo+ColiGuard", 
  fixed = TRUE)

# filter out heavy outliers

doggo <- doggo %>% 
  dplyr::filter(!unrooted_pd < 25)

doggo <- doggo %>% 
  dplyr::filter(!bwpd > 2.7)

# check
#hist(doggo$unrooted_pd)
#hist(doggo$bwpd)

# (no need to filter out pos, neg controls and Sows 
# as filtering by date is already doing it )

# reordering
doggo$Cohort <- factor(doggo$Cohort, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neo+D-Scour",
                                "Neo+ColiGuard"))

# reordering
doggo$collection_date <- factor(doggo$collection_date,
                                levels=c("t0",
                                         "t2",
                                         "t4",
                                         "t6"))

my_comparisons = list( c("t0", "t2"), 
                       c("t2", "t4"), 
                       c("t4", "t6"))


stat.test_unroo <- doggo %>%
  group_by(Cohort) %>%
  t_test(unrooted_pd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") %>%
  dplyr::mutate(y.position=rep(seq(370,480,length.out=6),6)) %>%
  dplyr::mutate_if(is.numeric, round, digits = 4)


unroo <- ggboxplot(doggo, x = "collection_date", y = "unrooted_pd",
                   color = "collection_date", palette = "jco",
                   add = "jitter",short.panel.labs = TRUE) +
  theme_bw()+
  facet_grid(~Cohort)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  ylim(200,500)+
  stat_pvalue_manual(stat.test_unroo, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 2)

###

stat.test_bwpd <- doggo %>%
  group_by(Cohort) %>%
  t_test(bwpd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") %>%
  mutate(y.position=rep(seq(2.3,2.65,length.out=6),6)) %>%
  mutate_if(is.numeric, round, digits = 4)

bw <- ggboxplot(doggo, x = "collection_date", y = "bwpd",
                color = "collection_date", palette = "jco",
                add = "jitter",short.panel.labs = TRUE) +
  facet_grid(~Cohort)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  ylim(1.5,2.70)+
  stat_pvalue_manual(stat.test_bwpd, label = "p.adj",
                     hide.ns=TRUE,
                     bracket.size = 0.3,
                     size = 2)


tosave <- ggarrange(unroo,bw,nrow=2,ncol=1,labels=c("A","B"),common.legend = TRUE)

pdf(paste0(out_dir,"time_alpha_cohorts.pdf"))
tosave
dev.off()

# add data to workbook 
both <- rbind(stat.test_bwpd,
              stat.test_unroo)
both$padj_method <- "bonferroni"

# addWorksheet(wb, "alpha_cohorts")
# writeData(wb, sheet = "alpha_cohorts", both, rowNames = FALSE)


#######alpha_time################

# comparing timepoints, irrespective of cohort

stat.test_unroo_all <- doggo %>%
  t_test(unrooted_pd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") 

stat.test_bwpd_all <- doggo %>%
  t_test(bwpd ~ collection_date) %>%
  adjust_pvalue(method="bonferroni") 

# add data to workbook 
both <- rbind(stat.test_bwpd_all,
              stat.test_unroo_all)
both$padj_method <- "bonferroni"

# addWorksheet(wb, "alpha_time")
# writeData(wb, sheet = "alpha_time", both, rowNames = FALSE)
fwrite(x=both, file=paste0(stats_dir,"alpha_time.csv"))

#######nursesow_alpha#####################################################################################

# 5   # plot ALPHA diversity (at pig trial start) 

# ALPHA diversity in piglets at the start of the trial 

startDF <- boggo 

startDF <- startDF %>% filter(
  collection_date == "t0")  %>%
  dplyr::select(phylo_entropy,quadratic,unrooted_pd,rooted_pd,bwpd,isolation_source)

# as we have 160 samples for 126 piglets at the startof the trial. this is because we have duplicates
# samples have been taken from the same animals twice
length(startDF$isolation_source)
length(unique(startDF$isolation_source))
# we need to average 
head(startDF)
# aggregate by taking the mean when subject and collection date is identical 
cols <- 1:5
startDF <- setDT(startDF)[, lapply(.SD, mean), by=c(names(startDF)[6]), .SDcols=cols]
# now we have the right number: 
length(startDF$isolation_source)
length(unique(startDF$isolation_source))


####

startDF1 <- merge(startDF,details, by="isolation_source")

startDF1 <- startDF1 %>%
  dplyr::select(phylo_entropy,unrooted_pd,bwpd,isolation_source,nurse_sow,maternal_sow,
                cross_breed,BIRTH_DAY,LINE)


# plots

cw_summary <- startDF1 %>% 
  group_by(nurse_sow) %>% 
  tally()

a <- ggplot(startDF1, aes(x=nurse_sow, y=bwpd, group=nurse_sow)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by nurse sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

b <- ggplot(startDF1, aes(x=nurse_sow, y=unrooted_pd, group=nurse_sow)) + 
  labs(title = "Piglets alpha diversity (unrooted)",
       subtitle = "Grouped by nurse sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(nurse_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(200,400))

###

cw_summary <- startDF1 %>% 
  group_by(maternal_sow) %>% 
  tally()

c <- ggplot(startDF1, aes(x=maternal_sow, y=bwpd, group=maternal_sow)) + 
  labs(title = "Piglets alpha diversity (BWPD)",
       subtitle = "Grouped by maternal sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(maternal_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11))

d <- ggplot(startDF1, aes(x=maternal_sow, y=unrooted_pd, group=maternal_sow)) + 
  labs(title = "Piglets alpha diversity (unrooted PD)",
       subtitle = "Grouped by maternal sow") +
  geom_boxplot() +
  geom_text(data = cw_summary,
            aes(maternal_sow, Inf, label = n), vjust="inward") +
  theme(axis.text.x=element_text(hjust=0, angle=90),
        plot.title = element_text(lineheight = 0.9, size=12),
        plot.subtitle = element_text(lineheight = 0.9, size=11)) +
  scale_y_continuous(limits=c(200,400))

###

# nurse_Sows and maternal Sows on the same plot, dividing BWPD from unrooted
pdf(paste0(out_dir,"nursesow_maternalsow_BWPD.pdf"))
ggarrange(c, a, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

pdf(paste0(out_dir,"nursesow_maternalsow_unrooted.pdf"))
ggarrange(d, b, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()


###

# look just at maternal effect

# keep only rows where mom and nurse cols are the same (pigs that stayed with their mother, no different nurse)
mteq <- startDF1[startDF1$maternal_sow==startDF1$nurse_sow, ]

# here you can see just the maternal effect, as boxplots are colored by breed and divided by birth day. 
pdf(paste0(out_dir,"maternal_effect.pdf"))
ggplot(mteq, aes(x=maternal_sow,y=unrooted_pd, color=cross_breed))+
  geom_boxplot()+
  geom_point(size=0.5)+
  theme(axis.text.x=element_text(angle=90))+
  stat_n_text(size = 3, angle=90)+
  facet_wrap( ~ as.Date(BIRTH_DAY), scales = 'free')
ggplot(mteq, aes(x=maternal_sow,y=bwpd, color=cross_breed))+
  geom_boxplot()+
  geom_point(size=0.5)+
  theme(axis.text.x=element_text(angle=90))+
  stat_n_text(size = 3, angle=90)+
  facet_wrap( ~ as.Date(BIRTH_DAY), scales = 'free')
dev.off()

########breed_bday_line_alpha######################################


# do piglets at arrival cluster by unrooted and bwpd
# based on cross_breed/line/birth day? 

startDF1$BIRTH_DAY <- as.character(startDF1$BIRTH_DAY)
startDF1$BIRTH_DAY <- factor(startDF1$BIRTH_DAY, 
                             levels=c("2017-01-06", 
                                      "2017-01-07", 
                                      "2017-01-08",
                                      "2017-01-09",
                                      "2017-01-10",
                                      "2017-01-11"))
startDF1$LINE <- as.character(startDF1$LINE)

# plots

# by cross_breed

####### get sample size within each cross_breed group:

cw_summary <- startDF1 %>% 
  group_by(cross_breed) %>% 
  tally()

# cross_breed - unrooted 
cross_breed_unrooted_plot <- ggboxplot(startDF1, x = "cross_breed", y = "unrooted_pd",
                                       color = "cross_breed", palette = "jco",
                                       add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(200,400)+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# cross_breed - bwpd

cross_breed_bwpd_plot <- ggboxplot(startDF1, x = "cross_breed", y = "bwpd",
                                   color = "cross_breed", palette = "jco",
                                   add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(),legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward")+
  stat_compare_means(method = "kruskal.test", label.y=1.5) 

tosave <- ggarrange(cross_breed_unrooted_plot, cross_breed_bwpd_plot, nrow = 2, labels=c("A","B"))
ggsave(file = paste0(out_dir,"breed_alpha.pdf"), tosave)


# by birth day

####### get sample size within each bday group:

cw_summary <- startDF1 %>% 
  group_by(BIRTH_DAY) %>% 
  tally()

# bday - unrooted 
bday_unrooted_plot <- ggboxplot(startDF1, x = "BIRTH_DAY", y = "unrooted_pd",
                                color = "BIRTH_DAY", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(200,400)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# bday - bwpd 
bday_bwpd_plot <- ggboxplot(startDF1, x = "BIRTH_DAY", y = "bwpd",
                            color = "BIRTH_DAY", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.7)+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value


tosave <- ggarrange(bday_unrooted_plot, bday_bwpd_plot, nrow = 2, labels=c("A","B"))

pdf(paste0(out_dir,"bday_alpha.pdf"))
tosave
dev.off()


# by line

####### get sample size within each bday group:

cw_summary <- startDF1 %>% 
  group_by(LINE) %>% 
  tally()

# line - unrooted 
LINE_unrooted_plot <- ggboxplot(startDF1, x = "LINE", y = "unrooted_pd",
                                color = "LINE", palette = "jco",
                                add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  ylim(200,400)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=50)  # Add pairwise comparisons p-value

# line - bwpd 
LINE_bwpd_plot <- ggboxplot(startDF1, x = "LINE", y = "bwpd",
                            color = "LINE", palette = "jco",
                            add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="none")+
  ylim(1.5,2.75)+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.y=1.5)  # Add pairwise comparisons p-value


tosave <- ggarrange(LINE_unrooted_plot, LINE_bwpd_plot, nrow = 2, labels=c("A","B"))

pdf(paste0(out_dir,"line_alpha.pdf"))
tosave
dev.off()


###
# bday AND cross_breed:
# with age unrooted pd decreases and bwpd increases 
# conclusion comes from two cross_breeds as each bday is not represented 
# equally by the 4 cross_breeds 
my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

startDF1_sub <- startDF1 %>%
  dplyr::filter(!cross_breed=="Landrace x Cross bred (LW x D)")
startDF1_sub <- startDF1_sub %>%
  dplyr::filter(!cross_breed=="Large white x Duroc")


p1 <- ggboxplot(startDF1_sub, x = "BIRTH_DAY", y = "unrooted_pd",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(200,500)+
  stat_compare_means(comparisons = my_comparisons)

p2 <- ggboxplot(startDF1_sub, x = "BIRTH_DAY", y = "bwpd",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  theme_bw()+
  theme(axis.text.x=element_blank())+
  ylim(1.6,3.4)+
  stat_compare_means(comparisons = my_comparisons)


tosave <- ggarrange(p1,p2, nrow = 2, labels=c("A","B"))

pdf(paste0(out_dir,"bday_bybreed_alpha.pdf"))
tosave
dev.off()




########nursesow_beta#########################################################################

# 6   # plot BETA diversity (at pig trial start) 

# BETA diversity in piglets at the start of the trial 

startDF <- coggo %>% filter(
  collection_date == "t0" )  %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source)

startDF1 <- merge(startDF,details, by="isolation_source")

startDF1 <- startDF1 %>%
  dplyr::select(pc1,pc2,pc3,pc4,pc5,isolation_source,nurse_sow,maternal_sow)

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="right",
             legend.text=element_text(size=8),
             legend.title=element_text(size=8),
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

# nurse_Sows

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse_sow is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(nurse_sow) | duplicated(nurse_sow, fromLast=TRUE))

# order alphabetically by nurse_sow
startDF1_unique <- startDF1_unique[order(startDF1_unique$nurse_sow),]

rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$nurse_sow))

# subsetting to have 5 nurse_Sows in each 
startDF11 <- startDF1_unique[1:27,]
startDF12 <- startDF1_unique[28:47,]
startDF13 <- startDF1_unique[48:63,]
startDF14 <- startDF1_unique[64:80,]
startDF15 <- startDF1_unique[81:98,]
startDF16 <- startDF1_unique[99:122,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme #+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme #+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=nurse_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)


pdf(paste0(out_dir,"nursesow_PC1PC2.pdf"))
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()

####

# maternal Sows

startDF1_unique <- startDF1 %>% group_by(isolation_source) %>% slice(1)

# remove rows were nurse_sow is unique (can't be plotted a pca)
startDF1_unique <- subset(startDF1_unique,duplicated(maternal_sow) | duplicated(maternal_sow, fromLast=TRUE))

# order alphabetically by nurse_sow
startDF1_unique <- startDF1_unique[order(startDF1_unique$maternal_sow),]
rownames(startDF1_unique) <- 1:nrow(startDF1_unique)
length(unique(startDF1_unique$maternal_sow))

# subsetting to have 5 maternal Sows in each 
startDF11 <- startDF1_unique[1:30,]
startDF12 <- startDF1_unique[31:58,]
startDF13 <- startDF1_unique[59:72,]
startDF14 <- startDF1_unique[73:85,]
startDF15 <- startDF1_unique[86:108,]
startDF16 <- startDF1_unique[109:123,]

p11<-ggplot(startDF11,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p12<-ggplot(startDF12,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p13<-ggplot(startDF13,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p14<-ggplot(startDF14,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p15<-ggplot(startDF15,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)
p16<-ggplot(startDF16,aes(x=pc1,y=pc2,color=maternal_sow ))+
  geom_point()+
  theme#+
  #stat_ellipse(inherit.aes = TRUE, level = 0.80)


pdf(paste0(out_dir,"maternalsow_PC1PC2.pdf"))
grid.arrange(p11,p12,p13,p14,p15,p16, nrow=3,ncol=2)
dev.off()


#######breed_bday_line_BETA################################

# do piglets at arrival cluster by PCA
# based on cross_breed/line/birth day? 

startDF2 <- merge(startDF,details, by="isolation_source")

startDF2$BIRTH_DAY <- as.character(startDF2$BIRTH_DAY)
startDF2$BIRTH_DAY <- factor(startDF2$BIRTH_DAY, 
                             levels=c("2017-01-06", 
                                      "2017-01-07", 
                                      "2017-01-08",
                                      "2017-01-09",
                                      "2017-01-10",
                                      "2017-01-11"))

startDF2$LINE <- as.character(startDF2$LINE)


# by cross_breed

####### get sample size within each cross_breed group:

cw_summary <- startDF2 %>% 
  group_by(cross_breed) %>% 
  tally()

# cross_breed - PCA
cross_breed_PC1_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc1",
                                  color = "cross_breed", palette = "jco",
                                  add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.5)  # Add pairwise comparisons p-value

cross_breed_PC2_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc2",
                                  color = "cross_breed", palette = "jco",
                                  add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value

cross_breed_PC3_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc3",
                                  color = "cross_breed", palette = "jco",
                                  add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value

cross_breed_PC4_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc4",
                                  color = "cross_breed", palette = "jco",
                                  add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.2)  # Add pairwise comparisons p-value

cross_breed_PC5_plot <- ggboxplot(startDF2, x = "cross_breed", y = "pc5",
                                  color = "cross_breed", palette = "jco",
                                  add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(cross_breed, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value



#tosave <- ggarrange(cross_breed_PC1_plot, cross_breed_PC2_plot, cross_breed_PC3_plot, cross_breed_PC4_plot, 
#          cross_breed_PC5_plot, nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
#          common.legend = TRUE)
#ggsave(file = "out/cross_breed_beta.pdf", tosave)

# by birth day

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(BIRTH_DAY) %>% 
  tally()

# bday - PCA
bday_PC1_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc1",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value

bday_PC2_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc2",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=1)  # Add pairwise comparisons p-value

bday_PC3_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc3",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.5)  # Add pairwise comparisons p-value

bday_PC4_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc4",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value

bday_PC5_plot <- ggboxplot(startDF2, x = "BIRTH_DAY", y = "pc5",
                           color = "BIRTH_DAY", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(BIRTH_DAY, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value



tosave <- ggarrange(bday_PC1_plot, bday_PC2_plot, bday_PC3_plot, bday_PC4_plot, 
                    bday_PC5_plot, nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
                    common.legend = TRUE)

pdf(paste0(out_dir,"bday_beta.pdf"))
tosave
dev.off()



# by line

####### get sample size within each bday group:

cw_summary <- startDF2 %>% 
  group_by(LINE) %>% 
  tally()

# line - PCA
line_PC1_plot <- ggboxplot(startDF2, x = "LINE", y = "pc1",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=0.2)  # Add pairwise comparisons p-value

line_PC2_plot <- ggboxplot(startDF2, x = "LINE", y = "pc2",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=.2)  # Add pairwise comparisons p-value

line_PC3_plot <- ggboxplot(startDF2, x = "LINE", y = "pc3",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1.4)  # Add pairwise comparisons p-value

line_PC4_plot <- ggboxplot(startDF2, x = "LINE", y = "pc4",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=-1)  # Add pairwise comparisons p-value

line_PC5_plot <- ggboxplot(startDF2, x = "LINE", y = "pc5",
                           color = "LINE", palette = "jco",
                           add = "jitter")+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position="top")+
  geom_text(data = cw_summary,
            aes(LINE, Inf, label = n), vjust="inward") +
  ylim(4.5,6.3)+
  stat_compare_means(method = "kruskal.test", label.x=1, label.y=4.5)  # Add pairwise comparisons p-value



# tosave <- ggarrange(line_PC1_plot, line_PC2_plot, line_PC3_plot, line_PC4_plot, 
#                     line_PC5_plot,nrow = 3, ncol=2, labels = c("A","B","C","D","E"),
#                     common.legend = TRUE)
# ggsave(file = "out/line_beta.pdf", tosave)


###

# bday AND cross_breed:

my_comparisons <- list(  c("2017-01-07", "2017-01-09"),  
                         c("2017-01-08", "2017-01-10"), 
                         c("2017-01-09", "2017-01-11"), 
                         c("2017-01-10", "2017-01-11") )

# need to exclude two cross_breeds as these two cross_breeds don't have enough
# age groups to be plotted and compared 
startDF2_sub <- startDF2 %>%
  filter(!cross_breed=="Landrace x Cross bred (LW x D)")
startDF2_sub <- startDF2_sub %>%
  filter(!cross_breed=="Large white x Duroc")


p1 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc1",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(0,4)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p2 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc2",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(0,7)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p3 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc3",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(-2,3)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p4 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc4",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(-1,2.2)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)
p5 <- ggboxplot(startDF2_sub, x = "BIRTH_DAY", y = "pc5",
                color = "BIRTH_DAY", palette = "jco",
                add = "jitter",
                facet.by = "cross_breed", short.panel.labs = FALSE) +
  ylim(4.7,6.6)+
  theme_bw()+
  theme(axis.text.x=element_blank(), legend.position = 'right')+
  stat_compare_means(comparisons = my_comparisons)

pdf(paste0(out_dir,"bday_bybreed_beta.pdf"))
grid.arrange(
  p1,p2, nrow = 2, ncol = 1
)
grid.arrange(
  p3,p4, nrow = 2, ncol = 1
)
grid.arrange(
  p5, nrow = 2, ncol = 1
)
dev.off()


#######distribution_breeds_bdays###########################################################################



# distribution of cross_breeds, birth days across cohorts

finalDF <- inner_join(boggo,coggo)
NROW(finalDF)

df <- merge(finalDF,details, by="isolation_source")
head(df)

# distribution of cross_breeds and bdays across cohorts
df1 <- setDT(df)[, .(Freq = .N), by = .(BIRTH_DAY,cross_breed,Cohort)]
df1[order(df1$cross_breed)]

p1 <- ggplot(df1, aes(fill=BIRTH_DAY, y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  theme(legend.position="none")+
  facet_wrap(.~cross_breed,scales="free")


# distribution of Cohorts across bdays
df1 <- setDT(df)[, .(Freq = .N), by = .(Cohort,BIRTH_DAY)]
df1[order(df1$Cohort)]

p2 <- ggplot(df1, aes(fill=sort(BIRTH_DAY), y=Freq, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Cohort",
       y = "number of samples",
       fill = "birth day") +
  theme_bw()+
  theme(legend.position="right",
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(),
        axis.title.y=element_text())


pdf(paste0(out_dir,"distribution_breeds&bday_cohorts.pdf"))
ggarrange(p1,p2,nrow=2,labels=c("A","B"))
dev.off()



#######pvalues######################################################################



# 8   # p-values

# PVALUES - 


# merge alpha&beta div (finalDF) to details and details metadata

df <- merge(finalDF,details, by="isolation_source")
head(df)

unique(df$collection_date)


NROW(df)
df <- na.omit(df, cols = c("Cohort","collection_date"))
NROW(df)

head(df)

df1 <- df %>%
  dplyr::select(unrooted_pd,bwpd,pc1,pc2,pc3,pc4,pc5,Cohort,collection_date,isolation_source,
                BIRTH_DAY,cross_breed,LINE,maternal_sow,nurse_sow)
NROW(df1)

# for some reasons df1$pc2 is character and not numeric. convert: 
df1$pc2 <- as.numeric(df1$pc2)

# aggregating by avg (unique samples kept)
cols <- 1:7
df1 <- setDT(df1)[, lapply(.SD, mean), by=c(names(df1)[8:15]), .SDcols=cols]
NROW(df1)


df_cross_breed <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$cross_breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$cross_breed)$p.value,
      pc1=kruskal.test(.$pc1, .$cross_breed)$p.value,
      pc2=kruskal.test(.$pc2, .$cross_breed)$p.value,
      pc3=kruskal.test(.$pc3, .$cross_breed)$p.value,
      pc4=kruskal.test(.$pc4, .$cross_breed)$p.value,
      pc5=kruskal.test(.$pc5, .$cross_breed)$p.value,
      grouping=paste0("cross_breed"),
      stringsAsFactors=FALSE)
  })

df_cross_breed_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$cross_breed)$p.value,
      bwpd=kruskal.test(.$bwpd, .$cross_breed)$p.value,
      pc1=kruskal.test(.$pc1, .$cross_breed)$p.value,
      pc2=kruskal.test(.$pc2, .$cross_breed)$p.value,
      pc3=kruskal.test(.$pc3, .$cross_breed)$p.value,
      pc4=kruskal.test(.$pc4, .$cross_breed)$p.value,
      pc5=kruskal.test(.$pc5, .$cross_breed)$p.value,
      grouping=paste0("cross_breed"),
      stringsAsFactors=FALSE)
  }) 

df_line <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$LINE)$p.value,
      bwpd=kruskal.test(.$bwpd, .$LINE)$p.value,
      pc1=kruskal.test(.$pc1, .$LINE)$p.value,
      pc2=kruskal.test(.$pc2, .$LINE)$p.value,
      pc3=kruskal.test(.$pc3, .$LINE)$p.value,
      pc4=kruskal.test(.$pc4, .$LINE)$p.value,
      pc5=kruskal.test(.$pc5, .$LINE)$p.value,
      grouping=paste0("line"),
      stringsAsFactors=FALSE)
  })

df_line_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$LINE)$p.value,
      bwpd=kruskal.test(.$bwpd, .$LINE)$p.value,
      pc1=kruskal.test(.$pc1, .$LINE)$p.value,
      pc2=kruskal.test(.$pc2, .$LINE)$p.value,
      pc3=kruskal.test(.$pc3, .$LINE)$p.value,
      pc4=kruskal.test(.$pc4, .$LINE)$p.value,
      pc5=kruskal.test(.$pc5, .$LINE)$p.value,
      grouping=paste0("line"),
      stringsAsFactors=FALSE)
  }) 

df_bday <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day"),
      stringsAsFactors=FALSE)
  })


# cross_breeds
"Landrace x Cross bred (LW x D)"
"Duroc x Landrace"
"Duroc x Large white"
"Large white x Duroc"

df_bday_DurocxLandrace <- df1[df1$cross_breed=="Duroc x Landrace",] %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day - Duroc x Landrace"),
      stringsAsFactors=FALSE)
  })

df_bday_DurocxLw <- df1[df1$cross_breed=="Duroc x Large white",] %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day - Duroc x Large white"),
      stringsAsFactors=FALSE)
  })

df_bday_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$BIRTH_DAY)$p.value,
      bwpd=kruskal.test(.$bwpd, .$BIRTH_DAY)$p.value,
      pc1=kruskal.test(.$pc1, .$BIRTH_DAY)$p.value,
      pc2=kruskal.test(.$pc2, .$BIRTH_DAY)$p.value,
      pc3=kruskal.test(.$pc3, .$BIRTH_DAY)$p.value,
      pc4=kruskal.test(.$pc4, .$BIRTH_DAY)$p.value,
      pc5=kruskal.test(.$pc5, .$BIRTH_DAY)$p.value,
      grouping=paste0("birth day"),
      stringsAsFactors=FALSE)
  }) 

df_maternal_sow <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$maternal_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$maternal_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$maternal_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$maternal_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$maternal_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$maternal_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$maternal_sow)$p.value,
      grouping=paste0("maternal_sow"),
      stringsAsFactors=FALSE)
  }) 

df_maternal_sow_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$maternal_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$maternal_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$maternal_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$maternal_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$maternal_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$maternal_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$maternal_sow)$p.value,
      grouping=paste0("maternal_sow"),
      stringsAsFactors=FALSE)
  }) 

df_nurse_sow <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse_sow)$p.value,
      grouping=paste0("nurse_sow"),
      stringsAsFactors=FALSE)
  })

df_nurse_sow_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$nurse_sow)$p.value,
      bwpd=kruskal.test(.$bwpd, .$nurse_sow)$p.value,
      pc1=kruskal.test(.$pc1, .$nurse_sow)$p.value,
      pc2=kruskal.test(.$pc2, .$nurse_sow)$p.value,
      pc3=kruskal.test(.$pc3, .$nurse_sow)$p.value,
      pc4=kruskal.test(.$pc4, .$nurse_sow)$p.value,
      pc5=kruskal.test(.$pc5, .$nurse_sow)$p.value,
      grouping=paste0("nurse_sow"),
      stringsAsFactors=FALSE)
  })

df_Cohort <- df1 %>%
  group_by(collection_date) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$Cohort)$p.value,
      bwpd=kruskal.test(.$bwpd, .$Cohort)$p.value,
      pc1=kruskal.test(.$pc1, .$Cohort)$p.value,
      pc2=kruskal.test(.$pc2, .$Cohort)$p.value,
      pc3=kruskal.test(.$pc3, .$Cohort)$p.value,
      pc4=kruskal.test(.$pc4, .$Cohort)$p.value,
      pc5=kruskal.test(.$pc5, .$Cohort)$p.value,
      grouping=paste0("cohorts"),
      stringsAsFactors=FALSE)
  })

df_Cohort_all <- df1 %>%
  do({
    data.frame(
      collection_date=paste0("all"),
      sample_size=NROW(.),
      unrooted_pd=kruskal.test(.$unrooted_pd, .$Cohort)$p.value,
      bwpd=kruskal.test(.$bwpd, .$Cohort)$p.value,
      pc1=kruskal.test(.$pc1, .$Cohort)$p.value,
      pc2=kruskal.test(.$pc2, .$Cohort)$p.value,
      pc3=kruskal.test(.$pc3, .$Cohort)$p.value,
      pc4=kruskal.test(.$pc4, .$Cohort)$p.value,
      pc5=kruskal.test(.$pc5, .$Cohort)$p.value,
      grouping=paste0("cohorts"),
      stringsAsFactors=FALSE)
  })

df_cross_breed_all <- as.data.frame(df_cross_breed_all)
df_cross_breed <- as.data.frame(df_cross_breed)
df_line_all <- as.data.frame(df_line_all)
df_line <- as.data.frame(df_line)
df_bday_all <- as.data.frame(df_bday_all)
df_bday <- as.data.frame(df_bday)
df_bday_DurocxLandrace <- as.data.frame(df_bday_DurocxLandrace)
df_bday_DurocxLw <- as.data.frame(df_bday_DurocxLw)
df_maternal_sow_all <- as.data.frame(df_maternal_sow_all)
df_maternal_sow <- as.data.frame(df_maternal_sow)
df_nurse_sow_all <- as.data.frame(df_nurse_sow_all)
df_nurse_sow <- as.data.frame(df_nurse_sow)
df_Cohort_all <- as.data.frame(df_Cohort_all)

all_pvalues <- rbind(df_cross_breed_all, df_cross_breed,
                     df_line_all, df_line, 
                     df_bday_all, df_bday, df_bday_DurocxLandrace, df_bday_DurocxLw, 
                     df_maternal_sow_all, df_maternal_sow, 
                     df_nurse_sow_all, df_nurse_sow, 
                     df_Cohort_all)

all_pvalues$test <- "Kruskal-Wallis"


# write out in workbook
# addWorksheet(wb, "all_pvalues")
# writeData(wb, sheet = "all_pvalues", all_pvalues, rowNames = FALSE)
fwrite(x=all_pvalues, file=paste0(stats_dir,"all_pvalues.csv"))

padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

df_cross_breed <- df_cross_breed %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_line <- df_line %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday <- df_bday %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLandrace <- df_bday_DurocxLandrace %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_bday_DurocxLw <- df_bday_DurocxLw %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_nurse_sow <- df_nurse_sow %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 

df_maternal_sow <- df_maternal_sow %>%
  dplyr::mutate_at(c("unrooted_pd","bwpd","pc1","pc2","pc3","pc4","pc5"),padj_function) 


all_padj_Hommel <- rbind(df_cross_breed,
                         df_line, 
                         df_bday, 
                         df_bday_DurocxLandrace, 
                         df_bday_DurocxLw, 
                         df_maternal_sow, 
                         df_nurse_sow)


all_padj_Hommel$test <- "Kruskal-Wallis"
all_padj_Hommel$padj_method <- "Hommel"

# write out in workbook
# addWorksheet(wb, "all_padj_Hommel")
# writeData(wb, sheet = "all_padj_Hommel", all_padj_Hommel, rowNames = FALSE)
fwrite(x=all_padj_Hommel, file=paste0(stats_dir,"all_padj_Hommel.csv"))


# adjusted pvalues

# chosen method is Tukey: 

# When you do Tukeys test, the variance is estimated from the whole set of data 
# as a pooled estimate. If the population variances are the 
# same in all groups, such a pooled estimate is much more robust and precise 
# than the individual estimated from just a part of the whole set of data. 
# Further, Tukeys procedure adjusts the p-values for multiple testing, so that 
# the family-wise error rate is controlled (probability to get at least one false 
# positive among the family of tests performed).


# by cross_breed

aov.out = aov(unrooted_pd ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ cross_breed, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$cross_breed)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_cross_breed <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_cross_breed$group = "cross_breed"


# by LINE

# to character otherwise considered numeric
df1$LINE <- as.character(df1$LINE)

aov.out = aov(unrooted_pd ~ LINE, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ LINE, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$LINE)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_LINE <- rbind(aov.out1,
                 aov.out2,
                 aov.out3,
                 aov.out4,
                 aov.out5,
                 aov.out6,
                 aov.out7)
by_LINE$group = "LINE"



# by BIRTH_DAY

# to character otherwise considered numeric
df1$BIRTH_DAY <- as.character(df1$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_BIRTH_DAY$group = "BIRTH_DAY"


# by BIRTH_DAY for cross_breed "Duroc x Landrace"

df1_sub <- df1[df1$cross_breed=="Duroc x Landrace",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Duroc_x_Landrace <- rbind(aov.out1,
                                       aov.out2,
                                       aov.out3,
                                       aov.out4,
                                       aov.out5,
                                       aov.out6,
                                       aov.out7)
by_BIRTH_DAY_Duroc_x_Landrace$group = "BIRTH_DAY_Duroc_x_Landrace"


# by BIRTH_DAY for cross_breed "Duroc x Large white"

df1_sub <- df1[df1$cross_breed=="Duroc x Large white",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Duroc_x_Large_white <- rbind(aov.out1,
                                          aov.out2,
                                          aov.out3,
                                          aov.out4,
                                          aov.out5,
                                          aov.out6,
                                          aov.out7)
by_BIRTH_DAY_Duroc_x_Large_white$group = "BIRTH_DAY_Duroc_x_Large_white"


# by BIRTH_DAY for cross_breed "Large white x Duroc"

df1_sub <- df1[df1$cross_breed=="Large white x Duroc",]

# to character otherwise considered numeric
df1_sub$BIRTH_DAY <- as.character(df1_sub$BIRTH_DAY)

aov.out = aov(unrooted_pd ~ BIRTH_DAY, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ BIRTH_DAY, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$BIRTH_DAY)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_BIRTH_DAY_Large_white_x_Duroc <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)
by_BIRTH_DAY_Large_white_x_Duroc$group = "BIRTH_DAY_Large_white_x_Duroc"


# not enough timepoint for "Landrace x Cross bred (LW x D)"


# by nurse_sow

# to character otherwise considered numeric
df1$nurse_sow <- as.character(df1$nurse_sow)


aov.out = aov(unrooted_pd ~ nurse_sow, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ nurse_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$nurse_sow)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_nurse_sow <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_nurse_sow$group = "nurse_sow"

# by maternal_sow

# to character otherwise considered numeric
df1$maternal_sow <- as.character(df1$maternal_sow)


aov.out = aov(unrooted_pd ~ maternal_sow, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ maternal_sow, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$maternal_sow)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_maternal_sow <- rbind( 
  aov.out1,
  aov.out2,
  aov.out3,
  aov.out4,
  aov.out5,
  aov.out6,
  aov.out7)

by_maternal_sow$group = "maternal_sow"

#######AOV_by_cohort####

# by Cohort

aov.out = aov(unrooted_pd ~ Cohort, data=df1)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort <- rbind(aov.out1,
                   aov.out2,
                   aov.out3,
                   aov.out4,
                   aov.out5,
                   aov.out6,
                   aov.out7)
by_Cohort$group = "Cohort"


# by Cohort t0

df1_sub <- df1[df1$collection_date=="t0",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t0 <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_Cohort_t0$group = "Cohort_t0"


# by Cohort t1

df1_sub <- df1[df1$collection_date=="t1",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t1 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t1$group = "Cohort_t1"



# by Cohort t2

df1_sub <- df1[df1$collection_date=="t2",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t2 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t2$group = "Cohort_t2"



# by Cohort t3

df1_sub <- df1[df1$collection_date=="t3",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t3 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t3$group = "Cohort_t3"



# by Cohort t4

df1_sub <- df1[df1$collection_date=="t4",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t4 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t4$group = "Cohort_t4"



# by Cohort t5

df1_sub <- df1[df1$collection_date=="t5",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t5 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t5$group = "Cohort_t5"



# by Cohort t6

df1_sub <- df1[df1$collection_date=="t6",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t6 <- rbind(aov.out1,
                        aov.out2,
                        aov.out3,
                        aov.out4,
                        aov.out5,
                        aov.out6,
                        aov.out7)
by_Cohort_t6$group = "Cohort_t6"



# by Cohort t7

df1_sub <- df1[df1$collection_date=="t7",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t7 <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_Cohort_t7$group = "Cohort_t7"



# by Cohort i6

df1_sub <- df1[df1$collection_date=="t8",]

aov.out = aov(unrooted_pd ~ Cohort, data=df1_sub)   
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out1 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out1$type="unrooted_pd"
#
aov.out = aov(bwpd ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out2 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out2$type="bwpd"
#
aov.out = aov(pc1 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out3 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out3$type="pc1"
#
aov.out = aov(pc2 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out4 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out4$type="pc2"
# 
aov.out = aov(pc3 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out5 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out5$type="pc3"
#
aov.out = aov(pc4 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out6 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out6$type="pc4"
#
aov.out = aov(pc5 ~ Cohort, data=df1_sub)
res <- TukeyHSD(aov.out)
aov.out <- as.data.frame(res$Cohort)
aov.out7 <- tibble::rownames_to_column(aov.out, "comparison")
aov.out7$type="pc5"

by_Cohort_t8 <- rbind(aov.out1,
                      aov.out2,
                      aov.out3,
                      aov.out4,
                      aov.out5,
                      aov.out6,
                      aov.out7)
by_Cohort_t8$group = "Cohort_t8"

all_padj_Tukey <- rbind(by_cross_breed,
                        by_LINE, 
                        by_BIRTH_DAY, 
                        by_BIRTH_DAY_Duroc_x_Landrace,
                        by_BIRTH_DAY_Duroc_x_Large_white,
                        by_BIRTH_DAY_Large_white_x_Duroc, 
                        by_nurse_sow,
                        by_maternal_sow,
                        by_Cohort, 
                        by_Cohort_t0,
                        by_Cohort_t1,
                        by_Cohort_t2,
                        by_Cohort_t3,
                        by_Cohort_t4,
                        by_Cohort_t5,
                        by_Cohort_t6,
                        by_Cohort_t7,
                        by_Cohort_t8
)

all_padj_Tukey$test <- "anova"
all_padj_Tukey$padj_method <- "TukeyHSD"

# write out in workbook
# addWorksheet(wb, "all_padj_Tukey")
# writeData(wb, sheet = "all_padj_Tukey", all_padj_Tukey, rowNames = FALSE)



# plot p-values for start factors

piglets_factors <- all_pvalues %>%
  dplyr::filter(grouping != "cohorts" &
                  grouping != "ctrl_neo" &
                  grouping != "DScour_ColiGuard" &
                  grouping != "NeoD_NeoC" &
                  collection_date != "all") 

piglets_factors$grouping <- gsub("birth day - Duroc x Large white",
                                 "birth day - DxLW",piglets_factors$grouping)
piglets_factors$grouping <- gsub("birth day - Duroc x Landrace",
                                 "birth day - DxL",piglets_factors$grouping)

piglets_factors2 <- all_padj_Hommel %>%
  filter(grouping != "cohorts" &
           grouping != "ctrl_neo" &
           grouping != "Dscour_ColiGuard" &
           grouping != "NeoD_NeoC" &
           collection_date != "all") 

piglets_factors2$grouping <- gsub("birth day - Duroc x Large white",
                                  "birth day - DxLW",piglets_factors2$grouping)
piglets_factors2$grouping <- gsub("birth day - Duroc x Landrace",
                                  "birth day - DxL",piglets_factors2$grouping)


unique(piglets_factors$grouping)
unique(piglets_factors2$grouping)

# order to show facets:
piglets_factors$grouping <- factor(piglets_factors$grouping,
                                   levels=c("birth day",
                                            "birth day - DxL",
                                            "birth day - DxLW",
                                            "cross_breed",
                                            "maternal_sow",
                                            "nurse_sow",
                                            "line"))
piglets_factors2$grouping <- factor(piglets_factors2$grouping,
                                    levels=c("birth day",
                                             "birth day - DxL",
                                             "birth day - DxLW",
                                             "cross_breed",
                                             "maternal_sow",
                                             "nurse_sow",
                                             "line"))

# alpha
df <- piglets_factors %>%
  pivot_longer(cols=c('unrooted_pd','bwpd'))
colnames(df)[colnames(df)=="name"] <- "parameter"
df2 <- piglets_factors2 %>%
  pivot_longer(cols=c('unrooted_pd','bwpd'))
colnames(df2)[colnames(df2)=="name"] <- "parameter"
df$grouping <- gsub(" - ","\n",df$grouping)
df2$grouping <- gsub(" - ","\n",df2$grouping)
df$grouping <- gsub("_","\n",df$grouping)
df2$grouping <- gsub("_","\n",df2$grouping)
# re-order dates
df2$collection_date <- factor(df2$collection_date, 
                           levels=c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))
alpha_plot <- ggplot(df, aes(x=collection_date,y=value)) + 
  ylim(0,0.06)+
  labs(y="alpha diversity - p-value",
       x="")+
  theme_bw()+
  geom_point(data=df2,aes(shape=parameter), color="red", size=2)+
  geom_point(aes(shape=parameter), color="black", size=2)+
  theme(axis.text.x=element_text(size=3),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        legend.position="right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))+
  facet_grid(~grouping)+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)
# beta 
df <- piglets_factors %>%
  pivot_longer(cols=c(contains('pc')))
colnames(df)[colnames(df)=="name"] <- "parameter"
df2 <- piglets_factors2 %>%
  pivot_longer(cols=c(contains('pc')))
colnames(df2)[colnames(df2)=="name"] <- "parameter"
df$grouping <- gsub(" - ","\n",df$grouping)
df2$grouping <- gsub(" - ","\n",df2$grouping)
df$grouping <- gsub("_","\n",df$grouping)
df2$grouping <- gsub("_","\n",df2$grouping)
# re-order dates
df2$collection_date <- factor(df2$collection_date, 
                              levels=c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))
beta_plot <- ggplot(df, aes(x=collection_date,y=value)) + 
  ylim(0,0.06)+
  labs(y="beta diversity - p-value",
       x="")+
  theme_bw()+
  geom_point(data=df2,aes(shape=parameter), color="red", size=2)+
  geom_point(aes(shape=parameter), color="black", size=2)+
  theme(axis.text.x=element_text(size=3),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),
        legend.position="right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))+
  facet_grid(~grouping)+
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "black", size=0.5)

start_factors_pvalues_plot <- ggarrange(
  alpha_plot, beta_plot,ncol=1,nrow=2, labels = c("A","B")
)

pdf(paste0(out_dir,"start_factors_pvalues.pdf"))
start_factors_pvalues_plot
dev.off()


###########################################################################################

# comment PhD Examiner 3: 
# p95 A plot where only p-values are shown. 
# Plotting estimates of the diversity indices 
# (or the difference in diversity indices) 
# along with confidence intervals is vastly preferable.

# EXTRA : as asked by Examiner 3: 

head(all_padj_Hommel)
out <- all_padj_Hommel %>%
  pivot_longer(cols=3:9) %>%
  dplyr::filter(value <= 0.05) %>%
  dplyr::filter(!grouping=="birth day - Duroc x Landrace")
head(out)

out$grouping <- gsub(pattern = "line","LINE", out$grouping)
out$grouping <- gsub(pattern = "birth day","BIRTH_DAY", out$grouping)


head(df1)
df1$BIRTH_DAY <- as.character(df1$BIRTH_DAY)
df_final <- data.frame()
for (row in 1:nrow(out)) {
  
  var_to_sel <- out[row,6]
  coll_date <- out[row,1]
  group_sel <- as.data.frame(out[row,3])
  
  test <- df1[df1$collection_date %in% coll_date,]
  
  test <- test %>%
    dplyr::select(group_sel$grouping,isolation_source,collection_date,var_to_sel$name)
  test <- as.data.frame(test)
  
  test$id <- paste0("set_",row)
  test$grouping <- paste0(as.character(group_sel$grouping))
  test$variable=paste0(var_to_sel)
  
  colnames(test) <- c("spec","isolation_source","collection_date","variable_values","id","grouping","variable")
  
  df_final <- rbind(df_final,test)
  
}


# splitting into multiple dataframes (by set)
multi_df <- split( df_final , f = df_final$id )

plot_me <- function(split_df) {
  
  # df as dataframe
  df0 <- as.data.frame(split_df)
  
  # plot
  return(print(df0 %>% 
                 ggplot(., aes(x=spec,y=variable_values))+
                 geom_boxplot(lwd=0.2, outlier.size = 0.5)+
                 stat_n_text(size = 1.5)+
                 coord_flip()+
                 ylab(paste0(df0$variable))+
                 xlab(paste0(df0$grouping))+
                 theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
                       axis.text.x=element_text(size=4),
                       axis.title.x=element_text(size=8),
                       axis.title.y=element_text(size=8))+
                 facet_grid(~collection_date, scale="free")))
}


View(df_final)


p1<-plot_me(df_final %>% dplyr::filter(id=="set_1")) # 
p2<-plot_me(df_final %>% dplyr::filter(id=="set_2")) # 
p3<-plot_me(df_final %>% dplyr::filter(id=="set_3")) # 
p4<-plot_me(df_final %>% dplyr::filter(id=="set_4")) # 
p5<-plot_me(df_final %>% dplyr::filter(id=="set_5")) # 
p6<-plot_me(df_final %>% dplyr::filter(id=="set_6")) # 
p7<-plot_me(df_final %>% dplyr::filter(id=="set_7")) # 
p8 <- plot_me(df_final %>% dplyr::filter(id=="set_8"))  # 
p9 <- plot_me(df_final %>% dplyr::filter(id=="set_9"))  # 


# extra plot: 

out2 <- all_padj_Hommel %>%
  pivot_longer(cols=3:9) %>%
  dplyr::filter(value <= 0.05) %>%
  dplyr::filter(grouping=="birth day - Duroc x Landrace")
head(out2)

p10<-df1 %>%
  dplyr::filter(cross_breed=="Duroc x Landrace") %>%
  dplyr::filter(collection_date=="t0") %>%
  dplyr::select(isolation_source,unrooted_pd,BIRTH_DAY,collection_date,cross_breed) %>%
  ggplot(., aes(x=BIRTH_DAY,y=unrooted_pd))+
  geom_boxplot(lwd=0.2, outlier.size = 0.5)+
  stat_n_text(size = 1.5)+
  coord_flip()+
  theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
        axis.text.x=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  xlab("BIRTH_DAY DxL")+
  facet_grid(~collection_date, scale="free")

p11<-df1 %>%
  dplyr::filter(cross_breed=="Duroc x Landrace") %>%
  dplyr::filter(collection_date=="t0") %>%
  dplyr::select(isolation_source,bwpd,BIRTH_DAY,collection_date,cross_breed) %>%
  ggplot(., aes(x=BIRTH_DAY,y=bwpd))+
  geom_boxplot(lwd=0.2, outlier.size = 0.5)+
  stat_n_text(size = 1.5)+
  coord_flip()+
  theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
        axis.text.x=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  xlab("BIRTH_DAY DxL")+
  facet_grid(~collection_date, scale="free")

p12<-df1 %>%
  dplyr::filter(cross_breed=="Duroc x Landrace") %>%
  dplyr::filter(collection_date=="t0") %>%
  dplyr::select(isolation_source,pc2,BIRTH_DAY,collection_date,cross_breed) %>%
  ggplot(., aes(x=BIRTH_DAY,y=pc2))+
  geom_boxplot(lwd=0.2, outlier.size = 0.5)+
  stat_n_text(size = 1.5)+
  coord_flip()+
  theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
        axis.text.x=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  xlab("BIRTH_DAY DxL")+
  facet_grid(~collection_date, scale="free")



# put them together:

sign_plots <- ggarrange(p2,
                        p3,
                        p10,
                        p11,
                        p4,
                        p5,
                        p12,
                        p1,
                        ncol = 4, nrow=2)

pdf(paste0(out_dir,"start_factors_pvalues_boxplots.pdf"))
sign_plots
dev.off()



merged_sign_plots <- ggarrange(start_factors_pvalues_plot,
                               sign_plots, labels = c("","C"),
                               ncol = 1, nrow=2)

pdf(paste0(out_dir,"start_factors_pvalues_merged.pdf"))
merged_sign_plots
dev.off()



# 
# 
# CI_me <- function(split_df) {
#   
#   # df as dataframe
#   df0 <- as.data.frame(split_df)
#   
#   # plot
#   return(df0 %>%
#            group_by(id,collection_date,grouping,variable,spec) %>%
#            dplyr::summarise(mean.variable_values = mean(variable_values, na.rm = TRUE),
#                             sd.variable_values = sd(variable_values, na.rm = TRUE),
#                             n.variable_values = n()) %>%
#            dplyr::mutate(se.variable_values = sd.variable_values / sqrt(n.variable_values),
#                          lower.ci.variable_values = mean.variable_values - qt(1 - (0.05 / 2), n.variable_values - 1) * se.variable_values,
#                          upper.ci.variable_values = mean.variable_values + qt(1 - (0.05 / 2), n.variable_values - 1) * se.variable_values)
#   )
# }
# 
# CIs <- rbind(as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_1"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_2"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_3"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_4"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_5"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_6"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_7"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_8"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_9"))),
#              as.data.frame(CI_me(df_final %>% dplyr::filter(id=="set_10"))))



#######cite packages#######################################################################

sink(paste0(out_dir,"metapigs_phylodiv_packages_citations.bib"))
out <- sapply(names(sessionInfo()$otherPkgs), 
              function(x) print(citation(x), style = "Bibtex"))
sink()


###########################################################################################