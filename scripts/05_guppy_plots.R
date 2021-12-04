######################################################################################################

# steps before this: 

##### previous steps: 

# 1 # groups for guppy are made in guppy_group.R 

# 2 # guppy is run 

# 3 # .xml conversion to .txt:        <-  MUST RUN THIS BEFORE RUNNING THIS SCRIPT 

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy file.xml file.txt
# in a loop: 
# for fpath in  /Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output/*.xml; 
# do java -cp /Users/danielagaio/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy "$fpath" 
# "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output/$(basename "$fpath").txt"; 
# done

# 4 # .xml(txt) files are read in and parsed in guppy_XML_process.R

##### HERE : 

# 1 # .jplace files are read in and parsed

# 2 # bach effect removal 

# 3 # merges metadata

# 4 # principal component are plotted, where xml data is used to completement the taxa underlying the variation


######################################################################################################

######################################################################################################
######################################################################################################

# functions to find PC of interest


find_PC1 <- function(x) {
  PC <- x[[1]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC2 <- function(x) {
  PC <- x[[2]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC3 <- function(x) {
  PC <- x[[3]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC4 <- function(x) {
  PC <- x[[4]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC5 <- function(x) {
  PC <- x[[5]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}

######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

PC_down <- function(x) {
  down <- paste(as.list((x$taxa_simple[x$PC_position=="down"]),"\n"),collapse="\n")
  return(down)
}


PC_up <- function(x) {
  up <- paste(as.list((x$taxa_simple[x$PC_position=="up"]),"\n"),collapse="\n")
  return(up)
}

# interrogate the function by typing 
# PC_up(anydataframeyouwant)


######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

get_var <- function(x) {
  var <- paste(as.list((x$var_explained)[1]))   # first item only 
  return(var)
}

# interrogate the function by typing 
# get_var(anydataframeyouwant)


######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################


library(vcd)
library(summarytools)
library(readr)
library(splitstackshape)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(sva)   # this is ComBat , careful not to install COMBAT instead which is another thing
library(openxlsx)
library(gridExtra)
#install.packages("wordspace")
library(wordspace)

library(pheatmap)

#install.packages("taxize")
library(taxize)

library(ggplot2)
library(plotrix)

source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
guppyout_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output" # local 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/" # local 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 


###########################################################################################

# manual settings 
removebatcheffect_allowed <- "yes"     # if yes, it removes the batch effect only where detected

######################################################################################################


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

mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)

####

# these NA rows are from negative controls
mdat <- mdat[!is.na(mdat$collection_date),]
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
  dplyr::select(isolation_source,BIRTH_DAY,breed,stig,nurse)

###########################################################################################

# load XML data
simplified <- read_csv(paste0(middle_dir,"guppy_xml_simplified.df"))
simplified <- as.data.frame(simplified)

unique(simplified$sample_type)
unique(simplified$sample_type)

###########################################################################################


jplace_files = list.files(guppyout_dir,pattern=".proj")

# construct an empty dataframe to build on 
jplace_df <- data.frame(
  DNA_plate = character(),
  DNA_well = character(),
  file = character(),
  PC1 = character(),
  PC2 = character(),
  PC3 = character(),
  PC4 = character(),
  PC5 = character(),
  stringsAsFactors = FALSE
)

for (jplace_file in jplace_files) {
  
  # read in file 
  pcadat <- read_csv(file.path(guppyout_dir,jplace_file), col_names = FALSE)
  
  pcadat <- cSplit(pcadat, "X1","_")
  
  pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
  pcadat$DNA_well <- pcadat$X1_3
  colnames(pcadat)[1:5] <- c("PC1","PC2","PC3","PC4","PC5")
  pcadat <- pcadat %>%
    dplyr::select(DNA_plate,DNA_well,PC1,PC2,PC3,PC4,PC5)
  
  pcadat$file <- basename(jplace_file)
  
  jplace_df <- rbind(
    jplace_df, 
    pcadat
  )
  
}

# convert PC columns to numeric class 
jplace_df <- jplace_df %>%
  mutate_at('PC1',as.numeric) %>% 
  mutate_at('PC2',as.numeric) %>% 
  mutate_at('PC3',as.numeric) %>% 
  mutate_at('PC4',as.numeric) %>% 
  mutate_at('PC5',as.numeric) 

# clean the file names
jplace_df$file <- gsub('.proj', '', jplace_df$file)
unique(jplace_df$file)
jplace_df$file <- gsub('_sel.txt', '', jplace_df$file)
unique(jplace_df$file)

##############################
##############################

# run guppy_XML_process.R to get simplified df
# (I tried to load it with read.csv, read_csv, read.csv2, would not keep the same format!!!!! grrrrrrr)

##############################
##############################


# function to adjust p-value
padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

# determine if and where there is a batch effect
checkbatch_before <- jplace_df %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("before"),
      stringsAsFactors=FALSE)
  }) %>%
  mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 


# look at which groups have a batch effect (based on adj p-value < 0.05)
#batch_affected <- checkbatch_before %>%
#  dplyr::filter(PC1<0.05|PC2<0.05|PC3<0.05|PC4<0.05|PC5<0.05) %>%
#  dplyr::select(file)

# splitting into multiple dataframes (by file name)
multiple_DFs <- split( jplace_df , f = jplace_df$file )

# construct an empty dataframe to build on 
unbatched <- data.frame(
  DNA_plate = character(),
  DNA_well = character(),
  file = character(),
  PC1 = character(),
  PC2 = character(),
  PC3 = character(),
  PC4 = character(),
  PC5 = character(),
  stringsAsFactors = FALSE
)


##############################
# this loop is entered if manually allowed (top of script)


if (removebatcheffect_allowed=="yes") {
  for (single_DF in multiple_DFs) {
    
    DNA_plate <- single_DF$DNA_plate
    DNA_well <- single_DF$DNA_well
    file <- single_DF$file
    
    single_DF<- data.matrix(single_DF[,3:7], rownames.force = NA)
    
    single_DF<-ComBat(dat=t(as.matrix(single_DF)),DNA_plate,mod=NULL)
    
    single_DF <- t(single_DF)
    single_DF <- as.data.frame(single_DF)
    #single_DF <- unfactor(single_DF[])
    single_DF <- cbind(DNA_plate,DNA_well,file,single_DF)
    
    unbatched <- rbind(
      unbatched,
      single_DF
    )
  }
} else {
  print("No batch effect removal allowed")
}

##############################

# check batch effect AFTER batch effect removal 
checkbatch_unbatched <- unbatched %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("after"),
      stringsAsFactors=FALSE)
  }) %>%
  dplyr::mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 


##############################

# take unbatched dataframe 
jplace_df_final <- unbatched

##############################


# Time to plot! 


###########################################################################################


#settings for plots
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             plot.title = element_text(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="top",
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.6, main, adj=c(0,0), cex=1.3)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.8)
}
rbow <- rainbow(40, end=0.7, alpha=0.7)

##############################


# add one level of grouping (e.g.: all group_A* files, all timepoints, belong together)

jplace_df_final$group <- jplace_df_final$file


unique(jplace_df_final$group)

jplace_df_final$group <- gsub('pca_piggies_group_A', 'groupA', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_group_B', 'groupB', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLNEO', 'groupC', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_NEONEOD', 'groupD', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_NEONEOC', 'groupE', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLDs', 'groupF', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLC', 'groupG', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t0', 'all_t0', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t2', 'all_t2', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t4', 'all_t4', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t6', 'all_t6', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t8', 'all_t8', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_t9', 'all_t9', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies', 'all_tALL', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_pos_controls', 'pos_tNONE', jplace_df_final$group)

unique(jplace_df_final$group)

jplace_df_final <- cSplit(jplace_df_final, "group","_")
unique(jplace_df_final$group_2)

jplace_df_final <- setnames(jplace_df_final, old = c('group_1','group_2'), new = c('sample_type','guppied_date'))
head(jplace_df_final)

unique(jplace_df_final$sample_type)
unique(jplace_df_final$guppied_date)

###

##############################

# merge metadata with details (breed,bday,nurse,...)
mdat_deets <- left_join(mdat,details)

unique(jplace_df_final$file)
unique(jplace_df_final$sample_type)
unique(jplace_df_final$guppied_date)



jplace_df_final %>% 
  dplyr::filter(sample_type=="all") %>%
  group_by(guppied_date) %>%
  tally()

# merge metadata with beta diversity data
multi_coggo <- inner_join(jplace_df_final,mdat_deets, by=c("DNA_plate","DNA_well"))

multi_coggo <- multi_coggo %>%
  dplyr::select(DNA_plate,DNA_well,sample_name,isolation_source,collection_date,Cohort,breed,BIRTH_DAY,nurse,stig,sample_type,guppied_date,PC1,PC2,PC3,PC4,PC5)

##############################

# give some order to the variables 
multi_coggo$sample_type <- factor(multi_coggo$sample_type, 
                                  levels=c("pos","all","groupA","groupB",
                                           "groupC","groupD","groupE","groupF","groupG"))

unique(multi_coggo$guppied_date)
multi_coggo$guppied_date <- factor(multi_coggo$guppied_date, 
                                   levels=c("tALL",
                                            "t0",
                                            "t2",
                                            "t4",
                                            "t6",
                                            "t8",
                                            "t9",
                                            "tNONE"))
unique(multi_coggo$guppied_date)

multi_coggo$BIRTH_DAY <- factor(multi_coggo$BIRTH_DAY, 
                                levels=c("2017-01-06", 
                                         "2017-01-07", 
                                         "2017-01-08",
                                         "2017-01-09",
                                         "2017-01-10",
                                         "2017-01-11"))

##############################

# splitting into multiple dataframes (by sample_type name)
unique(multi_coggo$sample_type)

multi_coggo <- split( multi_coggo , f = multi_coggo$sample_type )
NROW(multi_coggo)

##############################


# get dataframes 


##############################
##############################

DF_positive_controls <- as.data.frame(multi_coggo$pos)

##############################
##############################

# dataframe containing the data from the single guppy run where all the pig samples where run
DF_piggies_allguppied <- as.data.frame(multi_coggo$all) %>%
  filter(guppied_date=="tALL") 

unique(DF_piggies_allguppied$sample_type)
unique(DF_piggies_allguppied$guppied_date)
unique(DF_piggies_allguppied$collection_date)

##############################
##############################

# dataframe containing the data from separate (per time point) guppy runs 
DF_piggies_time <- as.data.frame(multi_coggo$all) %>%
  filter(!guppied_date=="tALL") %>%
  filter(!guppied_date=="tNONE") 

unique(DF_piggies_time$sample_type)
unique(DF_piggies_time$guppied_date)
unique(DF_piggies_time$collection_date)

###

# separation by guppy run (each run  = 1 time point)

DF_piggies_time_t0 <- DF_piggies_time %>%
  filter(guppied_date=="t0")
DF_piggies_time_t2 <- DF_piggies_time %>%
  filter(guppied_date=="t2")
DF_piggies_time_t4 <- DF_piggies_time %>%
  filter(guppied_date=="t4")
DF_piggies_time_t6 <- DF_piggies_time %>%
  filter(guppied_date=="t6")
DF_piggies_time_t8 <- DF_piggies_time %>%
  filter(guppied_date=="t8")
DF_piggies_time_t9 <- DF_piggies_time %>%
  filter(guppied_date=="t9")


##############################
##############################

groupA <- as.data.frame(multi_coggo$groupA)
groupA_t0 <- groupA %>%
  filter(guppied_date=="t0")
groupA_t2 <- groupA %>%
  filter(guppied_date=="t2")
groupA_t4 <- groupA %>%
  filter(guppied_date=="t4")
groupA_t6 <- groupA %>%
  filter(guppied_date=="t6")
groupA_t8 <- groupA %>%
  filter(guppied_date=="t8")
groupA_t9 <- groupA %>%
  filter(guppied_date=="t9")

##############################
##############################

groupB <- as.data.frame(multi_coggo$groupB)
groupB_t0 <- groupB %>%
  filter(guppied_date=="t0")
groupB_t2 <- groupB %>%
  filter(guppied_date=="t2")
groupB_t4 <- groupB %>%
  filter(guppied_date=="t4")
groupB_t6 <- groupB %>%
  filter(guppied_date=="t6")
groupB_t8 <- groupB %>%
  filter(guppied_date=="t8")
groupB_t9 <- groupB %>%
  filter(guppied_date=="t9")

##############################
##############################

groupC <- as.data.frame(multi_coggo$groupC)
groupC_t0 <- groupC %>%
  filter(guppied_date=="t0")
groupC_t2 <- groupC %>%
  filter(guppied_date=="t2")
groupC_t4 <- groupC %>%
  filter(guppied_date=="t4")
groupC_t6 <- groupC %>%
  filter(guppied_date=="t6")
groupC_t8 <- groupC %>%
  filter(guppied_date=="t8")
groupC_t9 <- groupC %>%
  filter(guppied_date=="t9")

##############################
##############################

groupD <- as.data.frame(multi_coggo$groupD)
groupD_t0 <- groupD %>%
  filter(guppied_date=="t0")
groupD_t2 <- groupD %>%
  filter(guppied_date=="t2")
groupD_t4 <- groupD %>%
  filter(guppied_date=="t4")
groupD_t6 <- groupD %>%
  filter(guppied_date=="t6")
groupD_t8 <- groupD %>%
  filter(guppied_date=="t8")
groupD_t9 <- groupD %>%
  filter(guppied_date=="t9")

##############################
##############################

groupE <- as.data.frame(multi_coggo$groupE)
groupE_t0 <- groupE %>%
  filter(guppied_date=="t0")
groupE_t2 <- groupE %>%
  filter(guppied_date=="t2")
groupE_t4 <- groupE %>%
  filter(guppied_date=="t4")
groupE_t6 <- groupE %>%
  filter(guppied_date=="t6")
groupE_t8 <- groupE %>%
  filter(guppied_date=="t8")
groupE_t9 <- groupE %>%
  filter(guppied_date=="t9")

##############################
##############################


groupF <- as.data.frame(multi_coggo$groupF)
groupF_t0 <- groupF %>%
  filter(guppied_date=="t0")
groupF_t2 <- groupF %>%
  filter(guppied_date=="t2")
groupF_t4 <- groupF %>%
  filter(guppied_date=="t4")
groupF_t6 <- groupF %>%
  filter(guppied_date=="t6")
groupF_t8 <- groupF %>%
  filter(guppied_date=="t8")
groupF_t9 <- groupF %>%
  filter(guppied_date=="t9")

##############################
##############################

groupG <- as.data.frame(multi_coggo$groupG)
groupG_t0 <- groupG %>%
  filter(guppied_date=="t0")
groupG_t2 <- groupG %>%
  filter(guppied_date=="t2")
groupG_t4 <- groupG %>%
  filter(guppied_date=="t4")
groupG_t6 <- groupG %>%
  filter(guppied_date=="t6")
groupG_t8 <- groupG %>%
  filter(guppied_date=="t8")
groupG_t9 <- groupG %>%
  filter(guppied_date=="t9")

##############################
##############################


# PLOT! 


##############################
##############################

# positive controls

DF_positive_controls$Cohort <- factor(DF_positive_controls$Cohort, 
                                      levels=c("MockCommunity",
                                               "PosControl_D-Scour",
                                               "PosControl_ColiGuard"))

unique(simplified$sample_type)
xmldata <- simplified %>%
  filter(sample_type=="pos") %>%
  group_split(component) 

PC1PC2_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC2,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC3PC4_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC3,y=PC4,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC1PC5_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC5,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))


pdf(paste0(out_dir,"guppy_pos_controls.pdf"))
### plot PC3PC4 
PC1PC2_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC3PC4 
PC3PC4_pos_controls
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC4(xmldata)),"%"), x = unit(0.1, "npc"), 
          y = unit(0.55, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC1PC5 
PC1PC5_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC5
grid.text(PC_down(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
dev.off()

##############################
##############################


# piggies (all time points)

a <- "all"

# df for plots
DF_piggies_allguppied
head(DF_piggies_allguppied)

unique(simplified$sample_type)
unique(simplified$guppied_date)

# df for xml data 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="tALL") %>%
  group_split(component)


rbow <- rainbow(41,end=1, alpha=0.7)
legvec <- c(0,10,20,30,40)
# legend
color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.6, main, adj=c(0,0), cex=0.9)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.7)
}

DF_piggies_allguppied$collection_date = strptime(DF_piggies_allguppied$collection_date, format = "%Y-%m-%d") # convert to datetime objects
DF_piggies_allguppied$BIRTH_DAY = strptime(DF_piggies_allguppied$BIRTH_DAY, format = "%Y-%m-%d") # convert to datetime objects
unique(DF_piggies_allguppied$sample_type)
unique(DF_piggies_allguppied$guppied_date)

pdf(paste0(out_dir,"time_beta.pdf"))
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies_allguppied$PC1,DF_piggies_allguppied$PC2,
     xlab=paste0("PC1  (",get_var(find_PC1(xmldata)),"%)"),
     ylab=paste0("PC2  (",get_var(find_PC2(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies_allguppied$collection_date)-as.Date("2017-01-29 00:00:00")])
color_legend(min(DF_piggies_allguppied$PC1), max(DF_piggies_allguppied$PC2)-0.7,
             2.8, 0.5, "time post-weaning (days)", legvec, rbow)
mtext(paste0(PC_down(find_PC1(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_up(find_PC1(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_down(find_PC2(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_up(find_PC2(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="black", outer=TRUE)
# # PC3PC4
par(oma=c(6,6,6,6)) # all sides have 4 lines of space
par(mar=c(4,4,0.01,0.01))
plot(DF_piggies_allguppied$PC3,DF_piggies_allguppied$PC4,
     xlab=paste0("PC3  (",get_var(find_PC3(xmldata)),"%)"),
     ylab=paste0("PC4  (",get_var(find_PC4(xmldata)),"%)"),
     type="p",cex=0.8,cex.axis=0.6,cex.lab=0.6,
     col=rbow[as.Date(DF_piggies_allguppied$collection_date)-as.Date("2017-01-29 00:00:00")])
color_legend(min(DF_piggies_allguppied$PC3), max(DF_piggies_allguppied$PC4)-0.15,
             2.7, 0.35, "time post-weaning (days)", legvec, rbow)
mtext(paste0(PC_down(find_PC3(xmldata))), side=1, line=2, adj=0.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_up(find_PC3(xmldata))), side=1, line=2, adj=1.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_down(find_PC4(xmldata))), side=2, line=1, adj=0.0, cex=0.5, col="black", outer=TRUE)
mtext(paste0(PC_up(find_PC4(xmldata))), side=2, line=1, adj=1.0, cex=0.5, col="black", outer=TRUE)
dev.off()


pdf(paste0(out_dir,"time_beta_cohorts_PC1PC2.pdf"))
#################### Cohorts separately
par(oma=c(0,0,0,0)) # resetting the outer margins to default for the next plots
# control groups (Control, D-Scour, ColiGuard)
par(mfrow=c(3,2), mai = c(0.4, 0.4, 0.4, 0.4))
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="Control"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="Control"],
     main="PC1 PC2 Control",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="Control"])-as.Date("2017-01-30 00:00:00")])
color_legend(min(DF_piggies_allguppied$PC1), max(DF_piggies_allguppied$PC2)-0.6, 
             3.5, 0.9, "", legvec, rbow)
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="D-Scour"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="D-Scour"],
     main="PC1 PC2 D-Scour",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="D-Scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="ColiGuard"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="ColiGuard"],
     main="PC1 PC2 ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="ColiGuard"])-as.Date("2017-01-30 00:00:00")])
# Neo groups (Neo, Neo+D, Neo+C)
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="Neomycin"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="Neomycin"],
     main="PC1 PC2 Neomycin",
     xlab="PC1",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="Neomycin"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="Neomycin+D-Scour"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="Neomycin+D-Scour"],
     main="PC1 PC2 Neomycin+D-Scour",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="Neomycin+D-Scour"])-as.Date("2017-01-30 00:00:00")])
plot(DF_piggies_allguppied$PC1[DF_piggies_allguppied$Cohort=="Neomycin+ColiGuard"],
     DF_piggies_allguppied$PC2[DF_piggies_allguppied$Cohort=="Neomycin+ColiGuard"],
     main="PC1 PC2 Neomycin+ColiGuard",
     xlab="",ylab="",
     xlim=c(-3.5,3),
     ylim=c(0.8,4.7),
     cex.axis=0.8,
     type="p",col=rbow[as.Date(DF_piggies_allguppied$collection_date
                               [DF_piggies_allguppied$Cohort=="Neomycin+ColiGuard"])-as.Date("2017-01-30 00:00:00")])
dev.off()


# TIME BETA DENSITIES 

mytheme <- theme(legend.position="none",
                 axis.text.x=element_text(size=4),
                 axis.title.x=element_text(size=6),
                 axis.text.y=element_text(size=4),
                 axis.title.y=element_text(size=5))

unique(DF_piggies_allguppied$collection_date)
p1 <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC1, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC1),max(DF_piggies_allguppied$PC1))  +
  theme_bw()+
  mytheme +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))
g1.1 <- text_grob(paste0(PC_down(find_PC1(xmldata))),size=4,lineheight = 1)
g1.2 <- text_grob(paste0(PC_up(find_PC1(xmldata))),size=4,lineheight = 1)

p2 <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC2, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC2),max(DF_piggies_allguppied$PC2))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
g2.1 <- text_grob(paste0(PC_down(find_PC2(xmldata))),size=4,lineheight = 1)
g2.2 <- text_grob(paste0(PC_up(find_PC2(xmldata))),size=4,lineheight = 1)

p3 <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC3, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC3),max(DF_piggies_allguppied$PC3))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))
g3.1 <- text_grob(paste0(PC_down(find_PC3(xmldata))),size=4,lineheight = 1)
g3.2 <- text_grob(paste0(PC_up(find_PC3(xmldata))),size=4,lineheight = 1)

p4 <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC4, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC4),max(DF_piggies_allguppied$PC4))  +
  theme_bw()+
  mytheme+
  xlab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
g4.1 <- text_grob(paste0(PC_down(find_PC4(xmldata))),size=4,lineheight = 1)
g4.2 <- text_grob(paste0(PC_up(find_PC4(xmldata))),size=4,lineheight = 1)

p5 <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC5, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC5),max(DF_piggies_allguppied$PC5)) +
  theme_bw()+
  mytheme+
  xlab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))
g5.1 <- text_grob(paste0(PC_down(find_PC5(xmldata))),size=4,lineheight = 1)
g5.2 <- text_grob(paste0(PC_up(find_PC5(xmldata))),size=4,lineheight = 1)


for_legend_only <- DF_piggies_allguppied %>%
  # filter(collection_date=="t0"|
  #          collection_date=="t2"|
  #          collection_date=="t4"|
  #          collection_date=="t6"|
  #          collection_date=="t8"|
  #          collection_date=="t9") %>%
  ggplot(., aes(x = PC5, fill = as.character(collection_date))) + 
  geom_density(alpha = 0.5) +
  xlim(min(DF_piggies_allguppied$PC5),max(DF_piggies_allguppied$PC5)) +
  theme(legend.position="right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))
leg <- get_legend(for_legend_only)

lay <- rbind(c(1,1,1,1,2,2,2,2),
             c(1,1,1,1,2,2,2,2),
             c(6,6,7,7,8,8,9,9),
             c(3,3,3,3,4,4,4,4),
             c(3,3,3,3,4,4,4,4),
             c(10,10,11,11,12,12,13,13),
             c(5,5,5,5,16,16,16,16),
             c(5,5,5,5,16,16,16,16),
             c(14,14,15,15,16,16,16,16))

pdf(paste0(out_dir,"time_beta_densities.pdf"), width=7,height=5)
grid.arrange(p1,p2,p3,p4,p5,
             g1.1,g1.2,
             g2.1,g2.2,
             g3.1,g3.2,
             g4.1,g4.2,
             g5.1,g5.2,
             leg,
             layout_matrix = lay)
dev.off()


##############################
##############################
##############################


# piggies (guppied by time point)

df <- DF_piggies_time # dataframe for plots to be used 
a <- "all" # setting for xml data extraction (only sample_type necessary) 


# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

########## guppied by time - clustering by cohort #####################
pdf(paste0(out_dir,"piggies_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component) 
t0 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
####
t0 <- PC3PC4_plots$plots[[1]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component) 
t0 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC3PC4_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC3PC4_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC3PC4_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC3PC4_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC3PC4_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
dev.off()
###############################

########### groupA (guppied by time point) ####################

df <- groupA # dataframe for plots to be used 
a <- "groupA" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupB (guppied by time point) ####################

df <- groupB # dataframe for plots to be used 
a <- "groupB" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupC (guppied by time point) ####################

df <- groupC # dataframe for plots to be used 
a <- "groupC" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupD (guppied by time point) ####################

df <- groupD # dataframe for plots to be used 
a <- "groupD" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupE (guppied by time point) ####################

df <- groupE # dataframe for plots to be used 
a <- "groupE" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupF (guppied by time point) ####################

df <- groupF # dataframe for plots to be used 
a <- "groupF" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################

########### groupG (guppied by time point) ####################

df <- groupG # dataframe for plots to be used 
a <- "groupG" # setting for xml data extraction (only sample_type necessary) 
# re-order cohort 
df$Cohort <- factor(df$Cohort, 
                    levels=c("Control", 
                             "D-Scour", 
                             "ColiGuard",
                             "Neomycin",
                             "Neomycin+D-Scour",
                             "Neomycin+ColiGuard"
                    ))

PC1PC2_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC2, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC3PC4_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC3, y=PC4, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))
PC1PC5_plots <- df %>% 
  group_by(guppied_date) %>% 
  do(plots=ggplot(data=.) +
       aes(x=PC1, y=PC5, color=Cohort) + 
       geom_point(size=1.5) + 
       theme(legend.position="top",
             plot.margin=unit(c(0.2,0.2,2.9,2.9),"cm")) +
       ggtitle(unique(.$guppied_date))+
       stat_ellipse(inherit.aes = TRUE, level = 0.80))

pdf(paste0(out_dir,a,"_guppied_by_time.pdf"))
par(mar=c(4,4,0.01,0.01))
par(oma=c(6,6,6,6))
###
t0 <- PC1PC2_plots$plots[[1]]
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t0") %>%
  group_split(component)
t0 +
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"),
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"),
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t2 <- PC1PC2_plots$plots[[2]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t2") %>%
  group_split(component) 
t2 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t4 <- PC1PC2_plots$plots[[3]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t4") %>%
  group_split(component) 
t4 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t6 <- PC1PC2_plots$plots[[4]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t6") %>%
  group_split(component) 
t6 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t8 <- PC1PC2_plots$plots[[5]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t8") %>%
  group_split(component) 
t8 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###
t9 <- PC1PC2_plots$plots[[6]] 
xmldata <- simplified %>%
  filter(sample_type==a) %>%
  filter(guppied_date=="t9") %>%
  group_split(component) 
t9 + 
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.3, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.85, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 5))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.25, "npc"),
          gp = gpar(fontsize = 5))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 5))
###############################



###########################################################################################
###########################################################################################
###########################################################################################

# Checking significance and plot based on significance: 


###########################################################################################

###################################
###################################


# function 1. : takes dataframe, filters less than n observations, per component, per collection date, per cohort 

myfun_pvalue_df_prep1 <- function(df1_here,n) {
  # first, join the PCs to be in one columns, so you can easily apply the t.test function to one column only 
  a <- df1_here %>%
    select(PC1,PC2,PC3,PC4,PC5,Cohort,collection_date,guppied_date,sample_type) %>%
    pivot_longer(
      cols = 1:5,
      names_to = "component",
      values_to = "value",
      values_drop_na = FALSE
    )
  
  # look at distribution - filter out if less than n observations, per component, per collection date, per cohort 
  df1 <- setDT(a)[, .(Freq = .N), by = .(collection_date,Cohort,component)]
  df1 <- df1 %>% filter(!Freq<n)
  a <- merge(a,df1)
  a <- a %>% 
    select(Cohort,collection_date,guppied_date,sample_type,component,value)
  return(a)
}

###################################
###################################

# function 2. : takes dataframe, computes t.test of PC values between Cohorts (within same date,group)

myfun_pvalue_df_prep2 <- function(df2_here) {
  
  # create an empty df to build on
  df_pval <- data.frame(
    group_2 = character(),
    group_1 = character(),
    p_value = numeric(),
    which_PC = character(),
    which_colldate = character()
  )
  
  a <- df2_here
  
  # for each of the components (PC1 to PC5) ...
  listcompos <- unique(a$component)
  
  for (compo in listcompos) {
    df <- a %>% filter(component==compo)
    
    listcoldates <- unique(df$guppied_date)
    
    for (colldate in listcoldates) {
      
      z <- pairwise.t.test(df$value[df$guppied_date==colldate], df$Cohort[df$guppied_date==colldate], p.adjust = "none")$p.value
      z <- as.data.frame(z)
      z
      z$group_2 <- rownames(z)
      rownames(z) <- NULL
      z <- z %>%
        pivot_longer(
          cols = -group_2,
          names_to = "group_1",
          values_to = "p_value",
          values_drop_na = TRUE
        ) 
      
      z <- z %>% mutate(which_colldate = colldate,
                        which_PC = compo)
      df_pval <- rbind(df_pval,z) 
    }
  }
  return(df_pval)
}

###################################
###################################

# function 3. : takes p-values, removes unnecessary comparisons 

myfun_pvalues_filtering <- function(df_pval) {
  
  # Final adjustments 
  
  df_pval <- df_pval %>%
    # join the two groups for which the p-value has been computed
    mutate(comparison=paste0(group_1,"_vs_",group_2)) %>%
    select(p_value,which_colldate,which_PC,comparison,group_1,group_2)
  
  
  # filtering to keep only meaningful comparisons 
  # to be kept: 
  meaningfulcomparisons <- c("Control_vs_ColiGuard", "ColiGuard_vs_Control",
                             "Control_vs_D-Scour", "D-Scour_vs_Control",
                             "Control_vs_Neomycin", "Neomycin_vs_Control",
                             "Neomycin_vs_Neomycin+D-Scour", "Neomycin+D-Scour_vs_Neomycin",
                             "Neomycin_vs_Neomycin+ColiGuard", "Neomycin+ColiGuard_vs_Neomycin")
  
  # eliminate useless comparisons
  df_pval <- df_pval[df_pval$comparison %in% meaningfulcomparisons,]
  
  return(df_pval)
  
}

###################################
###################################

# binding all dataframes except the one where all samples (irrespective of dates) where guppied in one run 

# DF_piggies_time,groupA,groupB,groupC,groupD,groupE,groupF,groupG

all <- rbind(DF_piggies_time, 
             groupA, 
             groupB, 
             groupC, # neomycin vs control
             groupD, # neomycin vs neomycin+D-Scour
             groupE, # neomycin vs neomycin+ColiGuard
             groupF, # control vs D-scour
             groupG) # control vs ColiGuard

all$groupsplit <- paste0(all$sample_type,"_",all$guppied_date)

# splitting into multiple dataframes (by file name)
multi_DFs <- split( all , f = all$groupsplit )

# prep empty df to build on 
significant <- data.frame(
  group_1 = character(),
  group_2 = character(),
  p_value = numeric(),
  which_colldate = character(),
  which_PC = character(),
  comparison = character(),
  pval.adj = numeric(),
  groupsplit = character(),
  stringsAsFactors = FALSE
)


###################################
###################################


# runs single dataframe (one df : one guppy run) through functions; returns all p-values 


for (singl_DF in multi_DFs) {
  
  # function 1 
  a1 <- myfun_pvalue_df_prep1(singl_DF,2)
  
  # function 2
  a2 <- myfun_pvalue_df_prep2(a1)
  
  # function 3
  a3 <- myfun_pvalues_filtering(a2)
  
  # convert to class dataframe
  a4 <- as.data.frame(a3)
  
  # assign name of dataframe to dataframe (useful for later rbinding)
  a5 <- a4 %>% 
    mutate(groupsplit = singl_DF$groupsplit[1])
  
  # rbind all
  significant <- rbind(
    significant,
    a5)
  
}

###################################
###################################

significant$test <- "pairwise.t.test"
# save stats
# addWorksheet(wb, "guppy_pvalues")
# writeData(wb, sheet = "guppy_pvalues", significant, rowNames = FALSE)
fwrite(x=significant, file=paste0(stats_dir,"guppy_pvalues.csv"))

significant$test <- NULL
###################################
###################################


# Post-hoc correction: 
final <- significant %>% 
  group_by(groupsplit) %>%
  add_tally() %>%
  mutate(threshold = 0.05/n) %>%
  filter(p_value<0.05) # threshold
final
head(final)

###################################
###################################


final$padj_method <- "Bonferroni"

# save stats
# addWorksheet(wb, "guppy_padj")
# writeData(wb, sheet = "guppy_padj", final, rowNames = FALSE)
fwrite(x=final, file=paste0(stats_dir,"guppy_padj.csv"))

###################################
###################################

final <- cSplit(final, "groupsplit","_")
colnames(final)[colnames(final)=="groupsplit_1"] <- "dataframe"
colnames(final)[colnames(final)=="groupsplit_2"] <- "guppied_date"


# dummy df to associate guppied_date with collection_date
dummy <- data.frame(guppied_date = as.character(c("t0","t2","t4","t6","t8","t9")),
                    collection_date = as.character(c("t0","t2","t4","t6","t8","t9")))


unique(final$guppied_date)
# adding collection_date to dataframe
df <- inner_join(final,dummy) 
unique(df$collection_date)


# string replacement (as we haven't included DF_piggies (single guppy run), 
# for coherence we need to specify that this data comes from DF_piggies_time (multiple guppy runs))
df$dataframe <- gsub("piggies","DF_piggies_time",df$dataframe)



############

# XML data extract

# taking all along except the guppy run with all the samples from all time points 
unique(simplified$guppied_date)
simplified2 <- simplified %>%
  filter(!guppied_date == "tALL")  
simplified2 <- simplified2 %>%
  filter(!guppied_date == "tNONE") 
# now the only piggies are from the DF_piggies_time guppy runs
unique(simplified2$sample_type)
#simplified2$sample_type <- gsub("piggies","DF_piggies_time",simplified2$sample_type)
unique(simplified2$guppied_date)
unique(df$guppied_date)


###################################
###################################



mytheme <- theme(legend.position = "none",
                 axis.text.x=element_text(size=3),
                 axis.title.x=element_text(size=5),
                 axis.text.y=element_text(size=3),
                 axis.title.y=element_text(size=4),
                 plot.title = element_text(size = 6, face = "bold"))

# Plots only statistically significant observations,
# reporting xml data and dataframe (guppy run) it comes from

unique(df$collection_date)

mygrobs <- vector('list', nrow(df))
pdf(paste0(out_dir,"guppy_sign_cohorts.pdf"), onefile = TRUE)
for (A in rownames(df)) {
  A <- as.numeric(A) # dataframes rownames must be taken as numeric
  
  kk <- eval(as.name(paste(df$dataframe)))
  kk$collection_date <- gsub("2017-01-31","t0",kk$collection_date)
  kk$collection_date <- gsub("2017-02-07","t2",kk$collection_date)
  kk$collection_date <- gsub("2017-02-14","t4",kk$collection_date)
  kk$collection_date <- gsub("2017-02-21","t6",kk$collection_date)
  kk$collection_date <- gsub("2017-02-28","t8",kk$collection_date)
  kk$collection_date <- gsub("2017-03-03","t9",kk$collection_date)
  
  # subsetting of original dataframe based on what is statistically significant (rows of df)
  pp <- kk %>%
    filter(collection_date==as.character(df$collection_date[A])) %>%
    filter(sample_type==as.character(df$dataframe[A])) %>%
    filter(Cohort==as.character(df$group_1[A])|Cohort==as.character(df$group_2[A])) %>%
    select(PC1,PC2,PC3,PC4,PC5,Cohort,collection_date,guppied_date,sample_type) %>%
    pivot_longer(
      cols = 1:5,
      names_to = "component",
      values_to = "value",
      values_drop_na = FALSE
    ) %>%
    filter(component==as.name(paste(df$which_PC[A])))
  
  
  pp$Cohort <- factor(pp$Cohort, 
                      levels=c("Control", 
                               "D-Scour", 
                               "ColiGuard",
                               "Neomycin",
                               "Neomycin+D-Scour",
                               "Neomycin+ColiGuard"))
  
  # save some parameters to report on plot
  title1 <- unique(pp$guppied_date)
  title2 <- unique(pp$sample_type)
  PC_lab <- unique(pp$component)
  
  # add xml data 
  # subsetting the xml data based on the significant observations dataframe
  a <- df$dataframe[A]
  b <- df$guppied_date[A]
  c <- df$which_PC[A]
  xmldata <- simplified2 %>%
    filter(sample_type==a) %>%
    filter(guppied_date==b) %>%
    filter(component==c)
  
  # build plot 
  p <- ggplot(pp, aes(x=value, fill=Cohort)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'stack')+
    ggtitle(paste0(title1)) +    # this way it contains df info: ggtitle(paste0(title1,"_",title2))
    theme_bw()+
    mytheme+
    xlab(paste0(PC_lab," (",get_var(xmldata),"%)"))+
    scale_fill_discrete(drop=FALSE)
  
  g1 <- text_grob(paste0(PC_down(xmldata)),size=3,lineheight = 1)
  g2 <- text_grob(paste0(PC_up(xmldata)),size=3,lineheight = 1)
  
  lay <- rbind(c(1,1,1,1,1),
               c(1,1,1,1,1),
               c(2,2,2,3,3))
  
  mygrobs[[A]]  <- grid.arrange(p,g1,g2, layout_matrix = lay)
}
dev.off()

# re-order 
DF_piggies_time$Cohort <- factor(DF_piggies_time$Cohort, 
                            levels=c("Control", 
                                     "D-Scour", 
                                     "ColiGuard",
                                     "Neomycin",
                                     "Neomycin+D-Scour",
                                     "Neomycin+ColiGuard"))
DF_piggies_time <- DF_piggies[!is.na(DF_piggies_time$Cohort),]



# for legend only 
pp <- ggplot(DF_piggies_time, aes(x=PC1, fill=Cohort)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')+
  theme_bw()+
  mytheme+
  theme(legend.position="right",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4))+
  scale_fill_discrete(drop=FALSE)
leg <- get_legend(pp)

# a selection of plots from the figure generated above 
pdf(paste0(out_dir,"guppy_sign_cohorts_selection.pdf"))
lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,8))
grid.arrange(leg,
             mygrobs[[10]], # ctrl vs coliguard
             mygrobs[[11]], # ctrl vs coliguard
             mygrobs[[4]], # ctrl vs neo
             mygrobs[[5]], # neo vs neoD
             mygrobs[[8]], # neo vs neoD
             mygrobs[[9]], # neo vs neoD
             layout_matrix = lay)
dev.off()




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
complete <- read_csv(file = paste0(middle_dir,"guppy_xml_complete.df"))

############ VARIATION BEST EXPLAINED BY - DURING TIME - IN PIGGIES: 

# "complete" comes from guppy_XML_process.R
unique(complete$sample_type)
simplified2 <- complete %>%
  dplyr::select(sample_type, guppied_date,branch_width,var_explained,component,PC_position,taxa_simple) %>%
  # taking sample_type=="all": this is all the piglets 
  dplyr::filter(sample_type == "all") %>%
  # NOT taking guppied_date=="tALL" because I only want the data from individual guppy runs
  dplyr::filter(!guppied_date=="tALL")

# branch_width * var_explained = importance
unique(simplified2$sample_type)
unique(simplified2$guppied_date)
simplified2 <- simplified2 %>%
  dplyr::mutate(importance = as.numeric(branch_width)*var_explained) %>%
  select(guppied_date,taxa_simple,importance) 

both <- simplified2 %>% 
  group_by(taxa_simple,guppied_date) %>%
  dplyr::summarize(Sum_importance = sum(importance))

unique(both$taxa_simple) # 52 taxa
NROW(both)
both <- both[complete.cases(both), ] # remove NAs - necessary! 
NROW(both)


# keep only lineages described at family, genus and species level 
mytaxa <- both$taxa_simple
mytaxa <- gsub("_"," ",mytaxa)
mytaxa <- unique(mytaxa)
uids <- get_uid(mytaxa)
out <- classification(uids)
out2 <- do.call(rbind.data.frame, out)
NROW(out2)

unique(out2$rank)
out2 <- out2 %>% 
  dplyr::filter(rank=="order"|rank=="family"|rank=="genus"|rank=="species")
head(out2)

# change name to capital letters and replace space with _
out2$name <- toupper(out2$name)
out2$name <- gsub(" ","_",out2$name)

# store dataframe with ranking order
ranking_order <- out2

###

# subset dataframe based on list (in list we have species, genus and fam-level characterized taxa)
both2 <- both[both$taxa_simple %in% ranking_order$name,]

# normalize by date:
both3 <- both2 %>%
  group_by(guppied_date) %>%
  dplyr::mutate(Sum_importance=Sum_importance/sum(Sum_importance))
NROW(both3)
tail(both3)
head(both3)

both_wide <- both3 %>%
  pivot_wider(names_from = guppied_date, values_from = Sum_importance) %>%
  select(taxa_simple,t0,t2,t4,t6,t8,t9) %>%
  dplyr::mutate_all(~replace(., is.na(.), 0))

out5 <- as.data.frame(both_wide)

##############################

# to matrix conversion

rownames(out5) <- out5[,1]
out5[,1] <- NULL
out_m <- as.matrix(out5)

##############################

# normalization by date already done above. See here: 
colSums(out_m)

# remove cells with low rowSums
out_m <- as.matrix(out_m[rowSums(out_m)>0.01,])

time_beta_plot <- pheatmap(out_m, display_numbers = T, angle_col = 0,
                           cluster_rows = T, cluster_cols = F, fontsize_number = 6,
                           fontsize_row = 8)

pdf(paste0(out_dir,"time_beta_heatmap.pdf"))
time_beta_plot
dev.off()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



