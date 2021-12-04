
library(readr)
library(dplyr)
library(readxl)
library(splitstackshape)
library(data.table)

#install.packages("robCompositions")
#library(robCompositions)

library(tidyr)
library(tidyverse)
library(ggbiplot)
library(magrittr)
library(ggpubr)
library(grDevices)

#install.packages("colorRamps")
library(colorRamps)

library(EnvStats)

#install.packages("corrplot")
library(corrplot)

library(grid)
library(cowplot)

#install.packages("factoextra")
library(factoextra)

library(broom)
library(openxlsx)


source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/sortmerna/" # local 


###########################################################################################


# upload sortmerna output (from the original output we retained columns:plate_well, 16Sgene ID, e-value)
# and we filtered out all the 16S rRNA genes below e-30 threshold 
so <- read_table2(paste0(out_dir,"sortmeall_evaluefiltered.tsv"), col_names = FALSE)
so$X1 <- gsub("_S","", so$X1)
so <- so[,1:2]

#number of unique 16S rRNA genes in dataset (filtered with e-value < e-30)
NROW(unique(so$X2))

NROW(so)

so_work <- so

######################################################################

# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)


colnames(mdat)

colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'

mdat <- mdat %>%
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source)

######################################################################

# upload 16S names 

# 1. unzips the file (otherwise too large for github) and places it in out_dir
zipF<- paste0(source_data,"silva-bac-16s-id90_accession_taxonomy.txt.zip")
unzip(zipF,exdir=out_dir)

# 2. now it reads it from out_dir
silva <- read_table2(paste0(out_dir,"silva-bac-16s-id90_accession_taxonomy.txt"), 
                     col_names = FALSE)


######################################################################

# upload weight info 

weights <- read.csv(paste0(source_data,"weights.csv"))
colnames(weights) <- c("room","pen","pig","t0","t2","t4","t6","t8")

weights_final <- read.csv(paste0(source_data,"weights_final.csv"))
colnames(weights_final) <- c("room","pen","pig","date","weight")
weights_final$date <- gsub("6-Mar","t10",weights_final$date)
weights_final$date <- gsub("7-Mar","t10",weights_final$date)
weights_final$date <- gsub("8-Mar","t10",weights_final$date)
weights_final$date <- gsub("9-Mar","t10",weights_final$date)
weights_final <- weights_final %>%
  dplyr::select(pig,date,weight) %>%
  filter(!date=="10-Mar") # as it's NA

weights <- weights %>%
  dplyr::select(pig,t0,t2,t4,t6,t8) %>%
  pivot_longer(., cols = c(t0,t2,t4,t6,t8), names_to = "date", values_to = "weight")
weights <- as.data.frame(weights)

weights <- rbind(weights,weights_final)
NROW(weights)

weights$sample <- paste0(weights$date,"_",weights$pig)

######################################################################


# further reducing size: AIM: 

so_work$count <- as.numeric(paste0(1))

NROW(so_work)
so_done <- setDT(so_work)[,.(A = sum(count)), by = 'X1,X2']
NROW(so_done)
# much an improvement! 

#############

# once size of sortme file is reduced, we can parse it: 
so_done_temp <- cSplit(so_done,"X1","_")
so_done_temp$DNA_plate <- paste0(so_done_temp$X1_1,"_",so_done_temp$X1_2)
so_done_temp$DNA_well <- paste0(so_done_temp$X1_3)
head(so_done_temp)

so_done_temp <- so_done_temp %>%
  dplyr::select(DNA_plate,DNA_well,X2,A)

colnames(so_done_temp)[colnames(so_done_temp) == 'X2'] <- 'rRNA16S'
colnames(so_done_temp)[colnames(so_done_temp) == 'A'] <- 'count'
head(so_done_temp)

###############################

# time to merge to metadata!

NROW(so_done_temp)
NROW(mdat)

df <- left_join(so_done_temp,mdat)

NROW(df)
head(df)


# tM <- "2017-01-30"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-01-30",
  replacement = "tM",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-01-31",
  replacement = "t0",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-01",
  replacement = "t0",
  fixed = TRUE)

# t1 <- "2017-02-03"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-03",
  replacement = "t1",
  fixed = TRUE)

# t2 <- "2017-02-06" "2017-02-07" "2017-02-08"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-06",
  replacement = "t2",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-07",
  replacement = "t2",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-08",
  replacement = "t2",
  fixed = TRUE)

# t3 <- "2017-02-10"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-10",
  replacement = "t3",
  fixed = TRUE)

# t4 <- "2017-02-14"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-14",
  replacement = "t4",
  fixed = TRUE)

# t4 <- "2017-02-16" "2017-02-17"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-16",
  replacement = "t5",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-17",
  replacement = "t5",
  fixed = TRUE)

# t6 <- "2017-02-21"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-21",
  replacement = "t6",
  fixed = TRUE)

# t7 <- "2017-02-24"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-24",
  replacement = "t7",
  fixed = TRUE)

# t8 <- "2017-02-28"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-28",
  replacement = "t8",
  fixed = TRUE)

# t9 <- "2017-03-03"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-03",
  replacement = "t9",
  fixed = TRUE)

# t10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-06",
  replacement = "t10",
  fixed = TRUE)

df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-07",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-08",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-09",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-10",
  replacement = "t10",
  fixed = TRUE)

df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-08-14",
  replacement = "no-t-pos",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2018-01-24",
  replacement = "no-t-pos",
  fixed = TRUE)

# no-t-neg for negative control
df <- df %>%
  dplyr::mutate(collection_date = if_else(is.na(collection_date), "no-t-neg", collection_date))

unique(df$collection_date)

df$sample <- paste0(df$collection_date,"_",df$isolation_source)


###########

# parsing silva and joining info to my df: 

silva <- cSplit(silva,"X3",";")
colnames(silva) <- c("rRNA16S","access","King","Phy","Class","Order","Fam","Species")
silva_edit <- silva

# replace non-taxa with NA
silva_edit[ silva_edit == "uncultured" ] <- NA
silva_edit[ silva_edit == "Family" ] <- NA
silva_edit[ silva_edit == "Order" ] <- NA
silva_edit[ silva_edit == "Subgroup" ] <- NA
silva_edit[ silva_edit == "Incertae" ] <- NA
#
# extract first non-NA taxonomic definition --> to  rightmost column 
silva_edit_rightmost <- as.matrix(silva_edit)[cbind(1:nrow(silva_edit), 
                                          max.col(!is.na(silva_edit), "last"))] 
silva_new <- cbind(silva_edit,silva_edit_rightmost)
#

silva_new <- silva_new %>% dplyr::select(rRNA16S,Species)

# joining 
new <- inner_join(df,silva_new) 

NROW(new)
new <- na.omit(new)
new$rRNA16S_full <- paste0(new$rRNA16S,"_",new$Species)

###########

# normalization for library size
df1 <- new %>%
  dplyr::select(sample,rRNA16S_full,count) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(norm_count=count/sum(count))
NROW(df1)
head(df1)

# sum all the norm values that fall within same pig,date,taxa_2
df2 <- df1 %>%
  dplyr::group_by(sample,rRNA16S_full) %>%
  dplyr::summarise(last_count = sum(norm_count))
head(df2)

# compute correlation of weight with 16S species abundance at specific time points 

df3 <- df2

df3 <- cSplit(df3,splitCols = "rRNA16S_full",sep = "_") 

df3 <- na.omit(df3)

df4 <- df3 %>%
  dplyr::group_by(sample,rRNA16S_full_2) %>%
  dplyr::summarise(last_count = sum(last_count))
head(df4)


df4 <- as.data.frame(df4)

final <- left_join(df4,weights)
# omit where observations where no weight value is present
final <- na.omit(final)
final <- final %>% dplyr::select(sample,weight,last_count,rRNA16S_full_2)
final <- cSplit(final,splitCols = "sample",sep = "_") 
colnames(final) <- c("weight","last_count","species","date","pig")
head(final)
NROW(final)

unique(final$date)

# Spearman: 

# correlation weight with species abundance: stats:
myfun_findcorr <- function(df,timepoint) {
  
  df <- df %>% dplyr::filter(date==timepoint)
  
  # empty df
  all_rho <- data.frame(estimate=character(),
                        statistic=character(),
                        p.value=double(),
                        method=character(),
                        alternative=character(),
                        rRNA_species=character())
  
  # split df by species
  multiple_DFs <- split( df , f = df$species ,drop = TRUE)
  NROW(multiple_DFs)
  
  
  for (single_DF in multiple_DFs) {
    
    single <- as.data.frame(single_DF)
    
    if (NROW(single)>4) {
      
      species <- single$species[1]
      
      y <- cor.test(single$last_count,single$weight,method="sp")
      y <- tidy(y)
      y <- as.data.frame(y)
      y$species=species
      all_rho <- rbind(all_rho,y)
      
    }
    
    else (print("not enough observations"))
  }
  
  sub_all_rho <- all_rho %>% filter(p.value<0.05)
  sub_all_rho$date<-timepoint
  return(sub_all_rho)
}


# rbind all the results 
all_corr <- rbind(myfun_findcorr(final,"t0"),
                  myfun_findcorr(final,"t2"),
                  myfun_findcorr(final,"t4"),
                  myfun_findcorr(final,"t6"),
                  myfun_findcorr(final,"t8"),
                  myfun_findcorr(final,"t10"))

#####
# save the results
all_corr$taxa_origin <- "SortMeRNA"

# apply various p-adj methods 
all_corr$padj_hommel <- round(p.adjust(all_corr$p.value, "BH"), 4)
all_corr$padj_BY <- round(p.adjust(all_corr$p.value, "BY"), 4)
all_corr$padj_bonferroni <- round(p.adjust(all_corr$p.value, "bonferroni"), 4)

# addWorksheet(wb, "weight_taxa")
# writeData(wb, sheet = "weight_taxa", all_corr, rowNames = FALSE)
fwrite(x=all_corr, file=paste0(stats_dir,"weight_taxa.csv"))

#####


# create lists of species that at each time point are correlating with weight
mylist_t0 <- all_corr %>% dplyr::filter(date=="t0") %>% dplyr::select(species)
mylist_t2 <- all_corr %>% dplyr::filter(date=="t2") %>% dplyr::select(species)
mylist_t4 <- all_corr %>% dplyr::filter(date=="t4") %>% dplyr::select(species)
mylist_t6 <- all_corr %>% dplyr::filter(date=="t6") %>% dplyr::select(species)
mylist_t8 <- all_corr %>% dplyr::filter(date=="t8") %>% dplyr::select(species)
mylist_t10 <- all_corr %>% dplyr::filter(date=="t10") %>% dplyr::select(species)



# function to make boxplots
make_plot <- function(df,species_list,date) {
  
  timepoint <- as.character(date)
  # now plot as whisker plots only species that appeared to be correlatd with weight group
  toplot <-inner_join(df,species_list) %>% dplyr::filter(date==timepoint)
  # add weight groups to dataframe 
  toplot$weight_group <- cut(toplot$weight, 4)   # bins : equal size groups by weight
  print(ggplot(toplot,aes(x=weight_group,y=log(last_count),fill=weight_group))+
          geom_boxplot()+
          labs(y="log (norm. abundance)")+
          theme(legend.position="top",
                axis.text.x=element_blank(),
                strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
          stat_n_text(size = 3,angle=90,hjust=0)+
          facet_grid(cols = vars(species))+
          ggtitle(timepoint))
}



# making plots : for a couple of timepoints, subsetting to fit nicely

chunk <- 13
n <- nrow(mylist_t0)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
d <- split(mylist_t0,r)
plot_t0_1 <- make_plot(final,d$`1`,"t0")
plot_t0_2 <- make_plot(final,d$`2`,"t0")
plot_t0 <- ggarrange(plot_t0_1,plot_t0_2,nrow=2)


plot_t2 <- make_plot(final,mylist_t2,"t2")


n <- nrow(mylist_t4)
chunk <- 10
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
d <- split(mylist_t4,r)
plot_t4_1 <- make_plot(final,d$`1`,"t4")
plot_t4_2 <- make_plot(final,d$`2`,"t4")
plot_t4 <- ggarrange(plot_t4_1,plot_t4_2,nrow=2)


plot_t6 <- make_plot(final,mylist_t6,"t6")


plot_t8 <- make_plot(final,mylist_t8,"t8")


plot_t10 <- make_plot(final,mylist_t10,"t10")


pdf(paste0(out_dir,"sortmerna_weight.pdf"))
plot_t0
plot_t2
plot_t4
plot_t6
plot_t8
plot_t10
dev.off()




# PLOT PCA 

# filter out pos and neg controls (they add unwanted variance and we want to look at just piglet samples)
df2_5 <- df2 %>%
  dplyr::filter(., !grepl("no-t",sample))  %>%
  dplyr::filter(., !grepl("tM",sample)) 

df2_5 <- as.data.frame(df2_5)

# long to wide format
df3 <- df2_5 %>%
  pivot_wider(names_from = rRNA16S_full, values_from = last_count, values_fill = list(last_count = 0)) 
head(df3)

x <- as.data.frame(df3)
rownames(x) <- x$sample
x$sample <- NULL

# order left to right in descending order 
x <- x[,names(sort(colSums(x), decreasing = TRUE))]

rownames(x)



#############
# PCA of 20 most abundant rRNA16S genes 
mtcars.pca2 <- prcomp(x[,1:20], center = TRUE,scale. = TRUE)
dates <- substr(rownames(x[,1:20]), start = 1, stop = 3)  %<>%
  gsub('_$', '', .) # removes the last _
# reorder dates 
dates  = factor(dates, levels=c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"))



# Interesting things of the PCA
var <- get_pca_var(mtcars.pca2)
# Quality of representation
corrplot(var$cos2, is.corr=FALSE)
# Contributions of variables to PCs
corrplot(var$contrib, is.corr=FALSE)


# plot! 
my20 <- fviz_pca_ind(mtcars.pca2, 
                     geom.ind="point",
                     #fill.ind = dates, #col.ind = rainbow(n = 11),
                     pointshape = 21, pointsize = 2,
                     habillage = dates,
                     #geom.ind = "point", # show points only (nbut not "text") 
                     col.ind = dates, # color by groups
                     #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                     addEllipses = FALSE, # Concentration ellipses
                     title="")+
  scale_color_manual(name="time point", 
                     values=rainbow(n = 11))+
  theme(legend.position="none")+
  guides(color = guide_legend(nrow = 1))

# plot contributions
contr <- fviz_pca_var(mtcars.pca2, col.var = "contrib",
                      gradient.cols = c("#999999","#666666","#333333"),
                      repel = TRUE,
                      labelsize=3)+ # Avoid text overlapping
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())+
  ggtitle("")

# together
together <- fviz_pca_biplot(mtcars.pca2,
                geom.ind="point",
                pointshape = 21, pointsize = 2,
                habillage = dates,
                alpha.var ="contrib",
                repel = TRUE,labelsize=4) +
  scale_color_manual(name="time point",
                     values=rainbow(n = 11))+
  ggtitle("")+
  theme(legend.position="right",
        panel.border = element_rect(colour = "black", fill=NA, size=1))

pdf(paste0(out_dir,"sortmerna_20.pdf"))
together
dev.off()



# this is how many taxa correlated with weight 

head(all_corr)

z <- all_corr %>% 
  filter(p.value<=0.05) 

n_occur <- data.frame(table(z$species))
n_occur[n_occur$Freq > 1,] # tells you which ids occurred more than once.
these <- z[z$species %in% n_occur$Var1[n_occur$Freq > 1],]
these %>% arrange(species)


# have a look at the one species still significantly correlating after Bonferroni correction. 
final %>%
  dplyr::filter(species=="Coprothermobacter") %>%
  dplyr::filter(date=="t0") %>%
  ggplot(., aes(x=weight,y=log(last_count)))+
  geom_point()+
  geom_smooth(method="lm")

##########################

head(df)

NROW(unique(df$rRNA16S))

df$sample <- paste0(df$DNA_plate,"_",df$DNA_well)

moms <- df %>% dplyr::filter(Cohort=="Mothers") 
  
piggies <- df %>% dplyr::filter(Cohort=="Control"|
                       Cohort=="D-Scour"|
                       Cohort=="ColiGuard"|
                       Cohort=="Neomycin"|
                       Cohort=="Neomycin+D-Scour"|
                       Cohort=="Neomycin+ColiGuard") 

NROW(unique(moms$rRNA16S))
NROW(unique(piggies$rRNA16S))

moms <- moms %>% dplyr::select(sample) %>% distinct()
piggies <- piggies %>% dplyr::select(sample) %>% distinct()

# save without heading
write_csv(x=moms,file=paste0(middle_dir,"moms_list.txt"), col_names = FALSE)
write_csv(x=piggies,file=paste0(middle_dir,"piggies_list.txt"), col_names = FALSE)

NROW(moms)
NROW(piggies)



###########################################################################################
