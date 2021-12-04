library(readr)
library(tidyverse)
library(ggplot2)
library(splitstackshape)
library(ggpubr)


source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/out/" # local 


readcounts_samples <- read.csv(paste0(source_data,"readcounts_samples.csv"), header=FALSE)
NROW(readcounts_samples)

r <- cSplit(readcounts_samples, "V1","_")
r$DNA_plate <- paste0(r$V1_1,"_",r$V1_2)
r <- cSplit(r, "V1_3",".")

r <- r %>%
  dplyr::select(DNA_plate,V1_3_1, V3)
colnames(r) <- c("DNA_plate","DNA_well","read_count")
head(r)

alpha_1M <- read.csv(paste0(middle_dir,"all.alphadiv"), sep="")
alpha_100K <- read.csv(paste0(middle_dir,"PD_100k/all.alphadiv"), sep="")

get_cleaner_fun <- function(df) {
  
  colnames(df)[1] <- "provv"
  df <- cSplit(df, "provv","_")
  df$DNA_plate <- paste0(df$provv_1,"_",df$provv_2)
  df$DNA_well <- df$provv_3
  df <- df %>% dplyr::select(!starts_with("provv"))
  df_clean <- df
  return(df_clean)
  
}


alpha_100K <- get_cleaner_fun(alpha_100K)
alpha_1M <- get_cleaner_fun(alpha_1M)

alpha_100K <- alpha_100K %>%
  dplyr::select(unrooted_pd,bwpd,DNA_plate,DNA_well) %>%
  pivot_longer(cols=c(unrooted_pd,bwpd), names_to="metric") %>%
  dplyr::mutate(analysis = "100K")


alpha_1M <- alpha_1M %>%
  dplyr::select(unrooted_pd,bwpd,DNA_plate,DNA_well) %>%
  pivot_longer(cols=c(unrooted_pd,bwpd), names_to="metric") %>%
  dplyr::mutate(analysis = "1M")

alpha <- rbind(alpha_100K,
      alpha_1M) 

# median and sd 
alpha %>%
  group_by(metric) %>%
  dplyr::summarise(median_unrooted=median(unrooted_pd),
                   sd_unrooted=sd(unrooted_pd),
                   median_bwpd=median(bwpd),
                   sd_bwpd=sd(bwpd))


beta_100K <- read_csv(paste0(middle_dir,"PD_100K/pca_piggies_sel.txt.proj"), 
                    col_names = FALSE)
beta_1M <- read_csv(paste0(middle_dir,"pca_piggies_sel.txt.proj"), 
                      col_names = FALSE)

beta_100K <- get_cleaner_fun(beta_100K)
beta_1M <- get_cleaner_fun(beta_1M)

beta_100K <- beta_100K %>%
  dplyr::select(X2,X3,X4,X5,X6,DNA_plate,DNA_well) %>%
  pivot_longer(cols=c(X2,X3,X4,X5,X6), names_to="metric") %>%
  dplyr::mutate(analysis = "100K")

beta_1M <- beta_1M %>%
  dplyr::select(X2,X3,X4,X5,X6,DNA_plate,DNA_well) %>%
  pivot_longer(cols=c(X2,X3,X4,X5,X6), names_to="metric") %>%
  dplyr::mutate(analysis = "1M")

beta <- rbind(beta_100K,
               beta_1M) 

beta$metric <- gsub("X2","PC1",beta$metric)
beta$metric <- gsub("X3","PC2",beta$metric)
beta$metric <- gsub("X4","PC3",beta$metric)
beta$metric <- gsub("X5","PC4",beta$metric)
beta$metric <- gsub("X6","PC5",beta$metric)


both <- rbind(alpha,beta)
# give some order 
both$metric <- factor(both$metric, 
                                  levels=c("unrooted_pd","bwpd",
                                           "PC1","PC2","PC3","PC4","PC5"))


pdf(paste0(out_dir,"correlation_100k_vs_1M_analysis.pdf"))
both %>%
  pivot_wider(values_from = value, names_from = analysis) %>%
  ggplot(., aes(x=`100K`,y=`1M`))+
  geom_point(shape=1, size=1)+
  geom_smooth(method="lm", size=0.3)+
  facet_wrap(~metric, scales="free")+
  stat_cor(method = "pearson", r.accuracy = 0.0001, size=3)
dev.off()
