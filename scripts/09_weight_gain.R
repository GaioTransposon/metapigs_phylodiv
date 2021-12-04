
library(magick)


# new weight measurements


source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/" # local 
out_dir_git =  "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 


###########################################################################################

# tiffs (timelines)
timeline<- image_read("/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/out/timeline.tiff")



# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == 'isolation_source'] <- 'pig'
mdat <- mdat %>%
  dplyr::select(Cohort,pig) %>%
  distinct()


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
head(weights)
head(weights_final)

weights <- rbind(weights,weights_final)
NROW(weights)
weights <- weights %>% distinct()
NROW(weights)

weights$sample <- paste0(weights$date,"_",weights$pig)


# load details (breed, line, bday, mothers)
details <- read_excel(paste0(source_data, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)
details$BIRTH_DAY <- as.character(details$BIRTH_DAY)
details$LINE <- as.character(details$LINE)
details <- details %>%
  dplyr::select(pig,BIRTH_DAY,LINE,breed,stig,nurse)


weights$pig <- as.character(weights$pig)

df <- inner_join(weights,mdat)
df <- inner_join(df,details)
df

df_weight <- df



############

# DELTAS of weight

##############################################################################

# functions for stats

myfun_weight_cohort_stats <- function(my_df, my_dates) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  # reorder
  x$Cohort.x <- factor(x$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))
  
  res1 <- aov(diff ~ Cohort.x, data=x)
  res <- TukeyHSD(res1)
  aov.out <- as.data.frame(res$Cohort.x)
  aov.out <- tibble::rownames_to_column(aov.out, "comparison")
  aov.out$time_delta=as.character(my_dates)
  aov.out$test <- "anova"
  aov.out$padj_method <- "TukeyHSD"
  
  return(aov.out)
}

myfun_weight_breed_stats <- function(my_df, my_dates) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  res1 <- aov(diff ~ breed.x, data=x)
  res <- TukeyHSD(res1)
  aov.out <- as.data.frame(res$breed.x)
  aov.out <- tibble::rownames_to_column(aov.out, "comparison")
  aov.out$time_delta=as.character(my_dates)
  aov.out$test <- "anova"
  aov.out$padj_method <- "TukeyHSD"
  
  return(aov.out)
}

myfun_weight_BIRTH_stats <- function(my_df, my_dates) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  res1 <- aov(diff ~ BIRTH_DAY.x, data=x)
  res <- TukeyHSD(res1)
  aov.out <- as.data.frame(res$BIRTH_DAY.x)
  aov.out <- tibble::rownames_to_column(aov.out, "comparison")
  aov.out$time_delta=as.character(my_dates)
  aov.out$test <- "anova"
  aov.out$padj_method <- "TukeyHSD"
  
  return(aov.out)
}

##############################################################################
round(2.00648E-05,5)
# settings for plots: 
your_font_size <- 2 
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

# functions for plots


myfun_weight_cohort_plot <- function(my_df) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  # reorder
  x$Cohort.x <- factor(x$Cohort.x, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))
  
  cw_summary <- x %>% 
    distinct() %>%
    group_by(Cohort.x) %>% 
    tally()
  
  my_plot <- ggboxplot(x, x="Cohort.x", y="diff", fill = "Cohort.x", 
                       ylab="weight change (%)", legend = "none")+
    My_Theme+
    geom_text(data = cw_summary,
              aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
    stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 
  
  return(my_plot)
}

myfun_weight_breed_plot <- function(my_df) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  cw_summary <- x %>% 
    distinct() %>%
    group_by(breed.x) %>% 
    tally()
  
  my_plot <- ggboxplot(x, x="breed.x", y="diff", fill = "breed.x", 
                       ylab="weight change (%)", legend = "none")+
    My_Theme+
    geom_text(data = cw_summary,
              aes(breed.x, Inf, label = n), vjust="inward", size = your_font_size)+
    stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 
  my_plot <- set_palette(my_plot, "jco")
  
  return(my_plot)
}

myfun_weight_BIRTH_plot <- function(my_df) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$weight.y-x$weight.x)/x$weight.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<100) %>% dplyr::filter(diff>-100)
  
  cw_summary <- x %>% 
    distinct() %>%
    group_by(BIRTH_DAY.x) %>% 
    tally()
  
  my_plot <- ggboxplot(x, x="BIRTH_DAY.x", y="diff", fill = "BIRTH_DAY.x", 
                       ylab="weight change (%)", legend = "none")+
    My_Theme+
    geom_text(data = cw_summary,
              aes(BIRTH_DAY.x, Inf, label = n), vjust="inward", size = your_font_size)+
    stat_compare_means(method = "anova", label.x=1.5, size = your_font_size)
  my_plot <- set_palette(my_plot, "jco")
  return(my_plot)
}


##############################################################################

# subsetting
piggies_t0 <- df_weight %>% dplyr::filter(date == "t0") 
piggies_t2 <- df_weight %>% dplyr::filter(date == "t2") 
piggies_t4 <- df_weight %>% dplyr::filter(date == "t4") 
piggies_t6 <- df_weight %>% dplyr::filter(date == "t6") 
piggies_t8 <- df_weight %>% dplyr::filter(date == "t8") 

###########################

df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("pig"))
df_t2_t4 <- merge(piggies_t2,piggies_t4, by=c("pig"))
df_t4_t6 <- merge(piggies_t4,piggies_t6, by=c("pig"))
df_t6_t8 <- merge(piggies_t6,piggies_t8, by=c("pig"))
df_t2_t6 <- merge(piggies_t2,piggies_t6, by=c("pig"))
df_t4_t8 <- merge(piggies_t4,piggies_t8, by=c("pig"))
df_t0_t8 <- merge(piggies_t0,piggies_t8, by=c("pig"))

###########################

# plot: 

myfun_weight_cohort_plot(df_t0_t2)
myfun_weight_cohort_plot(df_t2_t4)
myfun_weight_cohort_plot(df_t4_t6)
myfun_weight_cohort_plot(df_t6_t8)
myfun_weight_cohort_plot(df_t2_t6)
myfun_weight_cohort_plot(df_t4_t8)
myfun_weight_cohort_plot(df_t0_t8)
#
myfun_weight_breed_plot(df_t0_t2)
myfun_weight_breed_plot(df_t2_t4)
myfun_weight_breed_plot(df_t4_t6)
myfun_weight_breed_plot(df_t6_t8)
myfun_weight_breed_plot(df_t2_t6)
myfun_weight_breed_plot(df_t4_t8)
myfun_weight_breed_plot(df_t0_t8)
#
myfun_weight_BIRTH_plot(df_t0_t2)
myfun_weight_BIRTH_plot(df_t2_t4)
myfun_weight_BIRTH_plot(df_t4_t6)
myfun_weight_BIRTH_plot(df_t6_t8)
myfun_weight_BIRTH_plot(df_t2_t6)
myfun_weight_BIRTH_plot(df_t4_t8)
myfun_weight_BIRTH_plot(df_t0_t8)

###########################

# get and gather all stats: 
weight_cohort_stats <- rbind(myfun_weight_cohort_stats(df_t0_t2,"t0_t2"), 
                            myfun_weight_cohort_stats(df_t2_t4,"t2_t4"), 
                            myfun_weight_cohort_stats(df_t4_t6,"t4_t6"), 
                            myfun_weight_cohort_stats(df_t6_t8,"t6_t8"), 
                            myfun_weight_cohort_stats(df_t2_t6,"t2_t6"), 
                            myfun_weight_cohort_stats(df_t4_t8,"t4_t8"), 
                            myfun_weight_cohort_stats(df_t0_t8,"t0_t8")
)

weight_breed_stats <- rbind(myfun_weight_breed_stats(df_t0_t2,"t0_t2"), 
                             myfun_weight_breed_stats(df_t2_t4,"t2_t4"), 
                             myfun_weight_breed_stats(df_t4_t6,"t4_t6"), 
                             myfun_weight_breed_stats(df_t6_t8,"t6_t8"), 
                             myfun_weight_breed_stats(df_t2_t6,"t2_t6"), 
                             myfun_weight_breed_stats(df_t4_t8,"t4_t8"), 
                             myfun_weight_breed_stats(df_t0_t8,"t0_t8")
)

weight_BIRTH_stats <- rbind(myfun_weight_BIRTH_stats(df_t0_t2,"t0_t2"), 
                             myfun_weight_BIRTH_stats(df_t2_t4,"t2_t4"), 
                             myfun_weight_BIRTH_stats(df_t4_t6,"t4_t6"), 
                             myfun_weight_BIRTH_stats(df_t6_t8,"t6_t8"), 
                             myfun_weight_BIRTH_stats(df_t2_t6,"t2_t6"), 
                             myfun_weight_BIRTH_stats(df_t4_t8,"t4_t8"), 
                             myfun_weight_BIRTH_stats(df_t0_t8,"t0_t8")
)

############

# convert rownames to first column
weight_cohort_stats <- setDT(weight_cohort_stats, keep.rownames = TRUE)[]
# add data to workbook 
# addWorksheet(wb, "weight_cohort_stats")
# writeData(wb, sheet = "weight_cohort_stats", weight_cohort_stats, rowNames = FALSE)
fwrite(x=weight_cohort_stats, file=paste0(stats_dir,"weight_cohort_stats.csv"))


# convert rownames to first column
weight_breed_stats <- setDT(weight_breed_stats, keep.rownames = TRUE)[]
# add data to workbook 
# addWorksheet(wb, "weight_breed_stats")
# writeData(wb, sheet = "weight_breed_stats", weight_breed_stats, rowNames = FALSE)

# convert rownames to first column
weight_BIRTH_stats <- setDT(weight_BIRTH_stats, keep.rownames = TRUE)[]
# add data to workbook 
# addWorksheet(wb, "weight_BIRTH_stats")
# writeData(wb, sheet = "weight_BIRTH_stats", weight_BIRTH_stats, rowNames = FALSE)

############

stats_filtered_output <- weight_cohort_stats %>% dplyr::filter(`p adj` <= 0.05) %>%
  dplyr::filter(comparison=="Neomycin+D-Scour-Neomycin"| 
                  comparison=="Neomycin+ColiGuard-Neomycin"|
                  comparison=="Control-Neomycin"|
                  comparison=="Control-D-Scour"|
                  comparison=="Control-ColiGuard")

###########################

# this is for extracting the legend only
df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("pig"))
df_t0_t2$diff = ((df_t0_t2$weight.y-df_t0_t2$weight.x)/df_t0_t2$weight.y)*100
df_t0_t2$Cohort.x <- factor(df_t0_t2$Cohort.x, 
                            levels=c("Control", 
                                     "D-Scour", 
                                     "ColiGuard",
                                     "Neomycin",
                                     "Neomycin+D-Scour",
                                     "Neomycin+ColiGuard"))
for_legend_only <- ggboxplot(df_t0_t2, x = "Cohort.x", y = "diff", fill = "Cohort.x", 
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
  guides(fill=guide_legend("Cohort")) +
  My_Theme
addSmallLegend <- function(myPlot, pointSize = 10, textSize = 8, spaceLegend = 0.6) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
# Apply on original plot
for_legend_only <- addSmallLegend(for_legend_only)
leg <- get_legend(for_legend_only)

###########################

# final plot: 

empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(myfun_weight_cohort_plot(df_t0_t2),
                    myfun_weight_cohort_plot(df_t2_t4),
                    myfun_weight_cohort_plot(df_t4_t6),
                    myfun_weight_cohort_plot(df_t6_t8),
                    ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("t0-t2","t2-t4","t4-t6","t6-t8"),
                    label_size = 10)
bottom_row = plot_grid(myfun_weight_cohort_plot(df_t2_t6),
                       myfun_weight_cohort_plot(df_t4_t8),
                       myfun_weight_cohort_plot(df_t0_t8),
                       leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("t2-t6","t4-t8","t0-t8",""),
                       label_size = 10)
weight_cohort_plots <- plot_grid(empty_space,
                       top_row,
                       bottom_row,
                       nrow=3)

pdf(paste0(out_dir,"weight_deltas_by_cohort.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(weight_cohort_plots)
dev.off()

##############################################################################
##############################################################################

# BREED


###########################

# this is for extracting the legend only
df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("pig"))
df_t0_t2$diff = ((df_t0_t2$weight.y-df_t0_t2$weight.x)/df_t0_t2$weight.y)*100

for_legend_only <- ggboxplot(df_t0_t2, x = "breed.x", y = "diff", fill = "breed.x", 
                             legend = "right")+
  guides(fill=guide_legend("Cohort")) +
  My_Theme
for_legend_only <- set_palette(for_legend_only, "jco")
addSmallLegend <- function(myPlot, pointSize = 10, textSize = 8, spaceLegend = 0.6) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
# Apply on original plot
for_legend_only <- addSmallLegend(for_legend_only)
leg <- get_legend(for_legend_only)

###########################

# final plot: 

empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(myfun_weight_breed_plot(df_t0_t2),
                    myfun_weight_breed_plot(df_t2_t4),
                    myfun_weight_breed_plot(df_t4_t6),
                    myfun_weight_breed_plot(df_t6_t8),
                    ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("t0-t2","t2-t4","t4-t6","t6-t8"),
                    label_size = 10)
bottom_row = plot_grid(myfun_weight_breed_plot(df_t2_t6),
                       myfun_weight_breed_plot(df_t4_t8),
                       myfun_weight_breed_plot(df_t0_t8),
                       leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("t2-t6","t4-t8","t0-t8",""),
                       label_size = 10)
weight_breed_plots <- plot_grid(empty_space,
                                 top_row,
                                 bottom_row,
                                 nrow=3)

pdf(paste0(out_dir,"weight_deltas_by_breed.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(weight_breed_plots)
dev.off()

##############################################################################
##############################################################################


# BIRTH


###########################

# this is for extracting the legend only
df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("pig"))
df_t0_t2$diff = ((df_t0_t2$weight.y-df_t0_t2$weight.x)/df_t0_t2$weight.y)*100
for_legend_only <- ggboxplot(df_t0_t2, x = "BIRTH_DAY.x", y = "diff", fill = "BIRTH_DAY.x", 
                             legend = "right")+
  guides(fill=guide_legend("BIRTH")) +
  My_Theme
for_legend_only <- set_palette(for_legend_only, "jco")
addSmallLegend <- function(myPlot, pointSize = 10, textSize = 8, spaceLegend = 0.6) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
# Apply on original plot
for_legend_only <- addSmallLegend(for_legend_only)
leg <- get_legend(for_legend_only)

###########################

# final plot: 

empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(myfun_weight_BIRTH_plot(df_t0_t2),
                    myfun_weight_BIRTH_plot(df_t2_t4),
                    myfun_weight_BIRTH_plot(df_t4_t6),
                    myfun_weight_BIRTH_plot(df_t6_t8),
                    ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("t0-t2","t2-t4","t4-t6","t6-t8"),
                    label_size = 10)
bottom_row = plot_grid(myfun_weight_BIRTH_plot(df_t2_t6),
                       myfun_weight_BIRTH_plot(df_t4_t8),
                       myfun_weight_BIRTH_plot(df_t0_t8),
                       leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("t2-t6","t4-t8","t0-t8",""),
                       label_size = 10)
weight_BIRTH_plots <- plot_grid(empty_space,
                                top_row,
                                bottom_row,
                                nrow=3)

pdf(paste0(out_dir,"weight_deltas_by_BIRTH.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(weight_BIRTH_plots)
dev.off()


##################################################################################
##################################################################################


# # save stats in existing workbook
# saveWorkbook(wb, paste0(out_dir_git,"stats.xlsx"), overwrite=TRUE)

##################################################################################
##################################################################################
##################################################################################
##################################################################################


