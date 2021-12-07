

source_data = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/out/" # local 


# tiffs (timelines)
timeline <- image_read(paste0(out_dir,"timeline.tiff"))

############

# DELTAS of alpha diversity

##############################################################################

# functions for stats

myfun_unrooted_stats <- function(my_df, my_dates) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$unrooted_pd.y-x$unrooted_pd.x)/x$unrooted_pd.y)*100
  
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
  aov.out$variable="unrooted_PD"
  
  return(aov.out)
}

myfun_bwpd_stats <- function(my_df, my_dates) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$bwpd.y-x$bwpd.x)/x$bwpd.y)*100
  
  # remove crazy outliers 
  x <- x %>% dplyr::filter(diff<1000) %>% dplyr::filter(diff>-1000)
  
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
  aov.out$variable="BWPD"
  
  return(aov.out)
}

##############################################################################

# settings for plots: 
your_font_size <- 2 
My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(), 
  axis.title.y = element_text(size = 8),
  axis.text.y = element_text(size = 8)) 

# functions for plots

myfun_unrooted_plot <- function(my_df) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$unrooted_pd.y-x$unrooted_pd.x)/x$unrooted_pd.y)*100
  
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
                       ylab="unrooted_pd change (%)", legend = "none")+
    My_Theme+
    geom_text(data = cw_summary,
              aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
    stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 
  
  return(my_plot)
}

myfun_bwpd_plot <- function(my_df) {
  
  x <- my_df
  
  # calculate change
  x$diff = ((x$bwpd.y-x$bwpd.x)/x$bwpd.y)*100
  
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
                       ylab="bwpd change (%)", legend = "none")+
    My_Theme+
    geom_text(data = cw_summary,
              aes(Cohort.x, Inf, label = n), vjust="inward", size = your_font_size)+
    stat_compare_means(method = "anova", label.x=1.5, size = your_font_size) 
  
  return(my_plot)
}

##############################################################################

# remove pigs that were born on 2017-01-06 and 2017-01-07
# so that max pigs diff in age is 3 days (2017-01-08 and 2017-01-11)


# load details (cross_breed, line, bday, Sows)
details <- read_excel(paste0(source_data, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")

# format details
colnames(details)[colnames(details) == 'STIG'] <- 'pig'
details$pig <- gsub("G","",details$pig)
details$pig <- gsub("T","",details$pig)

details <- details %>%
  dplyr::select(pig,BIRTH_DAY)
unique(details$BIRTH_DAY)

keep <- details %>%
  dplyr::filter(BIRTH_DAY > "2017-01-08 00:00:00") 
unique(keep$BIRTH_DAY)
NROW(keep)
keep %>%
  group_by(BIRTH_DAY) %>%
  tally()
keep_ID <- unique(keep$pig)

boggo_selection <- boggo[boggo$isolation_source %in% keep_ID,]

test <- boggo_selection %>%
  group_by(Cohort, collection_date) %>%
  tally()
View(test)

# subsetting
piggies_t0 <- boggo_selection %>% dplyr::filter(collection_date == "t0") 
piggies_t2 <- boggo_selection %>% dplyr::filter(collection_date == "t2") 
piggies_t4 <- boggo_selection %>% dplyr::filter(collection_date == "t4") 
piggies_t6 <- boggo_selection %>% dplyr::filter(collection_date == "t6") 
piggies_t8 <- boggo_selection %>% dplyr::filter(collection_date == "t8") 

###########################

df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("isolation_source"))
df_t2_t4 <- merge(piggies_t2,piggies_t4, by=c("isolation_source"))
df_t4_t6 <- merge(piggies_t4,piggies_t6, by=c("isolation_source"))
df_t6_t8 <- merge(piggies_t6,piggies_t8, by=c("isolation_source"))
df_t2_t6 <- merge(piggies_t2,piggies_t6, by=c("isolation_source"))
df_t4_t8 <- merge(piggies_t4,piggies_t8, by=c("isolation_source"))
df_t0_t8 <- merge(piggies_t0,piggies_t8, by=c("isolation_source"))

###########################

# plot: 

myfun_unrooted_plot(df_t0_t2)
myfun_unrooted_plot(df_t2_t4)
myfun_unrooted_plot(df_t4_t6)
myfun_unrooted_plot(df_t6_t8)
myfun_unrooted_plot(df_t2_t6)
myfun_unrooted_plot(df_t4_t8)
myfun_unrooted_plot(df_t0_t8)

myfun_bwpd_plot(df_t0_t2)
myfun_bwpd_plot(df_t2_t4)
myfun_bwpd_plot(df_t4_t6)
myfun_bwpd_plot(df_t6_t8)
myfun_bwpd_plot(df_t2_t6)
myfun_bwpd_plot(df_t4_t8)
myfun_bwpd_plot(df_t0_t8)

###########################

# get and gather all stats: 
alpha_deltas_cohorts <- rbind(myfun_unrooted_stats(df_t0_t2,"t0_t2"), 
                              myfun_bwpd_stats(df_t0_t2,"t0_t2"),
                              myfun_unrooted_stats(df_t2_t4,"t2_t4"), 
                              myfun_bwpd_stats(df_t2_t4,"t2_t4"),
                              myfun_unrooted_stats(df_t4_t6,"t4_t6"), 
                              myfun_bwpd_stats(df_t4_t6,"t4_t6"),
                              myfun_unrooted_stats(df_t6_t8,"t6_t8"), 
                              myfun_bwpd_stats(df_t6_t8,"t6_t8"),
                              myfun_unrooted_stats(df_t2_t6,"t2_t6"), 
                              myfun_bwpd_stats(df_t2_t6,"t2_t6"),
                              myfun_unrooted_stats(df_t4_t8,"t4_t8"), 
                              myfun_bwpd_stats(df_t4_t8,"t4_t8"),
                              myfun_unrooted_stats(df_t0_t8,"t0_t8"), 
                              myfun_bwpd_stats(df_t0_t8,"t0_t8")
)

# # add data to workbook 
# addWorksheet(wb, "alpha_deltas_cohorts")
# writeData(wb, sheet = "alpha_deltas_cohorts", alpha_deltas_cohorts, rowNames = FALSE)

fwrite(x=alpha_deltas_cohorts, file=paste0(stats_dir,"alpha_deltas_cohorts.csv"))

############
stats_filtered_output <- alpha_deltas_cohorts %>% dplyr::filter(`p adj` <= 0.05) %>%
  dplyr::filter(comparison=="Neomycin+D-Scour-Neomycin"| 
                  comparison=="Neomycin+ColiGuard-Neomycin"|
                  comparison=="Control-Neomycin"|
                  comparison=="Control-D-Scour"|
                  comparison=="Control-ColiGuard")
NROW(stats_filtered_output)


###########################

# this is for extracting the legend only
df_t0_t2 <- merge(piggies_t0,piggies_t2, by=c("isolation_source"))
df_t0_t2$diff = ((df_t0_t2$bwpd.y-df_t0_t2$bwpd.x)/df_t0_t2$bwpd.y)*100
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

#unrooted pd
empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(myfun_unrooted_plot(df_t0_t2),
                    myfun_unrooted_plot(df_t2_t4),
                    myfun_unrooted_plot(df_t4_t6),
                    myfun_unrooted_plot(df_t6_t8),
                    ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("t0-t2","t2-t4","t4-t6","t6-t8"),
                    label_size = 10)
bottom_row = plot_grid(myfun_unrooted_plot(df_t2_t6),
                       myfun_unrooted_plot(df_t4_t8),
                       myfun_unrooted_plot(df_t0_t8),
                       leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("t2-t6","t4-t8","t0-t8",""),
                       label_size = 10)
all_plots_unroo <- plot_grid(empty_space,
                             top_row,
                             bottom_row,
                             nrow=3)


#BWPD
empty_space = plot_grid(NULL, NULL, NULL, NULL, ncol=4)
top_row = plot_grid(myfun_bwpd_plot(df_t0_t2),
                    myfun_bwpd_plot(df_t2_t4),
                    myfun_bwpd_plot(df_t4_t6),
                    myfun_bwpd_plot(df_t6_t8),
                    ncol=4, 
                    rel_widths=c(0.25,0.25,0.25,0.25),
                    labels=c("t0-t2","t2-t4","t4-t6","t6-t8"),
                    label_size = 10)
bottom_row = plot_grid(myfun_bwpd_plot(df_t2_t6),
                       myfun_bwpd_plot(df_t4_t8),
                       myfun_bwpd_plot(df_t0_t8),
                       leg, ncol=4, 
                       rel_widths=c(0.25,0.25,0.25,0.25),
                       labels=c("t2-t6","t4-t8","t0-t8",""),
                       label_size = 10)
all_plots_bwpd <- plot_grid(empty_space,
                            top_row,
                            bottom_row,
                            nrow=3)

pdf(paste0(out_dir,"alpha_deltas_by_cohort.pdf"))
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(all_plots_unroo)
ggdraw() +
  draw_image(timeline, x = 0, y = 0.33) +
  draw_plot(all_plots_bwpd)
dev.off()


