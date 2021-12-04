
###########################################################################################

library(readr)
library(openxlsx)


stats_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/stats/" # git 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/out/" # git 


## 1. create workbook 
wb <- createWorkbook()


## 2. open csv and save as sheets of the workbook
filenames <- list.files(stats_dir, pattern="*.csv", full.names=TRUE)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  clean_name <- gsub(".*stats//\\s*|.csv.*", "", as.character(each_filename))
  
  addWorksheet(wb, clean_name)
  writeData(wb, sheet = clean_name, ldf, colNames = TRUE)
  
}


## 3. save workbook 
saveWorkbook(wb, paste0(out_dir_git,"stats.xlsx"), overwrite=TRUE)


###########################################################################################



