######################################################################################################

library(readxl)
library(data.table)
library(readr)
library(splitstackshape)
library(dplyr)
library(tidyr)


# rationale: 

##### previous steps: 

# 1 # guppy fat is run 

# 2 # .xml conversion to .txt:

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# download forester.jar
# move forester.jar to HPC this way: rsync forester_1050.jar u152324@138.25.37.51:
# on HPC run: 
# java -cp ~/forester_1050.jar org.forester.application.phyloxml_converter -f=dummy some_input.xml some_output.txt

# in a loop: 
#for fpath in /shared/homes/s1/pig_microbiome/phy_10M_202109/PS_temp/*.xml
#do java -cp ~/forester_1050.jar org.forester.application.phyloxml_converter -f=dummy "$fpath" "$(basename "$fpath").txt"
#done

# move content of forester_out to local ~/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output

##### HERE : 
# 1 # .xml files are read in and parsed
# 2 # output can be used in guppy_plots.R

######################################################################################################

middle_dir = "/Users/danielagaio/Gaio/github/metapigs_phylodiv/middle_dir/" # git 
guppyout_dir = "/Users/danielagaio/Desktop/metapigs_phylodiv/phylosift/guppy/guppy_output" # local 

######################################################################################################


my.files = list.files(guppyout_dir,pattern=".gz.xml.txt")
my.files <- my.files[grep("fat_plate", my.files)]  


# construct an empty dataframe to build on 
complete.df <- data.frame(
  taxa = character(),
  branch_length = character(),
  branch_width = character(),
  file = character()
)


for (textfile in my.files) {
  
  # read in file 
  my.df <- read_csv(file.path(guppyout_dir,textfile), col_names = FALSE)
  
  # extract file name 
  myfilename <- basename(textfile)
  myfilename <- gsub("_R1_001.fastq.gz.xml.txt","",myfilename)

  my.df  <- my.df %>% dplyr::filter(!grepl("clade", X1))
  my.df  <- my.df %>% dplyr::filter(!grepl("color", X1))
  my.df  <- my.df %>% dplyr::filter(!grepl("red", X1))
  my.df  <- my.df %>% dplyr::filter(!grepl("green", X1))
  my.df  <- my.df %>% dplyr::filter(!grepl("blue", X1))
  
  # start the grepping of useful tree info(we only care about the hits containing "width")
  w <- grep("width", my.df$X1)
  widths <- my.df[w,]
  # length and name (row numbers)
  l <- w-1
  lengths <- my.df[l,]
  n <- w-2
  names <- my.df[n,]
  
  new_df <- cbind(names,lengths,widths)
  colnames(new_df) <- c("n","l","w")
  
  # clean
  new_df$n <- sub("<name>","",new_df$n)
  new_df$n <- sub("</name>","",new_df$n)
  new_df$l <- sub("<branch_length>","",new_df$l)
  new_df$l <- sub("</branch_length>","",new_df$l)
  new_df$w <- sub("<width>","",new_df$w)
  new_df$w <- sub("</width>","",new_df$w)

  # save file name as an extra column
  new_df$file <- myfilename
  
  colnames(new_df) <- c("taxa","branch_length","branch_width","file")
  
  # bind all dfs from all files 
  complete.df <- rbind(
    complete.df, 
    new_df
  )
  
}



complete.df <- cSplit(complete.df, "file","_")
complete.df$DNA_plate <- paste0(complete.df$file_2,"_",complete.df$file_3)
complete.df$DNA_well <- complete.df$file_4
complete.df$name <- complete.df$taxa


complete.df <- complete.df %>%
  dplyr::select(name,branch_length,branch_width,DNA_plate,DNA_well)

simplified.df <- complete.df
length(unique(simplified.df$name))

simplified.df$name <- simplified.df$name %>%  
  gsub('\\{|\\}', '', .) %>% # removes curly brackets
  gsub('\\[|\\]', '', .) %>% # removes square brackets
  gsub('[0-9]+', '', .)  %>% # removes digits
  gsub('^_', '', .) %>% # removes the first _
  gsub('_$', '', .) %>%  # removes the last _
  gsub('_$', '', .) # removes the last _

# remove empty cells (about 800) bu first replacing empty with NA, then subsetting to non-NA containing rows 
simplified.df$name[simplified.df$name==""] <- NA
simplified.df <- na.omit(simplified.df)

fwrite(x = simplified.df, file = paste0(middle_dir,"guppyfat_simplified"))
fwrite(x = complete.df, file = paste0(middle_dir,"guppyfat_complete"))






