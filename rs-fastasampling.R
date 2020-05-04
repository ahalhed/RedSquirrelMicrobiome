#---------------------
# Fasta File Sampling
# Alicia Halhed
# December 2019
#---------------------
#current path to this file: ./OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/fastasampling.R 
# set working directory to where the data is saved
setwd("~/Desktop/Guelph/red-squirrel-data") #current path to files
# attach required packages (most as just loaded as needed in the sampling script)
library(tidyverse)
library(stringr)
# read in the metadata
rs_metadata <- readxl::read_xlsx("/Users/aliciahalhed/Desktop/Guelph/red-squirrel-data/original/metadata_all_samples.xlsx")
# so they actually typed in NAs for this particular file, replace those with actual NAs
rs_metadata <- rs_metadata %>% 
  mutate_all(funs(str_replace(., "NA", "replace")))
rs_metadata[ rs_metadata == "replace" ] <- NA
# take a sample of the metadata (in this case, 5% of the total sample size)
rs_metadata %>% 
  group_by(Grid) %>% 
  sample_frac(.05) -> rs_metadata_sample
# read the full sequence information to a data frame in R
rs_fasta <- Biostrings::readDNAStringSet("/Users/aliciahalhed/Desktop/Guelph/red-squirrel-data/rs-QIIME2-AH/RS_seqs.fasta")
seq_name <- names(rs_fasta) %>% str_replace("_", " ")
sequence <- paste(rs_fasta)
rs_fasta_df <- data.frame(seq_name, sequence)
# create a new column with only the sampleID information (without run information)
rs_fasta_df$SampleID <- word(rs_fasta_df$seq_name, 1, sep = " ")
# join the metadata to the sequence information by sampleID
rs_subset_df <- inner_join(rs_fasta_df, rs_metadata_sample, by = 'SampleID')
# pull out the sequence name and nucleotides to write to fasta
#rs_subset_df %>% select(seq_name, sequence) %>% seqRFLP::dataframe2fas("RS_sample.fasta")
#rename(seq_name = seq.name, seq.text = sequence) %>%
#phylotools::dat2fasta(outfile = "RS_sample.fasta")