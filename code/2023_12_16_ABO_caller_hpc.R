# <!-- This is a description for my R document. -->

## Version control
## Repository: https://github.com/Karnes-Lab/ABO_Allele_Caller
## Commit:

# Print session information
sessionInfo()

# Install and load packages
rm(list = ls(all = TRUE))

# Install and load packages
packages <- c("tidyverse","ggplot2","doParallel", "Biostrings","dplyr","XML")

# ------------ If packages are not already installed run next chuck ---------------

# you may need to change ( or remove) the 'lib' argument 
# Install packages not yet installed


# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages], lib =" ")
# }

#----------------------------------------------------------------------


# Packages loading
lapply(packages, library, character.only = TRUE)

# Load in ABO allele names
abo_alleles_names <- read_csv("ABO_alleles_final_csv_allele_names.csv", col_names = TRUE)

scoring_matrix_1 <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
scoring_matrix_2 <- nucleotideSubstitutionMatrix(match = 1, mismatch = -4, baseOnly = TRUE)

# Set up parallel processing
cores <- detectCores()
cl <- makeCluster(92) 
registerDoParallel(cl)
start_time <- Sys.time()

# Reading in the fasta files as command line arguments
fasta_file1 <- "/groups/kianalee/jgiles/ABO_alleles_fixed_asterisk_rmperiods.fasta"
reference <- readDNAStringSet(fasta_file1, format = "fasta")
reference

output_directory <- "/groups/kianalee/jgiles/results/2023_12_16/"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

file_name_prefix <- "2023_12_16_ABO_alleles_onehundred_analyzed_"


for (k in 1:10) {
  fasta_file2 <- paste0("/groups/kianalee/jgiles/simulated_genotypes/onehundred_bases_mutated/ABO_mutated_100_variants_per_alleles_",k,".fasta")
  fasta_file_simulated <- readDNAStringSet(fasta_file2, format = "fasta")
  onehundred_mut_align_df <- data.frame(matrix(nrow = 280, ncol = 280))
  colnames(onehundred_mut_align_df) <- c(abo_alleles_names$ABO_Allele)
  rownames(onehundred_mut_align_df) <- paste0("patients_", abo_alleles_names$ABO_Allele)
  
  for (jow in 1:280){
    onehundred_mut_align_df[jow,] <- foreach (i = 1:280, .packages = "Biostrings", .combine = "cbind") %dopar% {
      alighn <- pairwiseAlignment(reference[[jow]], fasta_file_simulated[[i]],  substitutionMatrix=scoring_matrix_2)
      #print(paste0("Simulated faste file #", k, "has been ran in patient", jow, "for column:", i))
      scoree <- score(alighn)
      return(scoree)
    }
  }
  #Finding highest score per row
  onehundred_mut_align_df$highest_column <- as.character(apply(onehundred_mut_align_df, 1, function(row) {
    colnames(onehundred_mut_align_df)[which.max(row)]
  }))
  onehundred_mut_align_df$highest_column
  
  #removing patient string to compare against allele
  onehundred_mut_align_df$Patient_id <- sub("patients_", "", rownames(onehundred_mut_align_df))
  #comparing allele name and patient name
  onehundred_mut_align_df$Incorrect_call <- if_else(onehundred_mut_align_df$Patient_id == onehundred_mut_align_df$highest_column, 0, 1)
  
  index_formatted <- sprintf("%03d", k)
  file_name <- paste0(file_name_prefix, index_formatted, ".csv")
  onehundred_mut_align_df %>% write_csv(file = paste0(output_directory, file_name))
}
print(paste("total time for all 100 fasta files:", format(Sys.time() - start_time)))
# Stop parallel processing
stopCluster(cl)

getwd()
#input_csv_file <- "Q:/PharmPractice/Karnes Workgroup/08-Personal-Folders/jasongiles/ABO_caller_jg/results/2023_12_23_hpc/2023_12_06_ABO_alleles_onehundred_analyzed_"


aggregated_results_dir <- output_directory
results_file_prefix <- file_name_prefix
input_csv_file <- paste0(aggregated_results_dir,results_file_prefix)


full_result_df <- data.frame(matrix(nrow = 280))
for (i in 1:10) {
  file_index <- sprintf("%03d", i)
  analyzed_df <- read_csv(paste0(input_csv_file, file_index,".csv"))
  # Select  columns
  selected_columns <- analyzed_df %>% dplyr::select(highest_column, Patient_id, Incorrect_call)
  colnames(selected_columns) <- paste(c("highest_column", "Patient_id", "Incorrect_call"), "_", file_index, sep="")
  
  # Add the selected columns to the result_dataframe
  full_result_df <- bind_cols(full_result_df, selected_columns)
}

full_result_df %>% select(-1) %>%  write_csv("2023_12_16_ABO_alleles_onehundred_analyzed_complete.csv")
