library(dplyr)
library(nQuack)
library(kableExtra)
library(stringr)
library(future)

# move files with ploidy mismatch to new folder ----

# Read the nQuack results CSV
nquack_results <- read.csv("./data/output/nQuack_v2.csv")

# Filter for rows where match is 'no'
to_fix <- nquack_results %>% 
  filter(tolower(wrong) == "yes") %>% 
  pull(sample_name)

# Create the destination directory if it doesn't exist
dir.create("./data/prepared/fix/", showWarnings = FALSE, recursive = TRUE)

# List all files in the prepared directory
all_files <- list.files("./data/prepared/", full.names = TRUE)

# Loop through each sample_name and move matching files
for (sample in to_fix) {
  matching_files <- all_files[grepl(sample, basename(all_files))]
  for (file in matching_files) {
    file.rename(file, file.path("./data/prepared/fix/", basename(file)))
  }
}

# Process ----

inpathtext <- "./data/prepared/fix/"
newfilelist <- list.files(path = inpathtext, pattern = "*.txt")

for(i in 1:length(newfilelist)){
  samp <- newfilelist[i]
  temp <- process_data(paste0(inpathtext, samp), 
                       min.depth = 20, 
                       max.depth.quantile.prob = 0.8, 
                       error = 0.05, 
                       trunc = c(0.15,0.85))
  
  
  write.csv(temp, 
            file = paste0("./data/processed/fix/", gsub(".txt", "", samp), ".csv"),
            row.names = FALSE)
}

# Model ----
# bootsrap
bout <- c()

samples <- list.files(path = "./data/processed/fix/", pattern = "*.csv")
samples <- gsub(".csv", "", samples)  # Remove .csv to get sample names

for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("./data/processed/fix/", samples[i], ".csv")))
  bout[[i]] <- quackNboots(temp, 
                           nboots = 100,
                           distribution = "normal",
                           type = "fixed_2",
                           uniform = 1, 
                           mixtures = c("diploid", "triploid", "tetraploid"), 
                           samplename = samples[i])
}

# save all?

saveRDS(bout, "./data/output/fix/bestquack_boots.rds")

bout.df <- bind_rows(bout)

write.csv(bind_rows(bout), "./data/output/fix/bestquack_boots.csv")