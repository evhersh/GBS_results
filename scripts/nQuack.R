library(dplyr)
library(nQuack)
library(kableExtra)
library(stringr)

availableCores()
progressr::handlers(global=TRUE)
future::plan("multisession", workers=15)

# PREPARE DATA ----

# Set in and out paths of files
inpath <- "./bam/"
outpath <- "./prepared/"

# List files in the inpath and remove their ending
filelist <- list.files(path = inpath, pattern = "*.bam" )
filelist <- gsub(".bam", "", filelist)

for( i in 1:length(filelist)){
  prepare_data(filelist[i], inpath, outpath)
}

# PROCESS DATA ----

inpathtext <- "./data/prepared/fix/"
newfilelist <- list.files(path = inpathtext, pattern = "*.txt")

for(i in 1:length(newfilelist)){
  samp <- newfilelist[i]
  temp <- process_data(paste0(inpathtext, samp), 
                       min.depth = 10, 
                       max.depth.quantile.prob = 0.9, 
                       error = 0.05, 
                       trunc = c(0.15,0.85))
  
  
  write.csv(temp, 
            file = paste0("./data/processed/fix/", gsub(".txt", "", samp), ".csv"),
            row.names = FALSE)
}

## more filtering ----

### move files with ploidy mismatch to new folder ----
# Read the nQuack results CSV
nquack_results <- read.csv("./data/output/nQuack_v1.csv")

# Filter for rows where match is 'no'
to_fix <- nquack_results %>% 
  filter(tolower(match) == "no") %>% 
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

### process without truncation for inspection ----

inpathtext <- "./data/prepared/fix/"
newfilelist <- list.files(path = inpathtext, pattern = "*.txt")

for(i in 1:length(newfilelist)){
  samp <- newfilelist[i]
  temp <- process_data(paste0(inpathtext, samp), 
                       min.depth = 2, 
                       max.depth.quantile.prob = 0.9, 
                       error = 0.01, 
                       trunc = c(0,0))
  
  
  write.csv(temp, 
            file = paste0("./data/processed/notrunc/", gsub(".txt", "", samp), ".csv"),
            row.names = FALSE)
}

### try bclean ----

inpathtext <- "./data/processed/notrunc/"
newfilelist <- list.files(path = inpathtext, pattern = "*.csv")

# Create output directory if it doesn't exist
dir.create("./data/processed/bclean/", recursive = TRUE, showWarnings = FALSE)

for(i in seq_along(newfilelist)) {
  samp <- as.matrix(read.csv(paste0(inpathtext, newfilelist[i])))
  temp <- Bclean(samp)
  
  # Remove ".csv" from original filename to use in output
  outname <- gsub(".csv", "", newfilelist[i])
  
  write.csv(temp, 
            file = paste0("./data/processed/bclean/", outname, "_bclean.csv"),
            row.names = FALSE)
}


### inspect ----

process_ratio <- function(name_part, plot = TRUE) {
  # Construct the full path using the name_part
  filename <- paste0(name_part, "-RG.csv")
  filepath <- file.path("./data/processed/fix", filename)
  
  # Read the file and convert to matrix
  xm <- as.matrix(read.csv(filepath))
  
  # Check if matrix has at least two columns
  if (ncol(xm) < 2) {
    stop("The input file must have at least two columns.")
  }
  
  # Calculate the ratio
  xm.ab <- xm[, 2] / xm[, 1]
  
  # Plot histogram if requested
  if (plot) {
    hist(xm.ab, main = paste("Histogram of", filename), xlab = "Ratio (Col2 / Col1)")
  }
  
  return(xm.ab)
}

process_ratio("C23-A_4")

samples <- list.files(path = "./data/processed/", pattern = "*.csv")
samples <- gsub(".csv", "", samples)  # Remove .csv to get sample names

xm <- as.matrix(read.csv("./data/processed/fix/C85-A_3-RG.csv"))
xm.ab <- xm[,2]/xm[,1]
hist(xm.ab)
xm.dn <- denoise_data(xm)
xm.bc <- Bclean(xm)

# MODEL INFERENCE ----
samples <- list.files(path = "./data/processed/fix", pattern = "*.csv")
samples <- gsub(".csv", "", samples)  # Remove .csv to get sample names

# skip this one?
for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("./data/processed/fix/", samples[i], ".csv")))
  out1 <- quackNormal(xm = temp, samplename = samples[i], cores = 15, parallel = TRUE)
  out2 <- quackBeta(xm = temp, samplename = samples[i], cores = 15, parallel = TRUE)
  out3 <- quackBetaBinom(xm = temp, samplename = samples[i], cores = 15, parallel = TRUE)
  allout <- rbind(out1, out2, out3)
  write.csv(allout, 
            file = paste0("./data/output/fix/", samples[i], ".csv"),
            row.names = FALSE)
}

## examine best model ----
inpathtext <- "./data/output/fix/"
samples <- c("C23-A_1-RG", "C23-A_2-RG", "C23-A_3-RG")

for(i in 1:length(samples)){
  temp <- read.csv(paste0(inpathtext, samples[i], ".csv"))
  summary <- quackit(model_out =  temp, 
                     summary_statistic = "BIC", 
                     mixtures = c("diploid", "triploid", "tetraploid"))
  write.csv(summary, 
            file = paste0("./data/output/fix/model/", samples[i], ".csv"),
            row.names = FALSE)
}

key <- data.frame(sample = c("C23-A_1-RG", "C23-A_2-RG", "C23-A_3-RG"), 
                  ploidal.level = c("diploid","triploid", "tetraploid"))

# Read in quackit() output
dfs <- lapply(list.files("./data/output/fix/model/", full.names = TRUE  ), read.csv)
alloutput <- do.call(rbind, dfs)

# Combined
alloutputcombo <- dplyr::left_join(alloutput, key)

# Check the accuracy
alloutputcombo <- alloutputcombo %>%
  dplyr::mutate(accuracy = ifelse(winnerBIC == ploidal.level, 1, 0))

## What distribution and model type should we use?
sumcheck <- alloutputcombo %>% 
  group_by(Distribution, Type) %>% 
  summarize(total = n(), correct = sum(accuracy))

kbl(sumcheck) %>%
  kable_paper("hover", full_width = F) 

# run only the best models?
out <- c()

for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("./data/processed/fix/", samples[i], ".csv")))
  out[[i]] <- bestquack(temp, 
                        distribution = "normal",
                        type = "fixed_2",
                        uniform = 1, 
                        mixtures = c("diploid", "triploid", "tetraploid"), 
                        samplename = samples[i])
}

saveRDS(out, "./data/output/fix/bestquack.rds")

## bootstrap ----
bout <- c()

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

# join with metadata
sample_names <- indNames(AllPops.gc) # From adegenet/poppr
ploidy_vals <- AllPops.gc@ploidy

meta.df <- data.frame(sample_name = sample_names,
                      ploidy = ploidy_vals,
                      stringsAsFactors = FALSE)

# 3. Extract short sample names from bout.df (supporting -S or -A)
bout.df <- bout.df %>%
  mutate(short_sample = str_extract(sample, "^\\w+-[SA]_\\d+"))

# 4. Do the same for meta.df
meta.df <- meta.df %>%
  mutate(short_sample = str_extract(sample_name, "^\\w+-[SA]_\\d+"))

# 5. Join
bout.joined <- left_join(bout.df, meta.df, by = "short_sample") %>%
  select(sample_name, diploid, triploid, tetraploid, ploidy) %>%
  rename(assumed_ploidy = ploidy)

write.csv(bout.joined, "./data/output/nQuack_v1.csv")
