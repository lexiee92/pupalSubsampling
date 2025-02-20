### Subsampling Pupal Samples for PU3 and PU5
# By Lexie Edwards (Feb 2025)

library(dplyr)
library(writexl)

################# USER-DEFINED PARAMETERS #########################
sampleBlock <-  "PU3"  # Change this
motherComm <-   "R"    # Change this
commNo <-       "1"    # Change this

num_samples <- c(T1 = 22*8+1, T2 = 15*8+4, T3 = 16*8+1)  
##################################################################

# Define plate dimensions
rows <- LETTERS[1:8]  # Rows A to H
cols <- 1:12          # Columns 1 to 12
wells_per_plate <- 96 

total_samples <- sum(num_samples)
num_plates <- ceiling(total_samples / wells_per_plate)  

# Generate storage plate grid
storagePlates <- expand.grid(storageRow = rows, storageColumn = cols, storagePlate = 1:num_plates) %>%
  arrange(storagePlate, storageColumn, storageRow) %>%
  mutate(storageLocation = paste0(storagePlate, "-", storageRow, sprintf("%02d", storageColumn)),
         timePoint = "Empty")  # Default to "Empty"

# Assign samples across plates
assign_samples <- function(storagePlates, num_samples) {
  current_column <- 1
  current_plate <- 1
  
  for (tp in names(num_samples)) {
    num_tp_samples <- num_samples[tp]
    filled_wells <- 0
    
    while (filled_wells < num_tp_samples) {
      available_wells <- which(storagePlates$storageColumn == current_column & 
                                 storagePlates$storagePlate == current_plate & 
                                 storagePlates$timePoint == "Empty")
      
      if (length(available_wells) == 0) {
        current_column <- ifelse(current_column == 12, 1, current_column + 1)
        if (current_column == 1) current_plate <- current_plate + 1
        if (current_plate > num_plates) stop("Error: Exceeded maximum number of plates available!")
        next
      }
      
      num_to_fill <- min(num_tp_samples - filled_wells, length(available_wells))
      storagePlates$timePoint[available_wells[1:num_to_fill]] <- tp
      filled_wells <- filled_wells + num_to_fill
    }
  }
  return(storagePlates)
}

storagePlates <- assign_samples(storagePlates, num_samples)

# Subsample function
subsample <- function(df, n) {
  if (nrow(df) == 0) return(df)
  return(df[sample(nrow(df), min(nrow(df), n), replace = FALSE), ])
}

# Subsampling per timepoint
subsampled_counts <- c(T1 = 108 * 1.05, T2 = 90 * 1.05, T3 = 87 * 1.05) # extra 5% for backup
T1_sub <- subsample(subset(storagePlates, timePoint == 'T1'), subsampled_counts['T1'])
T2_sub <- subsample(subset(storagePlates, timePoint == 'T2'), subsampled_counts['T2'])
T3_sub <- subsample(subset(storagePlates, timePoint == 'T3'), subsampled_counts['T3'])

# Identify unsampled wells
unsampled <- subset(storagePlates, !storageLocation %in% c(T1_sub$storageLocation, T2_sub$storageLocation, T3_sub$storageLocation))
if (nrow(unsampled) > 0) unsampled$timePoint <- paste0("unsampled-", unsampled$timePoint)

# Combine subsampled and unsampled
final_subsampled <- bind_rows(T1_sub, T2_sub, T3_sub) %>%
  dplyr::select(storagePlate, storageRow, storageColumn, storageLocation, timePoint) %>%
  arrange(timePoint, storagePlate, storageColumn, storageRow)

if (nrow(unsampled) > 0) {
  final_subsampled <- bind_rows(final_subsampled, unsampled)
}

# Import and update extraction plate data
extractionPlates <- read.csv("extractionPlates.csv", stringsAsFactors = FALSE)
extractionPlates <- extractionPlates %>%
  mutate(
    # Convert all columns to character type first, if necessary
    across(everything(), ~ as.character(.), .names = "{.col}"), 
    
    # Now, handle the columns with specific replacements
    samplingBlock = ifelse(is.na(samplingBlock) | samplingBlock == "", sampleBlock, samplingBlock),
    motherCommunity = ifelse(is.na(motherCommunity) | motherCommunity == "", motherComm, motherCommunity),
    communityNumber = ifelse(is.na(communityNumber) | communityNumber == "", commNo, communityNumber),
    flySpeciesID = ifelse(is.na(flySpeciesID), "", flySpeciesID),
    parasitoidPresence = ifelse(is.na(parasitoidPresence), "", parasitoidPresence),
    parasitoidID = ifelse(is.na(parasitoidID), "", parasitoidID),
    comment = ifelse(is.na(comment), "", comment)
  )
# Save to Excel
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") 
output_file <- paste0(sampleBlock, "-", motherComm, commNo, "_", timestamp, ".xlsx")
write_xlsx(list("subsampled" = final_subsampled, "extractionPlates" = extractionPlates), output_file)

# Determine next action
if (nrow(unsampled) > nrow(final_subsampled)) {
  cat("Subsample first")
} else if (nrow(unsampled) < nrow(final_subsampled)) {
  cat("Remove unsampled first")
} else {
  cat("Subsample count equals unsampled count")
}
