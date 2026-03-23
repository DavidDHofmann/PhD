################################################################################
#### Copy Directories for the Online Repository
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)

# Define output directory
inpdir <- "D:/SwitchDrive/02_Academia/02_PhD/Chapter_3"
outdir <- "D:/SeasonalDispersal"
setwd(inpdir)

################################################################################
#### Dryad and Github
################################################################################
# Identify all folders to copy
tocopy <- tibble(Dir = list.dirs(recursive = F, full.names = F))
tocopy <- subset(tocopy, grepl(pattern = "Scripts|Data|Manuscript", Dir))

# Identify files within those folders. Note that we don't want to copy all
# files from the Manuscript folder.
tocopy$Files <- lapply(tocopy$Dir, function(x) {
  if (x == "04_Manuscript") {
    files <- dir(path = x, full.names = F, recursive = T)
    files <- files[grepl(files, pattern = "png|jpg|odp|pptx|svg|Figures/.*tex|Figures/.*pdf")]
  } else if (x == "03_Data") {
    files <- dir(path = x, full.names = F, recursive = T)
    files <- files[!grepl(files, pattern = "01_RawData")]
  } else {
    files <- dir(path = x, full.names = F, recursive = T)
  }
  return(files)
})

# Unnest the files and create copy-path
tocopy         <- unnest(tocopy, Files)
# tocopy$Where   <- ifelse(grepl(tocopy$Dir, pattern = "02_R-Scripts"), "GitHub", "Dryad")
tocopy$OldPath <- file.path(tocopy$Dir, tocopy$Files)
tocopy$NewPath <- file.path(outdir, tocopy$OldPath)

# Remove scripts for presentation figures
remove <- grepl(tocopy$NewPath, pattern = "03_Visualization/Presentation")
tocopy <- tocopy[!remove, ]

# Adjust paths of remaining visualization scripts
tocopy$NewPath <- gsub(tocopy$NewPath
  , pattern     = "02_Visualization/Manuscript"
  , replacement = "02_Visualization"
)

# Further adjustments
# tocopy$NewPath <- gsub(tocopy$NewPath
#   , pattern     = "02_R-Scripts/"
#   , replacement = ""
# )
# tocopy$NewPath <- gsub(tocopy$NewPath
#   , pattern     = "03_Data/"
#   , replacement = ""
# )
tocopy$NewPath <- gsub(tocopy$NewPath
  , pattern     = "04_Manuscript/Figures"
  , replacement = "04_Manuscript"
)

# Drop manuscript figures
tocopy <- subset(tocopy, !grepl(NewPath, pattern = "04_Manuscript"))

# Create directories that do not exist yet
dirs <- unique(dirname(tocopy$NewPath))
dirs <- dirs[!dir.exists(dirs)]
if (length(dirs) > 0) {
  for (i in dirs) {
    dir.create(i, recursive = T, showWarnings = F)
  }
}

# Which files don't exist yet?
tocopy$NewPath[!file.exists(tocopy$NewPath)]

# Copy files
file.copy(
    from = tocopy$OldPath
  , to   = tocopy$NewPath
)

# # Now we need to go through all R-files and remove the working directory
# scripts <- dir(path = file.path(outdir, "02_R-Scripts"), pattern = ".R$", recursive = T, full.names = T)
# for (i in scripts) {
#   cont     <- readLines(i)
#   toremove <- which(grepl(x = cont, pattern = "setwd|wd <-|Change the working"))
#   if (length(toremove) > 0) {
#     cont <- cont[-toremove]
#   }
#   cont <- gsub(x = cont, pattern = "04_Manuscript/Figures/", replacement = "04_Manuscript/")
#   writeLines(cont, i)
# }

# Copy readme
file.copy("README.md", file.path(outdir, "README.md"))

#  Remove files that are not needed
file.remove(file.path(outdir, "03_Data/02_CleanData/Dispersers.csv"))
file.remove(file.path(outdir, "03_Data/02_CleanData/DispersersSubsampled.csv"))
file.remove(file.path(outdir, "03_Data/02_CleanData/ValidationDispersers.rds"))
file.remove(file.path(outdir, "03_Data/02_CleanData/SSF.csv"))
file.remove(file.path(outdir, "03_Data/02_CleanData/SSFExtracted.rds"))
