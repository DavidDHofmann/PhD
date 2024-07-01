################################################################################
#### Copy Directories for the Online Repository
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)

# Define output directory
inpdir <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
outdir <- "/home/david/Schreibtisch/FloodDispersal"
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
remove <- grepl(tocopy$NewPath, pattern = "02_Visualization/Presentation")
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

# Now we need to go through all R-files and remove the working directory
scripts <- dir(path = file.path(outdir, "02_R-Scripts"), pattern = ".R$", recursive = T, full.names = T)
for (i in scripts) {
  cont     <- readLines(i)
  toremove <- which(grepl(x = cont, pattern = "setwd|wd <-|Change the working"))
  if (length(toremove) > 0) {
    cont <- cont[-toremove]
  }
  cont <- gsub(x = cont, pattern = "04_Manuscript/Figures/", replacement = "04_Manuscript/")
  writeLines(cont, i)
}

# Copy readme
file.copy("README.md", file.path(outdir, "README.md"))

#  Remove files that are not needed
file.remove(file.path(outdir, "02_R-Scripts/01_Analysis/00_Setup.R"))
file.remove(file.path(outdir, "03_Data/02_CleanData/HumanDensity.tif"))
file.remove(file.path(outdir, "03_Data/02_CleanData/CattleDensity.tif"))
file.remove(file.path(outdir, "03_Data/02_CleanData/Agriculture.tif"))
file.remove(file.path(outdir, "03_Data/02_CleanData/Dispersers.csv"))
file.remove(file.path(outdir, "03_Data/02_CleanData/DispersersUpdated.csv"))

################################################################################
#### Global Change
################################################################################
# Define output directory
inpdir <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
outdir <- "/home/david/Schreibtisch/GlobalChange"
setwd(inpdir)

# Files to rename
from <- c(
    "Manuscript.tex"
  , "Manuscript_Support.pdf"
  , "Literature.bib"
  , "Figures/GraphicalAbstract.pdf"
  , "Figures/StudyArea2.png"
  , "Figures/FloodExtent.png"
  , "Figures/InterpatchConnectivityTable2.png"
  , "Figures/InterpatchConnectivityMap.png"
  , "Figures/Metrics.png"
)
to <- c(
    "Manuscript.tex"
  , "Manuscript_Support.pdf"
  , "Literature.bib"
  , "Figure_01.pdf"
  , "Figure_02.png"
  , "Figure_03.png"
  , "Figure_04.png"
  , "Figure_05.png"
  , "Figure_06.png"
)
tocopy <- data.frame(
    OldPath = file.path("04_Manuscript", from)
  , NewPath = file.path(outdir, to)
)

# Create directories that do not exist yet
dirs <- unique(dirname(tocopy$NewPath))
dirs <- dirs[!dir.exists(dirs)]
if (length(dirs) > 0) {
  for (i in dirs) {
    dir.create(i, recursive = T, showWarnings = F)
  }
}

# Copy files
file.copy(
    from = tocopy$OldPath
  , to   = tocopy$NewPath
)

# Adjust figure names
manu <- file.path(outdir, "Manuscript.tex")
cont <- readLines(manu)
for (i in 4:length(from)) {
  cont <- gsub(x = cont, pattern = from[i], replacement = to[i])
}
writeLines(text = cont, manu)
