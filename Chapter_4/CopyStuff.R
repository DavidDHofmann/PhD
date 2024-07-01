################################################################################
#### Copy Directories for the Online Repository
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)

# Define output directory
inpdir <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
outdir <- "/home/david/Schreibtisch"
outnam <- "IrregularSSF"
setwd(inpdir)

# Identify all folders to copy
tocopy <- tibble(Dir = list.dirs(recursive = F, full.names = F))
tocopy <- subset(tocopy, grepl(pattern = "Scripts|Data|Manuscript", Dir))

# Identify files within those folders. Note that we don't want to copy all
# files from the Manuscript folder.
tocopy$Files <- lapply(tocopy$Dir, function(x) {
  if (x == "04_Manuscript") {
    files <- dir(path = x, full.names = F, recursive = T)
    files <- files[grepl(files, pattern = "png|jpg|odp")]
  } else {
    files <- dir(path = x, full.names = F, recursive = T)
  }
  return(files)
})

# Unnest the files and create copy-path
tocopy         <- unnest(tocopy, Files)
tocopy$OldPath <- file.path(tocopy$Dir, tocopy$Files)
tocopy$NewPath <- file.path(outdir, outnam, tocopy$OldPath)
# head(tocopy$OldPath)
# head(tocopy$NewPath)

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

# Also copy the readme, glossary, and Rproj file
file.copy("README.md", file.path(outdir, outnam, "README.md"))
file.copy("Glossary.csv", file.path(outdir, outnam, "Glossary.csv"))
file.copy("IrregularSSF.Rproj", file.path(outdir, outnam, "IrregularSSF.Rproj"))

# Now we need to go through all R-files and remove the working directory
scripts <- dir(path = file.path(outdir, outnam, "02_R-Scripts"), recursive = T, full.names = T)
for (i in scripts) {
  cont     <- readLines(i)
  toremove <- which(grepl(x = cont, pattern = "setwd|wd <-|Set working"))
  if (length(toremove) > 0) {
    cont <- cont[-toremove]
  }
  writeLines(cont, i)
}

