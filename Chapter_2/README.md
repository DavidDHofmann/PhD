# Dispersal and Connectivity in Increasingly Extreme Climatic Conditions

[https://doi.org/10.5061/dryad.z34tmpgnm](https://doi.org/10.5061/dryad.z34tmpgnm)

This repository contains all `R`-code and data to reproduce the analyses and visualizations from Hofmann et al., 2024. It is recommended to explore the data through the provided `R`-Scripts. A general design principle was to compartmentalize all analyses and simulations to reduce computational requirements. As such, there are often parent files (in `.rds` format) that provide overviews and bundle further data-files as `tidyverse` data-tibbles. All `R`-codes are extensively documented, giving detailed insights into the processing and analytical steps.

## Description of the File Structure

The file structure of the repository is as follows:

```bash
├── 02_R-Scripts
│   ├── 00_Functions.cpp                            # Custom C++ functions
│   ├── 00_Functions.R                              # Custom R functions
│   ├── 01_Analysis                                 # R codes for the simulation analysis
│   │   ├── 01_LandCover.R                          # Preparation of land cover data
│   │   ├── 02_SourceAreas.R                        # Generation of source areas from which to simulate dispersal
│   │   ├── 03_Simulation.R                         # Dispersal simulation under a minimum and maximum flood
│   │   ├── 04_HeatmapsBetweenness.R                # Generation of heatmaps and betweenness maps from simulations
│   │   ├── 05_Interpatch.R                         # Computing interpatch connectivity from simulations
│   │   ├── 06_ExtentFlood.R                        # Computing flood extent
│   │   └── 07_Distance.R                           # Identify simulated individuals in the vicinity to humans
│   ├──  02_Visualization                           # R codes for visualizations
│   │   ├── Distance.R                              # Figure S9
│   │   ├── Egression.R                             # Figure S5
│   │   ├── Floodmaps.R                             # Figure 3
│   │   ├── GraphicalAbstract.R                     # Figure 1
│   │   ├── InterpatchTable.R                       # Figure 5, Figure S2
│   │   ├── Metrics.R                               # Figures 6, S3, S4, S6, S7, S8, S10
│   │   ├── Model.R                                 # Figure S1
│   │   └── StudyArea.R                             # Figure 2
│   └── 99_SessionInformation                       # Folder containing R session information
│       └── ...
├── 03_Data
│   ├── 01_RawData                                  # Unprocessed raw data
│   │   ├── FLOODMAPS                               # Remote sensed floodmaps in .tif raster format. Filename indicates remote sensing date (YYYY.MM.DD)
│   │   ├── GLOBELAND                               # Static water map (used to represent areas that are not affected by flooding) in .tif raster format
│   │   ├── MERIT                                   # MERIT river data in .tif raster format
│   │   └── MODIS                                   # Modis Continuous Vegetation Fields data in .tif raster format
│   ├── 02_CleanData                                # Cleaned and pre-processed data
│   │   ├── Africa.shp                              # Shapefile of African continent used for visualizations
│   │   ├── AreasOfInterest.shp                     # Shapefile of areas of interest within which we compared human wildlife conflict
│   │   ├── Cutlines.shp                            # Shapefile of cutlines used to delineate egression zones
│   │   ├── DistanceToHumans.tif                    # Distance to human influence (villages, roads, agriculture) as .tif raster format
│   │   ├── DistanceToWater.tif                     # Distance to water (major waters and rivers) in .tif raster format
│   │   ├── Faults.shp                              # Shapefile of major fault lines near the Okavango Delta for visualizations
│   │   ├── GammaDistribution.rds                   # Gamma distribution used to sample step-lengths for the step-selection simulation as an R-file (.rds)
│   │   ├── HumanInfluence.tif                      # Human influence estimates in .tif raster format
│   │   ├── KAZA.shp                                # Shapefile of KAZA-TFCA borders for visualizations
│   │   ├── MababeDepression.shp                    # Shapefile of MababeDepression lines for visualizations
│   │   ├── MajorRivers.shp                         # Shapefile of major river lines for visualizations
│   │   ├── MajorWaters.shp                         # Shapefile of major water bodies for visualizations
│   │   ├── MovementModel.rds                       # Movement model using which dispersal is simulated as an R-file (.rds)
│   │   ├── Protected.shp                           # Shapefile of protected areas for visualizations
│   │   ├── ReferenceRaster.tif                     # Reference .tif raster (used to crop other rasters)
│   │   ├── ReferenceShape.shp                      # Reference .shp shapefile (used to crop other shapes)
│   │   ├── Roads.shp                               # Shapefile of roads for visualizations
│   │   ├── Scaling.rds                             # Values used to scale extracted covariates
│   │   ├── ShrubCover.tif                          # Shrub cover in .tif raster format
│   │   ├── SourceAreas.shp                         # Shapefile of source areas from which dispersal was simulated
│   │   ├── TreeCover.tif                           # Tree cover in .tif raster format
│   │   ├── Villages.shp                            # Shapefile of villages for visualizations
│   │   └── WaterCover.tif                          # Water cover during the extreme scenarios in .tif raster format
│   └── 03_Results                                  # Model and simulation results
│       ├── 99_Betweenness                          # Folder containing betweenness estimates for each source area in .tif raster format
│       ├── 99_Distances                            # Folder containing estimates of human wildlife conflict (i.e. low distance to humans) for each source area in .tif raster format
│       ├── 99_Heatmaps                             # Folder containing heatmaps for each source area in .tif raster format
│       ├── 99_Simulations                          # Folder containing simulated dispersal trajectories as R-files (.rds)
│       ├── BetweennessGlobal.tif                   # Consolidated betweenness estimates (from 99_Betweenness) across different source areas in .tif raster format
│       ├── BetweennessLocal.tif                    # Consolidated betweenness estimates (from 99_Betweenness) per source area in .tif raster
│       ├── DispersalSimulations.rds                # Consolidated dispersal simulations (from 99_Simulations) as .rds R-files.
│       ├── Distance.rds                            # Consolidated estimates of human wildlife conflict (i.e. low distance to humans, from 99_Distances) as .rds R-file
│       ├── Distance.tif                            # Consolidated estimates of human wildlife conflict (i.e. low distance to humans, from 99_Distances) in .tif raster format
│       ├── DistanceAOI.rds                         # Human wildlife conflict within areas of interest as .rds R-file
│       ├── HeatmapsBetweennessGlobal.rds           # Consolidated heatmaps and betweenness maps across source areas as .rds R-files
│       ├── HeatmapsBetweennessLocal.rds            # Consolidated heatmaps and betweenness maps per source area as .rds R-files
│       ├── HeatmapsGlobal.tif                      # Consolidated heatmaps (from 99_Heatmaps) across different source areas in .tif raster format
│       ├── HeatmapsLocal.tif                       # Consolidated heatmaps (from 99_Heatmaps) per source area in .tif raster
│       └── InterpatchConnectivityBootstrapped.rds  # Estimates of interpatch connectivity as .rds R-files
└── README.md                                       # This readme
```

## Description of Data Files
### 01_RawData
Raw data was further processed

- `FLOODMAPS/YYYY.MM.DD.tif`
  - 0: Water
  - 127: Cloud cover
  - 255: Dryland

- `GLOBELAND/Water.tif`
  - 0: Dryland
  - 1: Water

- `MERIT/Rivers.tif`
  - 0: Dryland
  - 1: River (only rivers that are > 10 meters wide)

- `MODIS/Shrubs.tif`
  - 0-100: Percentage cover by shrubs
  - 200: Invalid / water

- `MODIS/Trees.tif`
  - 0-100: Percentage cover by trees
  - 200: Invalid / water

### 02_CleanData
Clean data was not further processed. Some of this data was obtained already cleaned from other data-sources (see below).

#### Shapefiles
Shapefiles are in EPSG:4326 projection and can be loaded using `terra::vect()`.

- `Africa.shp:`
  - ID: Running number to identify individual polygons
  - CODE: Country code
  - COUNTRY: Name of the country

- `AreasOfInterest.shp:`
  - ID: Running number to identify individual polygons
  - Name: Name of the area covered by each polygon

- `Cutlines.shp` (generated from the `SourceAreas.R` script):
  - FID: Running number to identify individual cutlines

- `Faults.shp`:
  - id: Running number to identify individual fault lines
  - Name: Name of the fault line

- `KAZA.shp`:
  - Name: Name of the area

- `MababeDepression.shp`:
  - Name: Name of the area

- `MajorRivers.shp`:
  - Name: Name of each river

- `MajorWaters.shp`:
  - Name: Name of each water source
  - Category: Category of the water source

- `Protected.shp`:
  - Name: Name of the protected area
  - IUCN: IUCN category of the protected area
  - Country: Country in which the protected area lies
  - Desig: Designation of the protected area (reclassified from the world database of protected areas into national park, protected area, forest reserve)
  - Values: Numerical representation of the designation (3 = national park, 2 = protected area, 1 = forest reserve)

- `SourceAreas.shp` (generated from the `SourceAreas.R` script):
  - ID: Running number to identify the source areas
  - Type: Category of the source area
    - Main: Areas within the main study area
    - Buffer: Egression zone

- `Villages.shp`:
  - name: name of the village
  - place: category of the village (Village or City)

#### Raster data
Raster data is in EPSG:4326 projection and can be loaded using `terra::rast()`.

- `DistanceToHumans.tif`:
  - Values: Distance (in meters) to the nearest human influenced grid cell (roads, settlement, agriculture)

- `DistanceToWater.tif`:
  - Values: Distance (in meters) to the nearest grid cell covered by water for a minimum, average, and maximum flood scenario.

- `HumanInfluence.tif`:
  - Values: Relative strength of human influence (roads, settlements, and agriculture). The derivation of this layer is described by Hofmann et al., 2021.

- `ShrubCover.tif`:
  - Values: Proportion of shrub cover (0 to 1) for a minimum, average, and maximum flood scenario.

- `TreeCover.tif`:
  - Values: Proportion of tree cover (0 to 1) for a minimum, average, and maximum flood scenario.

- `WaterCover.tif`:
  - Values: Proportion of water cover (0 or 1) for a minimum, average, and maximum flood scenario.

#### R-Data
R-Data files are in `.rds` format and can be loaded into R using `readr::read_rds()`.

- `MovementModel.rds` (i.e. the dispersal model):
  - Model of the class `glmmTMB`, which requires the `glmmTMB` `R`-package to be installed and loaded. This model was used to predict selection scores. Covariates extracted during the simulation need to be scaled using the scaling parameters stored in `Scaling.rds`.

- `GammaDistribution.rds` (used to propose random steps):
  - Shape: Estimated shape parameter of the step-length distribution
  - Scale: Estimated scale parameter of the step-length distribution

- `Scaling.rds` (used to scale extracted covariates):
  - center: value by which covariates are shifted
  - scale: value by which covariates are scaled

### 03_Results
The following files will all be generated when running through the analyses. Intermediate results are stored in the `99_...` folders, and only later consolidated into single files. The `.rds` files give overviews over the exact simulation parameters and associated files. For instance, the file `DispersalSimulations.rds` consolidates the individual simulations stored under `99_Simulations`.

Overviews / Consolidated files:

- `DispersalSimulation.rds`:
  - x: Simulated longitude
  - y: Simulated latitude
  - absta_: absolute turning angle (in radians)
  - ta_: relative turning angle (in radians)
  - sl_: step length (in meters)
  - Timestamp: POSIXct timestamp
  - BoundaryHit: Logical indicator if the simulated individual has hit a map boundary or not
  - inactive: Logical indicator during what time a step is taken
    - 0 = Inactive
    - 1 = Active
  - TrackID: Identifier of the simulated track
  - SimID: Identifier of the simulation number
  - FloodLevel: Flood scenario under which dispersal was simulated (either minimum or maximum flood, see `Water.tif`)
  - SourceArea: Source area from which individuals were simulated (i.e., the id of the `SourceArea.shp`)
  - StepNumber: Step number of the simulated step (1-2000)

- `Distance.rds`:
  - FloodLevel: Flood scenario (either minimum or maximum flood, see `Water.tif`)
  - SourceArea: Source area from which individuals were simulated (i.e., the id of the `SourceArea.shp`)
  - Filename: Filepath to where the raster files are stored.
  - Level: Whether the metric is computed for a specific source area (Local) or if it was compute across all source areas (Global).

- `DistanceAOI.rds`:
  - FloodLevel: Flood scenario (either minimum or maximum flood, see `Water.tif`)
  - Name: Name of the area of interest, see `AreasOfInterest.shp`
  - Number: Density of trajectories within the area of interest
  - SE: Standard error of the density of trajectories within the area of interest
  - Percent: Percent change of the density from the minimum to maximum flood
  - combined: Latex code for the number and standard error. This can be included in reports

- `HeatmapsBetweennessGlobal.rds`:
  - Steps: Number of steps after which the metric was computed (500, 1000, or 2000)
  - FloodLevel: Flood scenario (either minimum or maximum flood, see `Water.tif`)
  - FilenameHeatmap: Filepath to where the raster file of the resulting heatmap is stored
  - FilenameBetweenness: Filepath to where the raster file of the resulting betweenness map is stored
  - Heatmap: Heatmap as raster file
  - Betwenness: Betweenness map as raster file

- `HeatmapsBetweennessLocal.rds`:
  - Steps: Number of steps after which the metric was computed (500, 1000, or 2000)
  - SourceArea: Source area from which individuals were simulated (i.e., the id of the `SourceArea.shp`)
  - FloodLevel: Flood scenario (either minimum or maximum flood, see `Water.tif`)
  - FilenameHeatmap: Filepath to where the raster file of the resulting heatmap is stored
  - FilenameBetweenness: Filepath to where the raster file of the resulting betweenness map is stored
  - Heatmap: Heatmap as raster file
  - Betwenness: Betweenness map as raster file

- `InterpatchConnectivityBootstrapped.rds`:
  - SourceArea: Origin of the inter-patch connection (see `SourceAreas.shp`)
  - CurrentArea: Target of the inter-patch connection (see `SourceAreas.shp`)
  - FloodLevel: Considered flood level (min or max)
  - Type: Type of the interpatch connection
    - Dispersal = Movement between two main source areas
    - Egression = Movement from a main source area to an egression zone
  - DispersalSuccess: Number of simulated individuals moving between the specified areas
  - SDDispersalSuccess: SD of simulated individuals moving between the specified areas
  - DispersalDuration: Average minimum duration (in steps) it takes before simulated individuals successfully move from the source to the current area
  - DispersalDuration: SD of the minimum duration (in steps) it takes before simulated individuals successfully move from the source to the current area

Unconsolidated files (these files should not be directly accessed):

- `99_Betweenness/Betweenness_Steps_FloodLevel.tif`:
  - Values: Betweenness values calculated from simulated dispersers across all source areas

- `99_Betweenness/Betweenness_Steps_SourceArea_FloodLevel.tif`:
  - Values: Betweenness values calculated from simulated dispersers separately for each source area

- `99_Betweenness/Distance_SourceArea_FloodLevel.tif`:
  - Values: Number of simulated coordinates within 500 meters of a human influenced grid cell computed separately for each source area

- `99_Heatmaps/Heatmap_Steps_FloodLevel.tif`:
  - Values: Traversal frequency calculated from simulated dispersers across all source areas

- `99_Heatmaps/Heatmap_Steps_SoufceArea_FloodLevel.tif`
  - Values: Traversal frequency calculated from simulated dispersers separately for each source area

- `99_Simulations/FloodLevel_SourceArea_Replicate.rds`
  - x: Simulated longitude
  - y: Simulated latitude
  - absta_: absolute turning angle (in radians)
  - ta_: relative turning angle (in radians)
  - sl_: step length (in meters)
  - Timestamp: POSIXct timestamp
  - BoundaryHit: Logical indicator if the simulated individual has hit a map boundary or not
  - inactive: Logical indicator if the step is during wild dogs' active (1) or inactive (0) phase
  - TrackID: Identifier of the simulated track

- `99_BetweennessGlobal.tif`
  - Values: Betweenness across all source areas

- `99_BetweennessLocal.tif`
  - Values: Betweenness for each source area separately

- `99_HeatmapsGlobal.tif`
  - Values: Betweenness across all source areas

- `99_HeatmapsLocal.tif`
  - Values: Heatmaps for each source area separately

## Sharing/Access information

Links to other publicly accessible locations of the data:

- Description of the spatial data preparation: [Hofmann et al., 2021](https://doi.org/10.1111/1365-2664.13868)
- Description of the dispersal model: [Hofmann et al., 2023](https://doi.org/10.1007/s10980-023-01602-4)
