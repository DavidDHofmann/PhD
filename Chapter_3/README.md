# Seasonal dynamism or model realism: What drives better predictions of landscape connectivity?

DOI TO PAPER:

This repository contains all `R`-code and data to reproduce the analyses and visualizations from Hofmann et al., 2026. It is recommended to explore the data through the provided `R`-Scripts. A general design principle was to compartmentalize all analyses and simulations to reduce computational requirements. As such, there are often parent files (in `.rds` format) that provide overviews and bundle further data files as `tidyverse` data-tibbles. All `R`-codes are extensively documented, giving detailed insights into the processing and analytical steps. All data can be understood by investigating the scripts for data cleaning. R-scripts are numbered and expected to be run in corresponding sequence.

## Description of the File Structure

The file structure of the repository is as follows:

```bash
├── 02_R-Scripts
│   ├── 00_Functions.R                              # Custom R functions
│   ├── 01_DataCleaning                             # R codes for pre-processing data
│   │   ├── 00_StudyArea.R                          # Definition of study area. Define extent and create shapefile
│   │   ├── 01_GPSData.R                            # Preprocessing of GPS data. Subsetting to dispersing individuals
│   │   ├── 02_LandCover.R                          # Preparation of Globeland land cover data
│   │   ├── 03_HumanInfluence.R                     # Preparation of human influence layer
│   │   ├── 04_Floodmaps.R                          # Preparation of floodmaps for Okavango delta
│   │   ├── 05_Rivers.R                             # Preparation of river network across study area using MERIT
│   │   ├── 06_MergeWater.R                         # Combining static water, dynamic flood, and rivers into one
│   │   ├── 07_Vegetation.R                         # Preparation of continuous vegetation layers
│   │   ├── 08_DynamicWaterAndVegetation.R          # Preparation of temporal water + vegetation layers
│   │   ├── 09_ClimateData.R                        # Downloading and preparing climate data
│   │   ├── 10_PanMappingTraining.R                 # Preparing remote sensing classifier for mapping pans
│   │   ├── 11_SentinelPrerequisits.R               # Identification of files to bulk download from Sentinel
│   │   ├── 12_SentinelDownload.R                   # Bulk download of sentinel satellite imagery
│   │   ├── 13_SentinelProcessing.R                 # Process downloaded data and obtain bottom of atmosphere images
│   │   ├── 14_SentinelPrediction.R                 # Predicting pans from sentinel imagery (applying the classifier)
│   │   ├── 15_SentinelTesting.R                    # Testing the performance of the classifier
│   │   ├── 16_FinalizePanmaps.R                    # Combine classified maps into one, obtain distance to maps
│   │   ├── 17_Moonlight.R                          # Obtain moonlight information
│   │   ├── 18_SeasonalSummary.R                    # Compute "average" seasonal layers
│   │   ├── 19_CombineCovariates.R                  # Overview file to keep track of covariates and their dates
│   │   ├── 20_LookupTable.R                        # Lookup table to know which covariates to use for which date
│   │   └── 21_SourceAreas.R                        # Definition of source areas from which to simulate dispersal
│   ├── 01_Analysis                                 # R codes for the simulation analysis
│   │   ├── 00_SSF.R                                # Preparing random steps for issf
│   │   ├── 01_CovariateExtraction.R                # Extracting covariates along observed and random steps
│   │   ├── 02_MovementModel.R                      # Fitting the iSSF model
│   │   ├── 03_Simulation.R                         # Simulating from the fitted model
│   │   ├── 04_Connectivity.R                       # Obtaining connectivity from simulated trajectories
│   │   ├── 99_ExploratoryAnalysis.R                # Exploratory analysis of movement data
│   │   └── 99_GeneralMetrics.R                     # Summary statistics for paper
│   ├──  02_Visualization                           # R codes for visualizations
│   │   ├── Connectivity.R                          # Figure 6
│   │   ├── Covariates.R                            # Table 1
│   │   ├── DynamicVariables.R                      # Figure S1
│   │   ├── GeneralMetrics.R                        # Figures S2, S3
│   │   ├── Moonlight.R                             # Figure S4
│   │   ├── MovementModel.R                         # Figures 4, S12, Tables S3-S6
│   │   ├── MovementModelInterpretation.R           # Figure S13
│   │   ├── ND.R                                    # Table S2
│   │   ├── PanMapping.R                            # Figure S8, S9, S10, S11
│   │   ├── PansWetDry.R                            # Figure S1
│   │   ├── RankCorrelation.R                       # Figure 5
│   │   ├── Satellites.R                            # Table S1
│   │   ├── SpearmanNumberSteps.R                   # Figure S14
│   │   ├── StudyArea.R                             # Figure 2
│   │   └── TemporalResolution.R                    # Figure 3
│   └── 99_SessionInformation                       # Folder containing R session information in .rds format
│       └── ...
├── 03_Data
│   ├── 02_CleanData                                # Cleaned and pre-processed data
│   │   ├── 00_Floodmaps                            # tif files of floodmaps obtained from MODIS satellite imagery (YYYY-MM-DD)
│   │   ├── 00_NDVI                                 # tif files of NDVI obtained from MODIS through earth engine (YYYY-MM)
│   │   ├── 00_Panmaps                              # tif files of pans /distance to pans predicted from Sentinel II imagery  (YYYY-MM)
│   │   ├── 00_Rainmaps                             # tif files of precipitation from JAXA GSMAP through earth engine (YYYY-MM) 
│   │   ├── 00_Tempmaps                             # tif files of temperature from ERA5 through earth engine (YYYY-MM)
│   │   ├── 00_Vegmaps                              # tif files of continuous vegetation from MODIS (YYYY-MM-DD)
│   │   ├── Africa.gpkg                             # Shapefile of African continent used for visualizations
│   │   ├── Covariates.rds                          # Table of covariates -> For overview of all covariates stored as R-object
│   │   ├── DistanceToPansDynamic.tif               # Stack of layers indicating distance to pans dynamically (i.e. with seasonally)
│   │   ├── DistanceToPansStatic.tif                # Layer indicating distance to pans statically (i.e. without seasonality)
│   │   ├── DistanceToWaterDynamic.tif              # Stack of layers indicating distance to water dynamically
│   │   ├── DistanceToWaterDynamicAggregated.tif    # Stack of layers indicating distance to water dynamically across a typical year
│   │   ├── DistanceToWaterDynamicStatic.tif        # Stack of layers indicating distance to water statically
│   │   ├── Humans.tif                              # Human influence covariate layer
│   │   ├── LandCover.csv                           # Table of land cover categories for the layer below
│   │   ├── LandCover.tif                           # Categorical land cover classes obtained from Globeland30
│   │   ├── LookupTable.rds                         # Lookup table containing sequence of dates and associated seasonal covariate layer
│   │   ├── MajorRivers.gpkg                        # Shapefile of major rivers across the Okavango Delta for visualizations
│   │   ├── MajorWaters.gpkg                        # Shapefile of major water across the Okavango Delta for visualizations
│   │   ├── Moonlight.rds                           # R object containing a table of moonlight statistics
│   │   ├── NDVIDynamic.tif                         # Stack of layers indicating NDVI dynamically
│   │   ├── NDVIDynamicAggregated.tif               # Stack of layers indicating NDVI dynamically across a typical year
│   │   ├── NDVIStatic.tif                          # Stack of layers indicating NDVI statically
│   │   ├── PrecipitationDynamic.tif                # Stack of layers indicating precipitation dynamically
│   │   ├── PrecipitationDynamicAggregated.tif      # Stack of layers indicating precipitation dynamically across a typical year
│   │   ├── PrecipitationStatic.tif                 # Stack of layers indicating precipitation statically
│   │   ├── Protected.gpkg                          # Shapefile of protected areas for visualizations
│   │   ├── Raster.tif                              # Reference raster (used for reprojections and cropping)
│   │   ├── Rivers.tif                              # River raster from the MERIT hydro dataset
│   │   ├── Roads.gpkg                              # Road shapefile obtained from OSM
│   │   ├── Shapefile.gpkg                          # Reference shapefile (used for reprojections and cropping)
│   │   ├── ShrubsDynamic.tif                       # Stack of layers indicating shrubs dynamically
│   │   ├── ShrubsDynamicAggregated.tif             # Stack of layers indicating shrubs dynamically across a typical year
│   │   ├── ShrubsStatic.tif                        # Stack of layers indicating shrubs statically
│   │   ├── Sources.gpkg                            # Shapefile of source areas from which dispersal was simulated
│   │   ├── SSFExtractedRedacted.rds                # Redacted dispersal data (GPS and names omitted). Contains extract covariate data.
│   │   ├── TemperatureDynamic.tif                  # Stack of layers indicating temperature dynamically
│   │   ├── TemperatureDynamicAggregated.tif        # Stack of layers indicating temperature dynamically across a typical year
│   │   ├── TemperatureStatic.tif                   # Stack of layers indicating temperature statically
│   │   ├── TrainingClasses.gpkg                    # Training polygons for the pan mapping algorithm
│   │   ├── TreesDynamic.tif                        # Stack of layers indicating trees dynamically
│   │   ├── TreesDynamicAggregated.tif              # Stack of layers indicating trees dynamically across a typical year
│   │   ├── TreesStatic.tif                         # Stack of layers indicating trees statically
│   │   ├── Villages.gpkg                           # Shapefile of villages for visualizations
│   │   ├── WaterDynamic.tif                        # Stack of layers indicating water dynamically
│   │   ├── WaterDynamicAggregated.tif              # Stack of layers indicating water dynamically across a typical year
│   │   ├── WaterStatic.tif                         # Stack of layers indicating water statically
│   │   ├── WildDogs.gpkg                           # Shapefile of source areas from which dispersal was simulated
│   │   └── Windows.rds                             # Moving windows that were used to download Sentinel II satellite imagery
│   └── 03_Results                                  # Model and simulation results
│       ├── Connectivity                            # Folder containing R-objects with connectivity maps under different configurations
│       ├── Simulations                             # Folder containing R-objects with simulated trajectories under different configurations
│       ├── Validation                              # Folder containing R-objects with validation results under different configurations
│       ├── Connectivity.rds                        # R-object (tibble) providing an overview of connectivity metrics under configurations
│       ├── Formula.rds                             # Formula used for fitting the iSSF model and predicting in the simulations
│       ├── GomotiPred.tif                          # Predicted pan map for validation of pans in the gomoti area
│       ├── MbomaPred.tif                           # Predicted pan map for validation of pans in the Mboma island area
│       ├── MovementModel.rds                       # Movement model results (tibble) under different configurations
│       ├── PanMapping.rds                          # Results from the pan mapping validation
│       ├── RankFrequency.rds                       # Results from the validation of predictive performance under different configurations
│       ├── Scaling.rds                             # Scaling values used to center and standardize covariates
│       ├── SeasonalCovariates.rds                  # Values required to generate figure S1
│       ├── SentinelJoined.rds                      # Overview of all Sentinel II files that must be downloaded
│       ├── SentinelMetadata.rds                    # Metadata of Sentinel II files that are downloaded
│       ├── SentinelResults.rds                     # Sentinel II search results
│       ├── Simulation.rds                          # Simulated dispersal trajectories
│       ├── StepLengthDistribution.rds              # Estimated step length distribution for iSSF.
│       └── Validation.rds                          # Overview of validation results
└── README.md                                       # This readme
```

## Description of Data Files

### 02\_CleanData

Clean data were used for analyses and not further processed .

#### Shapefiles

Shapefiles (`.gpkg`) are in EPSG:4326 projection and can be loaded using `terra::vect()`. Some files may contain superfluous columns that are not further explained here.

- `Africa.gpkg`: Country boundaries obtained using the `rworldmap` R package.
  - ID: Running number identifying individual polygons  
  - CODE: Country code  
  - COUNTRY: Name of the country  

- `MajorRivers.gpkg`: Combination of rivers from OpenStreetMap (https://www.openstreetmap.org) and manually digitized rivers (see Hofmann et al., 2021). Used for visualization only.
  - Name: Name of the river

- `MajorWaters.gpkg`: Large water bodies obtained from OpenStreetMap (https://www.openstreetmap.org). Used for visualization only.
  - name: Name of the waterbody

- `Protected.gpkg`: Protected areas obtained from the World Database on Protected Areas (https://www.protectedplanet.net/en/thematic-areas/wdpa). Used for visualization only.
  - Name: Name of the protected area  
  - IUCN: IUCN category  
  - Country: Country in which the protected area lies  
  - Desig: Reclassified designation (national park, protected area, forest reserve)  
  - Values: Numerical designation  
    - 3: National park  
    - 2: Protected area  
    - 1: Forest reserve  

- `Roads.gpkg`: Road network obtained from OpenStreetMap (https://www.openstreetmap.org). Only major tar roads were used in analyses.
  - fclass: Road type according to OSM classification  

- `Shapefile.gpkg`: Reference shapefile used to crop all spatial data to the study area extent.
  - Name: Name of the study area  

- `Sources.gpkg`: Source areas from which simulated dispersers were released.
  - Name: Name of the source area  

- `TrainingClasses.gpkg`: Training polygons for the pan classification algorithm using Landsat and Sentinel-2 imagery (see paper Appendix). Polygons were manually digitized from Google Earth.
  - Date: Date of the Google Earth satellite image used for creating the training polygon
  - Class: Polygon classification (Dryland, Water, Wetpan)  

- `Villages.gpkg`: Locations of villages within the study area obtained from OpenStreetMap (https://www.openstreetmap.org). Used for visualization only.
  - Name: Name of the village  
  - Place: Type of settlement  

- `WildDogs.gpkg`: IUCN range of African wild dogs obtained from https://www.iucnredlist.org/. Used for visualization only.

#### Raster data
Raster (tif) objects are in EPSG:4326 projection and can be loaded using `terra::rast()`. Many raster covariates are provided in three temporal representations:
- Dynamic: Values correspond to specific timestamps and retain maximum possible temporal variability.
- Static: Values represent a long-term average across all available time steps matching the study period (2011-2022).
- Dynamic Aggregated: Values represent a typical seasonal cycle, obtained by averaging across years while retaining intra-annual (seasonal) variation.

- `WaterDynamic.tif`, `WaterStatic.tif`, `WaterDynamicAggregated.tif`: Binary surface water layers derived from flood mapping (using the `floodmapr` R-package), combined with static water from Globeland30 and river data from MERIT Hydro. Dynamic layers are updated weekly.
  - Values: Binary indicator of absence (0) and presence (1) of water.
    - 0: Dryland
    - 1: Water

- `DistanceToWaterDynamic.tif`, `DistanceToWaterStatic.tif`, `DistanceToWaterDynamicAggregated.tif`: Distance to the nearest non-pan water body. Dynamic layers are updated weekly.
  - Values: Distance (in meters) to the nearest body of water (except pans)

- `DistanceToPansDynamic.tif`, `DistanceToPansStatic.tif`: Distance to the nearest pan containing water. These layers were obtained by applying the trained pan mapping algorithm to seasonally updated Sentinel II imagery. Dynamic layers are updated every 5–10 days. This covariate was not used in the final model due to convergence issues; therefore, no aggregated layer is provided.
  - Values: Distance (in meters) to the nearest pan with water

- `NDVIDynamic.tif`, `NDVIStatic.tif`, `NDVIDynamicAggregated.tif`: Vegetation productivity derived from MODIS MOD13Q1.
  - Values: NDVI between -1 and 1

- `ShrubsDynamic.tif`, `ShrubsStatic.tif`, `ShrubsDynamicAggregated.tif`: Percentage shrub cover derived from MODIS vegetation continuous fields MOD44B. Dynamic layers are updated annually.
  - Values: Percentage shrub cover (between 0 and 1)

- `TreesDynamic.tif`, `TreesStatic.tif`, `TreesDynamicAggregated.tif`: Percentage tree cover derived from MODIS vegetation continuous fields MOD44B. Dynamic layers are updated annually.
  - Values: Percentage tree cover (between 0 and 1)

- `PrecipitationDynamic.tif`, `PrecipitationStatic.tif`, `PrecipitationDynamicAggregated.tif`: Precipitation derived from satellite-based rainfall products.
  - Values: Precipitation in mm

- `TemperatureDynamic.tif`, `TemperatureStatic.tif`, `TemperatureDynamicAggregated.tif`: Air temperature derived from ERA5.
  - Values: Temperature in °C

- `Humans.tif`: Human influence index representing anthropogenic pressure (details on the derivation of this layer are provided by Hofmann et al., 2021).
  - Values: Continuous index of human influence

- `LandCover.tif`: Globeland30 land cover classification, used to extend water mapping beyond dynamically mapped areas.
  - Values:
    - 1: Water
    - 2: Urban
    - 3: Cropland
    - 4: Forest
    - 5: Shrubs
    - 6: Grassland
    - 7: Bare

- `Rivers.tif`: River network derived from the MERIT Hydro dataset. Only rivers wider than 10 meters are retained.
  - Values:
    - 0: No river
    - 1: River (> 10 m)

- `Raster.tif`: Reference raster used for alignment, reprojection, and cropping of all spatial layers

- `Protected.tif`: Rasterized representation of protected areas
  - Values:
    - 0: Unprotected
    - 1: Forest Reserve
    - 2: Other Protected Area
    - 3: National Park
    
These folders contain the downloaded `.tif` files used as inputs to create the cleaned and merged covariate layers.

- `00_Floodmaps`: Input floodmap `.tif` files obtained using the `floodmapr` R package. These maps capture surface water extent derived from remote sensing imagery.  

- `00_NDVI`: Input NDVI `.tif` files obtained from MODIS MOD13Q1 via Earth Engine. Each file corresponds to a specific timestamp and represents vegetation greenness.  

- `00_Panmaps`: Input `.tif` files of surface water pans derived from Sentinel-2 imagery. These are the raw classification outputs from the pan mapping algorithm before any seasonal aggregation.  

- `00_Rainmaps`: Input precipitation `.tif` files obtained from JAXA GSMaP through Earth Engine. Each file represents rainfall (mm) at a specific timestamp.  

- `00_Tempmaps`: Input air temperature `.tif` files obtained from ERA5 through Earth Engine. Each file represents temperature (°C) at a specific timestamp.  

- `00_Vegmaps`: Input vegetation cover `.tif` files derived from MODIS continuous vegetation products (MOD44B for trees and shrubs). These layers represent fractional cover for each vegetation type at each timestamp.

#### R-Data

R-Data files are in `.rds` format and can be loaded into R using `readr::read_rds()`. 

- `Covariates.rds`: Table referencing to all covariate files and indicating under what configuration they ough to be used:
  - Type: Type of the configuration:
    - Static: Covariates used for the static configuration
    - Dynamic: Covariates used for the dynamic configuration
    - DynamicAggregated: Covaraites used for dynamic predictions (for a typical year)
  - Covariate: Name of the covariate
  - Filename: Relative path pointing to the location of that covariate on the hard drive
  - Dates: List of the dates associated with that covariate
  
- `LookupTable.rds`: Lookup table used to keep track of the indices from which covariates must be extracted depending on the configuration.
  - Timestamp: Timestamp for which covariates must be extracted. This is just a vector of all potential dates that might be encountered during the study
  - Type: Type of the configuration:
    - Static: Covariates used for the static configuration
    - Dynamic: Covariates used for the dynamic configuration
    - DynamicAggregated: Covaraites used for dynmic predictions (for a typical year)
  - Covariate: Name of the covariate
  - Layerdate: Date of the layer that best aligns with the above timestamp
  - Layerindex: Index of the layer that best aligns with the above timestamp
  
- `Moonlight.rds`: Moonlight metrics as derived from the `moonlit` R-package
  - Timestamp (UTC): Timestamp for which moonlight metrics are provided (in 4-hourly intervals, matching the collected GPS sampling scheme). The time refers to the beginning of a four-hour period, except for the 07:00 timestamp, which refers to a 8-hour period
  - meanMoonPhase: average moon phase over the 4-hour period (values between 0 (new moon) and 1 (full moon)
  - meanMoonAlt: average altitude above the horizon (in degrees) of the moon
  - meanSunAlt: average altitude above the horizon (in degrees) of the sun
  - meanMoonlight: average moonlight intensity (relative to full moon)
  - meanMoonlightLux: average moonlight intensity (in lux)
  - Night: percentage of the four-hour period that is considered true night
  - LightType: Categorization of the four hour period into a light type (see Figure S5 in Appendix of associated manuscript)
    - Dark
    - Bright
    
- `SSFExtractedRedacted.rds`: Data prepared for step selection functions (iSSF). Each row represents either an observed or a random step, with all covariates extracted along the step.
  - ID: Unique identifier of the animal
  - BurstID: Identifier of the burst to which a step belongs
  - Timestamp (UTC): Timestamp at which the step was recorded (beginning of the step)
  - TimestampRounded (UTC): Timestamp rounded to the nearest hour
  - sl: Step length of the step (meters)
  - absta: Absolute turning angle (heading) of the step
  - relta: Relative turning angle of the step
  - dt: Duration of the step (hours)
  - step_id: Stratum identifier (used for conditional logistic regression)
  - inactive: Binary indicator if a step falls outside the main African wild dog activity period
  - case: Binary indicator if a step is observed (1) or random (0)
  - step_id_within: Unique identifier of steps within each stratum
  - Points: Coordinates (longitude, latitude) of interpolated points along the step (used for covariate extraction)
  - SeasonClimate: Climatic season of the step (wet vs dry)
  - SeasonHerbivores: Seasonal classification based on herbivore aggregation (concentrated vs dispersed)
  
  **Static covariates** (long-term averages):
  - HumansStatic
  - TreesStatic
  - ShrubsStatic
  - WaterStatic
  - NDVIStatic
  - DistanceToWaterStatic
  - DistanceToPansStatic
  - TemperatureStatic
  - PrecipitationStatic

  **Dynamic covariates** (time-specific):
  - HumansDynamic
  - TreesDynamic
  - ShrubsDynamic
  - WaterDynamic
  - NDVIDynamic
  - DistanceToWaterDynamic
  - DistanceToPansDynamic
  - TemperatureDynamic
  - PrecipitationDynamic

  **Moonlight covariates** (derived from the `moonlit` R-package):
  - meanMoonPhase: Average moon phase (0 = new moon, 1 = full moon)
  - meanMoonAlt: Average altitude of the moon above the horizon (degrees)
  - meanSunAlt: Average altitude of the sun above the horizon (degrees)
  - meanMoonlight: Average moonlight intensity relative to full moon
  - meanMoonlightLux: Average moonlight intensity in lux
  - Night: Percentage of the step period considered true night
  - LightType: Light category of the step period
    - Dark
    - Bright
    
- `Windows.rds`: Moving windows used to download Sentinel II data
  - Year: year of the satellite imagery to be downloaded
  - Month: Month of the satellite imagery to be downloaded
  - Window: List of the shapefiles of the moving windows
  - Tiles: List of the shapefiles of the Footprint of the Sentinel II satellite images

### 03\_Results

R-Data files are in `.rds` format and can be loaded into R using `readr::read_rds()`.  The only exception is the `SeasonalCovariates.rds` file, which has to be read using `base::load()`:

- `Formula.rds`: Contains the model formulas used for fitting the iSSF and predicting movement. Formulas represent either the simple or the full (“realistic”) model specification.

- `MovementModel.rds`: Results of the integrated step selection function (iSSF) movement models under different configurations.
  - FittingCovariates: Whether static or dynamic covariates were used for model fitting (Static vs Dynamic)
  - ModelSeasons: Whether data were split or merged across seasons (Single vs. Multi)
  - Formula: Model formula applied (simple or realistic)
  - NumberRandomSteps: Number of random steps generated per observed step (10, 25, 50, 75, 100)
  - Season: Climatic season associated with each fitted coefficient (Dry, Wet, or All)
  - Covariate: Name of the covariate/predictor
  - Coefficient: Estimated coefficient for the covariate
  - SE: Standard error of the coefficient
  - zvalue: Z statistic for testing the coefficient
  - pvalue: P-value of the coefficient
  - RandomVariance: Estimated variance of random effects
  - RandomSD: Estimated standard deviation of random effects

- `PanMapping.rds`: Validation results of the pan-mapping algorithm comparing CART and RandomForest classifiers.
  - Satellite: Satellite data used for prediction (Sentinel-2 or Landsat)
  - Data: Training or validation dataset used
  - Model: Classifier type (CART or RandomForest)
  - Varimp: Variable importance measures from the model
  - Validation: Validation dataset results
  - Confusion: List of confusion matrices for model evaluation
  - Specificity: Specificity of predictions (true negative rate)
  - Sensitivity: Sensitivity of predictions (true positive rate)
  - Accuracy: Overall classification accuracy

- `RankFrequency.rds`: Results of the spearman rank correlation validation procedure
  - FittingCovariates: Whether static or dynamic covariates were used for model fitting (Static vs Dynamic)
  - ModelSeasons: Whether data were split or merged across seasons (Single vs. Multi)
  - PredictionsCovariates: Whether static or dynamic covariates were used for predictions
  - Formula: Model formula applied (simple or realistic)
  - Replicate: Replicate of the analysis
  - Preferences: Whether predictions wer made using observed (estimated) preferences or using randomized (null model) preferences
  - Spearman: Spearman's rank correlation emerging under that configuration
  - ModelCode: Code of the model
  
- `Scaling.rds`: Scaling values used to center and standardize covariates/predictors
  - Covariate: Predictor to which scaling is applied
  - center: numeric value by which covariate was centered
  - scale: numeric value by which covariate was scaled

- `SeasonalCovariates.rds`: Results from the generalized additive models used to show the trends in Figure S1.

- `SentinelJoined.rds`: Helper dataframe used to download sentinel data:
  - Year: Year for which Sentinel II should be downloaded
  - Month: Month for which Sentinel II should be downloaded
  - Window: Spatial extent for which Sentinel II should be downloaded
  - Tiles: Tiles that match the footprint to be downloaded
  - data: Search results from Sentinel II matching the above criteria
  - NumberFiles: Number of files to be downloaded for that row

- `SentinelMetadata.rds`: Helper table with metadata of the Sentinel‑2 data that were downloaded. Metadata are extracted from the Sentinel‑2 filenames and associated XML metadata files, and describe acquisition details, processing, and tile information.
  - filepath: Local path to the downloaded `.tif` or Sentinel product file
  - name: Original filename of the downloaded Sentinel product
  - validname: Filename cleaned or standardized for internal use
  - exists: Logical indicator if the file exists on disk
  - prod_type: Product type (e.g., MSI Level‑1C, MSI Level‑2A) extracted from the filename structure or metadata
  - version: Processing baseline or version identifier of the Sentinel product
  - xml_main: Name of the main XML metadata file associated with the product
  - xml_granules: Names of the XML files for each granule/tile within the product
  - mission: Sentinel mission identifier (e.g., S2A, S2B)
  - level: Processing level (e.g., “1C” or “2A”) indicating atmospheric or surface reflectance processing
  - sensing_datetime: Timestamp of when the image was acquired (UTC)
  - id_baseline: Processing baseline identifier extracted from metadata or filename
  - tiles: Satellite tile identifier(s) (MGRS tile codes) covered by the product
  - utm: UTM zone(s) for the tile(s)
  - res: Spatial resolution(s) available in the product (e.g., 10, 20, 60 m)
  - clouds: Estimated cloud coverage percentage from product metadata
  - direction: Orbit direction of the satellite (e.g., ascending/descending)
  - orbit_n: Orbit or relative orbit number associated with the acquisition
  - preview_url: URL or local path to preview or quicklook image (if available)
  - nodata_value: Value representing no data in the raster files
  - saturated_value: Value representing saturation limits in the raster files
  - footprint: Geospatial footprint (polygon or bounding box) of the product
  - Timestamp: Date/time of the scene (usually same as `sensing_datetime`)
  - Year: Year of acquisition
  - Month: Month of acquisition
  
- `SentinelResults.rds`: Helper table summarizing the Sentinel‑2 search results for each acquisition window. Provides an overview of available imagery before downloading or processing.  
  - Year: Year of acquisition  
  - Month: Month of acquisition  
  - Window: Spatial window used for searching images
  - Tiles: MGRS tile(s) covered by the search  
  - data: Summary of available products within the window  
  - From: Earliest acquisition date in the window  
  - To: Latest acquisition date in the window  
  - Files: List of Sentinel‑2 files identified for download within the window  
  - FilesMetaData: Corresponding metadata files (XML) associated with each identified product
  
- `StepLengthDistribution.rds`: Estimated shape and scale parameter for the tentative step-length distribution:
  - shape: estimated shape parameter used for fitting iSSF and simulating random steps
  - scale: estimated scale parameter used for fitting iSSF and simulating random steps

The following files are produced when running the simulation and connectivity analyses. Intermediate results are stored in subfolders. The `.rds` files provide overviews of simulation parameters, configurations, and associated output files. For instance, `Connectivity.rds` consolidates the results stored under the `Connectivity` folder.

- `Connectivity.rds`: Overview file for connectivity analyses under different configurations.
  - Source: Source area from which simulated dispersers were released
  - Formula: Formula used to fit the iSSF model and predict movement
    - Simple: Simplistic model
    - Full: Mechanistic/Realistic model
  - ModelSeasons: Whether the model considered seasonality
    - Single: Single-season model
    - Multi: Multi-season model
  - FittingCovariates: Covariates used to fit the iSSF model
    - Static: Fixed covariates (seasonality disregarded)
    - Dynamic: Seasonal covariates
  - PredictionCovariates: Covariates used for prediction from the iSSF model
    - Static: Fixed covariates (seasonality disregarded)
    - Dynamic: Seasonal covariates
  - ModelCode: Code indicating the precise configuration
    - SSS_S: Static-Single-Static_Simple
    - SSS_F: Static-Single-Static_Full
    - DMD_S: Dynamic-Multi-Dynamic_Simple
    - DMD_F: Dynamic-Multi-Dynamic_Full

- `Simulations.rds`: Overview of all simulated dispersal trajectories under different configurations.
  - Source: Source area from which simulated dispersers were released
  - Formula: Formula used for the underlying iSSF model
  - Replicate: Simulation replicate number
  - FittingCovariates: Covariates used to fit the iSSF model
  - ModelSeasons: Seasonal configuration of the model
  - PredictionCovariates: Covariates used for prediction
  - ModelCode: Code indicating the model configuration
  - x: Simulated x-coordinate of the disperser
  - y: Simulated y-coordinate of the disperser
  - Timestamp: Simulated timestamp (UTC)
  - Filename: Path to the stored trajectory file on disk
  - Done: Logical indicator whether the simulation was completed

- `Validation.rds`: Overview of validation results for the iSSF predictions.
  - FittingCovariates: Covariates used to fit the iSSF model
  - ModelSeasons: Seasonal configuration of the model
  - PredictionCovariates: Covariates used for prediction
  - Formula: Formula used in the fitted model
  - NumberRandomSteps: Number of random steps generated per observed step for validation
  - Replicate: Validation replicate number
  - Season: Climatic season (wet vs dry)
  - ModelCode: Code indicating the model configuration
  - Filename: Path to the validation file on disk
  - Done: Logical indicator whether the validation was completed

Unconsolidated files (these files should not be directly accessed):

- `Connectivity`: Folder with connectivity results. Names are following the convention of ModelCode_SourceArea.rds
- `Simulations`: Folder with simulated trajectories. Names are following the convention of ModelCode_SourceArea_Replicate
- `Validation`: folder with validation results. Names are following the convention of ModelCode_Replicate

## Covariate Data Sources

### Landscape Characteristics

- **Trees**  
  - Type: Continuous (C)  
  - Description: Percentage tree cover  
  - Temporal resolution: 1 year  
  - Spatial resolution: 250 m  
  - Source: MODIS MOD44B  
  - Download method: RGISTools  
  - Reference: (1)

- **Shrubs / Grassland**  
  - Type: Continuous (C)  
  - Description: Percentage non-tree vegetation  
  - Temporal resolution: 1 year  
  - Spatial resolution: 250 m  
  - Source: MODIS MOD44B  
  - Download method: RGISTools  
  - Reference: (1)

- **NDVI**  
  - Type: Continuous (C)  
  - Description: Normalized Difference Vegetation Index  
  - Temporal resolution: 16 days  
  - Spatial resolution: 250 m  
  - Source: MODIS MOD13Q1  
  - Download method: rgee  
  - Reference: (2)

- **Rivers**  
  - Type: Binary (B)  
  - Description: Presence of rivers  
  - Temporal resolution: Static  
  - Spatial resolution: 90 m  
  - Source: MERIT Hydro  
  - Download method: Website  
  - Reference: (3)

- **Permanent water**  
  - Type: Binary (B)  
  - Description: Presence of permanent water bodies  
  - Temporal resolution: Static  
  - Spatial resolution: 30 m  
  - Source: Globeland30  
  - Download method: Website  
  - Reference: (4)

- **Floodwater**  
  - Type: Binary (B)  
  - Description: Presence of seasonal flood water  
  - Temporal resolution: 8 days  
  - Spatial resolution: 500 m  
  - Source: MCD43A4  
  - Download method: floodmapr  
  - Reference: (5)

- **Distance to water**  
  - Type: Continuous (C)  
  - Description: Distance to nearest river, permanent water, or seasonal flood (m)  
  - Temporal resolution: 8 days  
  - Spatial resolution: 500 m  
  - Source: Derived from water layers  
  - Download method: floodmapr  
  - Reference: (5)

- **Pans**  
  - Type: Continuous (C)  
  - Description: Presence of ephemeral pans containing water  
  - Temporal resolution: 5–10 days  
  - Spatial resolution: 10 m  
  - Source: Sentinel-2  
  - Download method: sen2r  
  - Reference: (6)

- **Distance to pans**  
  - Type: Continuous (C)  
  - Description: Distance to nearest ephemeral pan containing water (m)  
  - Temporal resolution: 5–10 days  
  - Spatial resolution: 10 m  
  - Source: Sentinel-2  
  - Download method: sen2r  
  - Reference: (6)

### Climate

- **Temperature**  
  - Type: Continuous (C)  
  - Description: 2 m above-ground temperature (°C)  
  - Temporal resolution: 4 hours  
  - Spatial resolution: 1,000 m  
  - Source: ERA5  
  - Download method: rgee  
  - Reference: (7)

- **Precipitation**  
  - Type: Continuous (C)  
  - Description: Accumulated precipitation (mm/h)  
  - Temporal resolution: 4 hours  
  - Spatial resolution: 1,000 m  
  - Source: JAXA GSMaP  
  - Download method: rgee  
  - Reference: (8)

### Anthropogenic

- **Human density**  
  - Type: Continuous (C)  
  - Description: Estimated human density (inhabitants per km²)  
  - Temporal resolution: Static  
  - Spatial resolution: 30 m  
  - Source: Facebook  
  - Download method: Website  
  - Reference: (9)

- **Agriculture**  
  - Type: Binary (B)  
  - Description: Presence of agriculture  
  - Temporal resolution: Static  
  - Spatial resolution: 30 m  
  - Source: Globeland30 / Cropland  
  - Download method: Website  
  - Reference: (4) & (10)

- **Roads**  
  - Type: Binary (B)  
  - Description: Presence of roads  
  - Temporal resolution: Static  
  - Spatial resolution: Vector data  
  - Source: OpenStreetMap  
  - Download method: Website  
  - Reference: (11)

### Light Intensity

- **Night**  
  - Type: Binary (B)  
  - Description: Binary indicator (0 = day, 1 = night)  
  - Temporal resolution: 4 hours  
  - Spatial resolution: Not applicable  
  - Source: moonlit  
  - Reference: (12)

- **Moon illumination**  
  - Type: Continuous (C)  
  - Description: Estimated moonlight illumination (lux)  
  - Temporal resolution: 4 hours  
  - Spatial resolution: Not applicable  
  - Source: moonlit  
  - Reference: (12)

---

### Notes

- Covariates were classified as either **continuous (C)** or **binary (B)**.  
- Some covariates were combined into composite predictors (e.g., water availability, human influence, light intensity).  
- Certain variables were excluded from the final model due to convergence issues.  
- Download methods listed as code-style names (e.g., `rgee`, `sen2r`, `RGISTools`) refer to R packages.

### Data Source References

(1) MODIS Vegetation Continuous Fields (MOD44B)  
NASA LP DAAC.  
https://lpdaac.usgs.gov/products/mod44bv006/

(2) MODIS NDVI (MOD13Q1)  
NASA LP DAAC.  
https://lpdaac.usgs.gov/products/mod13q1v006/

(3) MERIT Hydro  
Yamazaki, D. et al.  
https://global-hydrodynamics.github.io/MERIT_Hydro/

(4) Globeland30 Land Cover  
National Geomatics Center of China.  
http://www.globallandcover.com/

(5) MODIS Surface Reflectance (MCD43A4) / Flood mapping inputs  
NASA LP DAAC.  
https://lpdaac.usgs.gov/products/mcd43a4v006/

(Note: Floodwater layers were derived using the `floodmapr` R package from MODIS products.)

(6) Sentinel-2 Satellite Imagery  
European Space Agency (Copernicus Open Access Hub / Copernicus Data Space).  
https://dataspace.copernicus.eu/

(7) ERA5 Climate Data  
ECMWF Copernicus Climate Data Store.  
https://cds.climate.copernicus.eu/

(8) JAXA GSMaP Precipitation  
Japan Aerospace Exploration Agency (JAXA).  
https://sharaku.eorc.jaxa.jp/GSMaP/

(9) Facebook High Resolution Population Density Maps  
Meta (Facebook Data for Good).  
https://dataforgood.facebook.com/dfg/tools/high-resolution-population-density-maps

(10) Global Cropland Data
https://croplands.org/

(11) OpenStreetMap Data  
OpenStreetMap Contributors.  
https://www.openstreetmap.org/

(12) moonlit R package (moonlight calculations)  
Śmielak, M.  
https://cran.r-project.org/package=moonlit

Further details on covariate preparation and aggregation are provided in Appendix A3 and Hofmann et al. (2021).

## Sharing/Access information

If you use these data, please cite: Hofmann et al. (2026) [DOI]

Links to other publicly accessible locations of the data:

- Description of the spatial data preparation: [Hofmann et al., 2021](https://doi.org/10.1111/1365-2664.13868)
