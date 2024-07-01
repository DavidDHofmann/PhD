This repository contains all files needed to reproduce the simulation analysis
and case study presented in the main manuscript. The repository also contains
two example analysis, showcasing application of the _dynamic+model_ approach to
simulated and real data. Finally, it contains all visualizations from the
manuscript are provided, together with the code used to produce them.

```bash
├── 02_R-Scripts                # Main folder containing all R code
│   ├── 01_SimulationAnalysis   # R codes for the simulation analysis
│   │   ├── 01_Simulation.R     # Simulation of landscapes and movement
│   │   ├── 02_Analysis.R       # iSSF analysis of simulated data
│   │   ├── 03_ValidSteps.R     # Example of how missingness influences the number of valid steps
│   │   ├── SessionInfo         # Folder with session information for the R-Scripts
│   │   └── Visualization       # Folder with codes to reproduce visualizations
│   ├── 02_CaseStudy            # R codes for the case study
│   │   ├── 01_CaseStudy.R      # R codes for the case study
│   │   ├── SessionInfo         # Folder with session information for the R-Scripts
│   │   └── Visualization       # Folder for the codes to reproduce visualizations
│   ├── 03_AppliedExamples      # Folder with coded example analyses
│   │   ├── Example_01.Rmd      # Codes for the simulation example
│   │   ├── Example_01.html     # Rendered HTML of simulation example
│   │   ├── Example_02.html     # Codes for the example with real data
│   │   └── Example_02.Rmd      # Rendered HTML for the example with real data
│   └── Functions.R             # Custom functions used throughout different R scripts
├── 03_Data
│   ├── Simulation.rds          # Overview table of simulations (links to Simulations)
│   ├── Simulations             # Folder with individual simulations
│   │   └── ...
│   ├── Analysis.rds            # Design of the simulation analysis (links to SimulationResults)
│   ├── SimulationResults       # Folder with analysis results of individual simulations
│   │   └── ...
│   ├── SimulationResultsConsolidated.rds # Consolidated and cleaned results of the simulation analysis
│   └── ValidSteps.rds          # Data for the plot showing how the number of valid steps decreases
│   ├── CaseStudyResults.rds    # File with results from the case study
│   ├── Hyena_Covariates.tif    # Spatial covariate layers for the case study (EPSG:32734)
│   ├── Hyena_GPS.rds           # GPS data for the case study
├── 04_Manuscript
│   └── ...                     # Figures shown in the manuscript
├── Glossary.csv                # Definitions of the terms in the ms
├── IrregularSSF.Rproj          # R-project file for this repository
└── README.md                   # Description of files in this repository
```
