# Wildlife Movement Forecasting System — Shiny prototype

This is a sponsor-facing demonstration of the proposed Wildlife Movement
Forecasting System. It uses simulated content and does **not** run the chapter's
analytical pipeline.

## Run locally

The prototype requires R and the `shiny` package.

```r
install.packages("shiny")
shiny::runApp("ShinyPrototype", launch.browser = TRUE)
```

Run the command from the `02_R-Scripts` directory. Alternatively, open
`ShinyPrototype/app.R` in RStudio and select **Run App**.

For a fixed local address, change into `ShinyPrototype` and run:

```powershell
Rscript run.R
```

The app will be available at <http://127.0.0.1:3838>.

## What is interactive

- Navigate through five left-to-right workflow stages.
- Filter movement data and configure track regularization and iSSF preparation.
- Configure spatial start/along/end and nearest-time covariate extraction.
- See extracted terms flow into model selection and forecast-layer translation.
- Fit a demonstration model and run a demonstration simulation.
- Switch between heatmap, betweenness and interpatch outputs.
- Preview upload, covariate and export controls.

All metrics and spatial outputs are illustrative. Buttons intentionally create
short simulated jobs so sponsors can see how long-running work will be presented.
