library(shiny)

# Wildlife Movement Forecasting System — sponsor demonstration prototype.
# All data and analytical results shown here are illustrative.

source("R/mock_data.R", local = TRUE)
addResourcePath("wfs-assets", normalizePath("www", mustWork = TRUE))

source_font <- tags$link(
  rel = "stylesheet",
  href = "https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@400;500;600;700&display=swap"
)

icon <- function(symbol) {
  tags$span(class = "ui-icon", `aria-hidden` = "true", symbol)
}

status_badge <- function(text, type = "ready") {
  tags$span(class = paste("status-badge", paste0("status-", type)), text)
}

section_header <- function(step, title, description, status = NULL) {
  div(
    class = "section-heading",
    div(
      div(class = "kicker", paste("STEP", step)),
      h1(title),
      p(description)
    ),
    if (!is.null(status)) status_badge(status, "ready")
  )
}

stage_button <- function(id, number, title, active = FALSE) {
  actionButton(
    id,
    label = tagList(
      span(class = "stage-number", number),
      span(class = "stage-title", title)
    ),
    class = paste("stage-button", if (active) "active" else "")
  )
}

summary_item <- function(label, value, detail = NULL) {
  div(
    class = "summary-item",
    span(label),
    strong(value),
    if (!is.null(detail)) tags$small(detail)
  )
}

default_model_formula <- function(terms) {
  predictors <- if (length(terms)) paste0("scaled_", terms) else character()
  paste(
    "case ~",
    paste(c("sl", "log_sl", "cos_ta", predictors, "(1 | step_id)"), collapse = " + ")
  )
}

prepared_alias <- function(term, form = "Linear", scaled = TRUE) {
  form_prefix <- switch(
    form,
    "Square root" = "sqrt_",
    "Log" = "log_",
    "Quadratic" = "quadratic_",
    "Factor" = "factor_",
    ""
  )
  scale_prefix <- if (scaled && form != "Factor") "scaled_" else ""
  paste0(scale_prefix, form_prefix, term)
}

spatial_covariates <- data.frame(
  key = c("humans", "water", "distance_water", "tree_cover"),
  file = c("Humans.tif", "Water.tif", "DistanceToWater.tif", "TreeCover.tif"),
  term = c("Humans", "Water", "DistanceToWater", "TreeCover"),
  type = c("Static", "Dynamic · 96 layers", "Dynamic · 96 layers", "Dynamic · 24 layers"),
  stringsAsFactors = FALSE
)

temporal_covariates <- data.frame(
  key = c("moonlight", "temperature", "rainfall"),
  file = c("Moonlight.csv", "Temperature.csv", "Rainfall.csv"),
  term = c("Moonlight", "Temperature", "Rainfall"),
  detail = c("Hourly timestamps", "Hourly timestamps", "Daily timestamps"),
  stringsAsFactors = FALSE
)

spatial_row <- function(key, filename, type, defaults = character()) {
  tags$tr(
    tags$td(div(class = "file-cell", icon("▧"), div(strong(filename), tags$small(type)))),
    tags$td(class = "center-cell", checkboxInput(paste0(key, "_start"), NULL, "Start" %in% defaults)),
    tags$td(class = "center-cell", checkboxInput(paste0(key, "_along"), NULL, "Along" %in% defaults)),
    tags$td(class = "center-cell", checkboxInput(paste0(key, "_end"), NULL, "End" %in% defaults))
  )
}

temporal_row <- function(key, filename, detail, selected = TRUE) {
  tags$tr(
    tags$td(div(class = "file-cell", icon("◷"), div(strong(filename), tags$small(detail)))),
    tags$td(class = "center-cell", checkboxInput(paste0(key, "_nearest"), NULL, selected)),
    tags$td(status_badge("Nearest in time", "neutral"))
  )
}

ui <- fluidPage(
  source_font,
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", type = "text/css", href = "wfs-assets/styles.css"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('setActiveStage', function(x) {
        document.querySelectorAll('.stage-button').forEach(function(el) {
          el.classList.remove('active');
        });
        var target = document.getElementById(x.id);
        if (target) target.classList.add('active');
        window.scrollTo({top: 0, behavior: 'smooth'});
      });
    "))
  ),

  div(
    class = "app-shell",
    tags$header(
      class = "topbar",
      div(
        class = "brand",
        div(class = "brand-mark", "W"),
        div(strong("Wildlife Movement Forecasting System"), span("Movement modelling and landscape connectivity"))
      ),
      div(class = "run-context", span("ACTIVE ANALYSIS"), strong("Run 2026-08-05-001"), status_badge("Demo", "demo"))
    ),

    tags$nav(
      class = "stage-nav",
      stage_button("nav_movement", "01", "Movement data", TRUE),
      div(class = "stage-connector"),
      stage_button("nav_covariates", "02", "Covariates"),
      div(class = "stage-connector"),
      stage_button("nav_model", "03", "Movement model"),
      div(class = "stage-connector"),
      stage_button("nav_forecast", "04", "Forecast"),
      div(class = "stage-connector"),
      stage_button("nav_outputs", "05", "Spatial outputs")
    ),

    tags$main(class = "main-content", uiOutput("stage_content"))
  )
)

server <- function(input, output, session) {
  active_stage <- reactiveVal("movement")
  model_fitted <- reactiveVal(FALSE)
  forecast_ready <- reactiveVal(FALSE)
  outputs_ready <- reactiveVal(FALSE)

  nav_ids <- c(
    movement = "nav_movement", covariates = "nav_covariates",
    model = "nav_model", forecast = "nav_forecast", outputs = "nav_outputs"
  )

  go_to <- function(stage) {
    active_stage(stage)
    session$sendCustomMessage("setActiveStage", list(id = nav_ids[[stage]]))
  }

  lapply(names(nav_ids), function(stage) {
    observeEvent(input[[nav_ids[[stage]]]], go_to(stage), ignoreInit = TRUE)
  })
  observeEvent(input$go_covariates, go_to("covariates"))
  observeEvent(input$go_model, go_to("model"))
  observeEvent(input$go_forecast, go_to("forecast"))
  observeEvent(input$go_outputs, go_to("outputs"))

  all_individuals <- paste0("WD", sprintf("%02d", 1:28))
  observeEvent(input$select_all_individuals, {
    updateCheckboxGroupInput(session, "animal_ids", selected = all_individuals)
  })
  observeEvent(input$clear_individuals, {
    updateCheckboxGroupInput(session, "animal_ids", selected = character())
  })

  output$individual_summary <- renderText({
    n <- length(input$animal_ids %||% all_individuals)
    if (n == length(all_individuals)) "All 28 individuals" else paste(n, "individuals selected")
  })

  flag <- function(id, default = FALSE) {
    value <- input[[id]]
    if (is.null(value)) default else isTRUE(value)
  }

  extracted_terms <- reactive({
    defaults <- list(
      humans = c(FALSE, TRUE, TRUE),
      water = c(TRUE, TRUE, TRUE),
      distance_water = c(FALSE, TRUE, TRUE),
      tree_cover = c(FALSE, FALSE, TRUE)
    )
    terms <- character()
    for (i in seq_len(nrow(spatial_covariates))) {
      key <- spatial_covariates$key[i]
      base <- spatial_covariates$term[i]
      positions <- c("Start", "Along", "End")
      for (j in seq_along(positions)) {
        if (flag(paste0(key, "_", tolower(positions[j])), defaults[[key]][j])) {
          terms <- c(terms, paste(base, positions[j], sep = "_"))
        }
      }
    }
    temporal_defaults <- c(moonlight = TRUE, temperature = TRUE, rainfall = FALSE)
    for (i in seq_len(nrow(temporal_covariates))) {
      key <- temporal_covariates$key[i]
      if (flag(paste0(key, "_nearest"), temporal_defaults[[key]])) {
        terms <- c(terms, paste0(temporal_covariates$term[i], "_Nearest"))
      }
    }
    unique(terms)
  })

  model_terms <- reactive({
    available <- extracted_terms()
    formula_text <- input$model_formula_text %||% default_model_formula(available)
    available[vapply(available, function(term) grepl(term, formula_text, fixed = TRUE), logical(1))]
  })

  output$stage_content <- renderUI({
    switch(
      active_stage(),
      movement = movement_ui(),
      covariates = covariates_ui(),
      model = model_ui(
        extracted_terms(),
        isolate(input$model_formula_text) %||% default_model_formula(extracted_terms()),
        model_fitted()
      ),
      forecast = forecast_ui(forecast_ready()),
      outputs = outputs_ui(outputs_ready())
    )
  })

  output$gps_overview <- renderPlot(draw_gps_overview(), res = 110)
  output$movement_map <- renderPlot(draw_landscape("observations"), res = 110)
  output$step_distribution <- renderPlot(draw_step_distribution(input$n_random %||% 25), res = 110)
  output$coefficient_plot <- renderPlot(draw_coefficients(model_terms()), res = 110)
  output$validation_plot <- renderPlot(draw_validation(), res = 110)
  output$forecast_map <- renderPlot(draw_landscape("simulation"), res = 110)
  output$output_map <- renderPlot(draw_connectivity(input$output_type %||% "Visitation heatmap"), res = 110)
  output$output_title <- renderText(input$output_type %||% "Visitation heatmap")

  output$ssf_artifact <- renderText({
    n <- input$n_random %||% 25
    paste(format(28472 * (n + 1), big.mark = ","), "rows")
  })

  output$extraction_summary <- renderUI({
    terms <- extracted_terms()
    tagList(
      strong(paste(length(terms), "model-ready variables")),
      span(paste(terms, collapse = " · "))
    )
  })

  output$predictor_settings <- renderUI({
    terms <- extracted_terms()
    if (!length(terms)) return(div(class = "empty-state", "Configure covariates in Step 02 first."))
    tagList(lapply(seq_along(terms), function(i) {
      div(
        class = "predictor-row",
        code(terms[i]),
        checkboxInput(paste0("predictor_include_", i), NULL, TRUE),
        selectInput(
          paste0("predictor_form_", i), NULL,
          c("Linear", "Square root", "Log", "Quadratic", "Factor"),
          width = "100%"
        ),
        selectInput(
          paste0("predictor_scaled_", i), NULL,
          c("Yes" = "yes", "No" = "no"),
          width = "100%"
        )
      )
    }))
  })

  transformed_predictors <- reactive({
    terms <- extracted_terms()
    prepared <- lapply(seq_along(terms), function(i) {
      include <- isTRUE(input[[paste0("predictor_include_", i)]] %||% TRUE)
      form <- input[[paste0("predictor_form_", i)]] %||% "Linear"
      scaled <- identical(input[[paste0("predictor_scaled_", i)]] %||% "yes", "yes")
      if (!include) return(NULL)
      prepared_alias(terms[i], form, scaled)
    })
    unlist(prepared, use.names = FALSE)
  })

  prepared_predictor_map <- reactive({
    terms <- extracted_terms()
    rows <- lapply(seq_along(terms), function(i) {
      form <- input[[paste0("predictor_form_", i)]] %||% "Linear"
      scaled <- identical(input[[paste0("predictor_scaled_", i)]] %||% "yes", "yes")
      data.frame(
        raw_term = terms[i],
        prepared_term = prepared_alias(terms[i], form, scaled),
        included = isTRUE(input[[paste0("predictor_include_", i)]] %||% TRUE),
        form = form,
        scaled = scaled && form != "Factor",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })

  observeEvent(input$build_formula, {
    formula_text <- paste(
      "case ~",
      paste(c("sl", "log_sl", "cos_ta", transformed_predictors(), "(1 | step_id)"), collapse = " + ")
    )
    updateTextAreaInput(session, "model_formula_text", value = formula_text)
  })

  output$forecast_mapping <- renderUI({
    terms <- model_terms()
    if (!length(terms)) return(div(class = "empty-state", "Select model covariates in Step 03 first."))
    rows <- lapply(seq_along(terms), function(i) {
      term <- terms[i]
      base <- sub("_(Start|Along|End|Nearest)$", "", term)
      temporal <- grepl("_Nearest$", term)
      fitted_source <- switch(base,
        Water = "Water.tif · dynamic",
        DistanceToWater = "DistanceToWater.tif · dynamic",
        Humans = "Humans.tif · static",
        TreeCover = "TreeCover.tif · dynamic",
        Moonlight = "Moonlight.csv",
        Temperature = "Temperature.csv",
        Rainfall = "Rainfall.csv",
        paste0(base, ".tif")
      )
      choices <- if (temporal) {
        c(paste0(base, ".csv"), paste0(base, "_Scenario.csv"))
      } else if (base %in% c("Water", "DistanceToWater", "TreeCover")) {
        c(paste0(base, "Static.tif"), paste0(base, "Dynamic.tif"), "Upload replacement…")
      } else {
        c(paste0(base, ".tif"), "Upload replacement…")
      }
      div(
        class = "mapping-row",
        div(strong(term), tags$small(paste("Fitted from", fitted_source))),
        div(class = "mapping-arrow", "→"),
        selectInput(paste0("forecast_source_", i), NULL, choices, width = "100%"),
        div(
          class = "alias-cell",
          tags$small("PREPARED MODEL ALIAS"),
          code({
            mapping <- prepared_predictor_map()
            alias <- mapping$prepared_term[mapping$raw_term == term]
            if (length(alias)) alias[1] else term
          })
        )
      )
    })
    tagList(rows)
  })

  observeEvent(input$fit_model, {
    withProgress(message = "Fitting demonstration model", value = 0, {
      for (i in seq_len(5)) {
        Sys.sleep(0.16)
        incProgress(0.2, detail = c("Preparing strata", "Optimizing likelihood", "Checking convergence", "Validating", "Saving model")[i])
      }
    })
    model_fitted(TRUE)
    showNotification("Demonstration model fitted successfully.", type = "message")
  })

  observeEvent(input$run_forecast, {
    withProgress(message = "Running forecast", value = 0, {
      for (i in seq_len(6)) {
        Sys.sleep(0.14)
        incProgress(1 / 6, detail = paste("Trajectory batch", i, "of 6"))
      }
    })
    forecast_ready(TRUE)
    showNotification("500 demonstration trajectories generated.", type = "message")
  })

  observeEvent(input$derive_outputs, {
    withProgress(message = "Deriving spatial outputs", value = 0, {
      for (i in seq_len(4)) {
        Sys.sleep(0.16)
        incProgress(0.25)
      }
    })
    outputs_ready(TRUE)
    showNotification("Spatial outputs are ready.", type = "message")
  })
}

movement_ui <- function() {
  tagList(
    section_header("01", "Movement data", "Import standardized tracking data, focus the analysis and prepare observed and available movement steps.", "Ready"),
    div(
      class = "content-card import-bar",
      div(class = "import-copy", h3("Import tracking data"), p("Movebank-style CSV or RDS with individual ID, timestamp, location and optional behavioural fields.")),
      fileInput("movement_file", NULL, accept = c(".csv", ".rds", ".parquet"), buttonLabel = "Choose movement file"),
      status_badge("184,392 fixes loaded", "ready")
    ),
    div(
      class = "content-card",
      div(class = "card-heading", h3("Filter observations"), p("Filters define the data used for regularization and model fitting.")),
      div(
        class = "control-grid four",
        dateRangeInput("movement_dates", "Date range", start = "2016-03-01", end = "2025-08-01"),
        div(
          class = "individual-field",
          tags$label("Individuals"),
          tags$details(
            class = "individual-dropdown",
            tags$summary(textOutput("individual_summary", inline = TRUE)),
            div(
              class = "individual-menu",
              div(class = "individual-actions", actionLink("select_all_individuals", "Select all"), actionLink("clear_individuals", "Clear")),
              checkboxGroupInput("animal_ids", NULL, choices = paste0("WD", sprintf("%02d", 1:28)), selected = paste0("WD", sprintf("%02d", 1:28)))
            )
          )
        ),
        selectInput("behaviour", "Behavioural mode", c("All modes", "Dispersal", "Resident", "Exploratory"), selected = "Dispersal"),
        selectInput("sex_filter", "Sex", c("All", "Female", "Male"))
      )
    ),
    div(class = "summary-strip", summary_item("Individuals", "28"), summary_item("Valid fixes", "184,392", "97.8% retained"), summary_item("Date coverage", "9.4 years"), summary_item("Median interval", "4.0 hours")),
    div(
      class = "two-column diagnostics-layout",
      div(class = "content-card plot-card", div(class = "card-heading", h3("GPS overview"), p("Diagnostics used in the chapter: fix totals, collection hours and successive-fix time lags.")), plotOutput("gps_overview", height = "340px")),
      div(class = "content-card plot-card", div(class = "card-heading", h3("Filtered tracks"), p("Spatial extent after the current filters.")), plotOutput("movement_map", height = "340px"))
    ),
    div(class = "subsection-heading", div(h2("Regularization and iSSF preparation"), p("Standardize fix schedules, create bursts and generate available steps."))),
    div(
      class = "two-column settings-layout",
      div(
        class = "content-card",
        div(class = "control-grid two",
          sliderInput("fix_rate", "Target fix interval", 1, 12, 4, step = 1, post = " hours"),
          sliderInput("tolerance", "Temporal tolerance", .1, 2, .5, step = .1, post = " hours"),
          numericInput("min_burst", "Minimum fixes per burst", 3, min = 3, max = 20),
          selectInput("step_distribution", "Step-length distribution", c("Gamma", "Log-normal", "Empirical")),
          sliderInput("n_random", "Available steps per observed step", 10, 100, 25, step = 5),
          numericInput("seed_steps", "Random seed", 20260805)
        ),
        div(class = "action-row", actionButton("prepare_steps", "Regularize and generate steps", class = "primary-action"), div(class = "artifact-label", icon("◇"), div(strong("SSF.parquet"), span(textOutput("ssf_artifact", inline = TRUE)))))
      ),
      div(class = "content-card plot-card", div(class = "card-heading", h3("Step-length distribution"), p("Observed steps and fitted availability distribution.")), plotOutput("step_distribution", height = "270px"))
    ),
    div(class = "next-row", actionButton("go_covariates", "Continue to covariates  →", class = "primary-action"))
  )
}

covariates_ui <- function() {
  tagList(
    section_header("02", "Covariates", "Configure where spatial values are sampled and associate time-indexed environmental data.", "Configured"),
    div(
      class = "content-card",
      div(class = "card-heading heading-with-action", div(h3("Spatial raster covariates"), p("A raster may contribute separate start, along-step and end-point variables.")), actionButton("add_spatial", "+ Add raster", class = "secondary-action")),
      div(class = "table-wrap", tags$table(
        class = "config-table",
        tags$thead(tags$tr(tags$th("Raster and structure"), tags$th("Extract at start"), tags$th("Extract along"), tags$th("Extract at end"))),
        tags$tbody(
          spatial_row("humans", "Humans.tif", "Static · 1 layer", c("Along", "End")),
          spatial_row("water", "Water.tif", "Dynamic · 96 layers", c("Start", "Along", "End")),
          spatial_row("distance_water", "DistanceToWater.tif", "Dynamic · 96 layers", c("Along", "End")),
          spatial_row("tree_cover", "TreeCover.tif", "Dynamic · 24 layers", c("End"))
        )
      ))
    ),
    div(
      class = "content-card along-settings",
      div(class = "card-heading", h3("Along-step extraction"), p("Points are interpolated along each step before their values are summarized.")),
      div(class = "control-grid three",
        numericInput("interpolation_distance", "Interpolation distance", 250, min = 10, step = 10),
        selectInput("interpolation_units", "Distance unit", c("metres", "kilometres")),
        selectInput("along_function", "Combine values using", c("Mean", "Minimum", "Maximum", "Median", "Sum", "Standard deviation"))
      )
    ),
    div(
      class = "content-card",
      div(class = "card-heading heading-with-action", div(h3("Temporal covariates"), p("Timestamped CSV data are joined to each step using the nearest available observation.")), actionButton("add_temporal", "+ Add time series", class = "secondary-action")),
      div(class = "table-wrap", tags$table(
        class = "config-table temporal-table",
        tags$thead(tags$tr(tags$th("Time series"), tags$th("Use"), tags$th("Association"))),
        tags$tbody(
          temporal_row("moonlight", "Moonlight.csv", "Hourly timestamps", TRUE),
          temporal_row("temperature", "Temperature.csv", "Hourly timestamps", TRUE),
          temporal_row("rainfall", "Rainfall.csv", "Daily timestamps", FALSE)
        )
      )),
      div(class = "control-grid two temporal-options", numericInput("nearest_tolerance", "Maximum time difference", 2, min = 0), selectInput("nearest_units", "Time unit", c("hours", "days", "minutes")))
    ),
    div(class = "artifact-strip", div(icon("◇"), uiOutput("extraction_summary", inline = TRUE)), actionButton("extract_covariates", "Extract configured variables", class = "primary-action")),
    div(class = "next-row", actionButton("go_model", "Continue to movement model  →", class = "primary-action"))
  )
}

model_ui <- function(choices, formula_value, fitted) {
  tagList(
    section_header("03", "Movement model", "Specify the iSSF from variables created in Step 02. Available terms update automatically with the extraction configuration.", if (fitted) "Model fitted" else "Ready to fit"),
    div(
      class = "two-column model-layout",
      div(
        class = "content-card",
        div(class = "card-heading", h3("Predictor preparation"), p("Set the functional form and scaling used when a predictor is added to the formula.")),
        div(class = "predictor-header", span("PREDICTOR"), span("INCLUDE"), span("FUNCTIONAL FORM"), span("SCALED")),
        uiOutput("predictor_settings"),
        actionButton("build_formula", "Prepare predictors and build formula", class = "secondary-action full-width"),
        div(
          class = "scaling-artifact",
          icon("◇"),
          div(
            strong("PredictorScaling.rds"),
            span("Transformations, centres and scales saved for forecasting")
          )
        ),
        selectInput("model_random_steps", "Available steps used per stratum", c(10, 25, 50, 75, 100), selected = 25),
        textAreaInput(
          "model_formula_text", "Model formula",
          value = formula_value,
          rows = 7,
          placeholder = "case ~ sl + log_sl + scaled_sqrt_DistanceToWater_Along + Water_End:Humans_Along + (1 | step_id)"
        ),
        p(class = "formula-help", "The formula uses prepared columns. Write or edit it directly to add interactions or specify random effects."),
        actionButton("fit_model", if (fitted) "Refit demonstration model" else "Fit demonstration model", class = "primary-action full-width")
      ),
      div(
        class = "model-results",
        div(class = "summary-strip compact", summary_item("Strata", "28,472"), summary_item("Individuals", "28"), summary_item("Convergence", if (fitted) "Successful" else "Not fitted")),
        div(class = "content-card plot-card", div(class = "card-heading", h3("Selection coefficients"), p("Demonstration estimates with 95% confidence intervals.")), plotOutput("coefficient_plot", height = "310px")),
        div(class = "content-card plot-card", div(class = "card-heading", h3("Predictive validation"), p("Used-step ranking across five validation folds.")), plotOutput("validation_plot", height = "220px"))
      )
    ),
    div(class = "next-row", actionButton("go_forecast", "Continue to forecast  →", class = "primary-action"))
  )
}

forecast_ui <- function(ready) {
  tagList(
    section_header("04", "Forecast", "Apply the fitted model to a new landscape by translating forecast layers into the covariate names expected by the model.", if (ready) "Forecast complete" else "Configure landscape"),
    div(
      class = "content-card landscape-choice",
      div(class = "card-heading", h3("Forecast landscape"), p("Use the fitting covariates by default, or provide replacement layers for a new landscape or scenario.")),
      checkboxInput("replace_landscape", "Use different covariate sources for this forecast", FALSE),
      conditionalPanel(
        condition = "input.replace_landscape === false",
        div(class = "same-landscape-note", icon("✓"), "Fitting covariate sources will be reused; no translation is required.")
      ),
      conditionalPanel(
        condition = "input.replace_landscape === true",
        div(
          class = "mapping-table",
          div(class = "mapping-header", span("MODEL VARIABLE / FITTING SOURCE"), span(""), span("FORECAST SOURCE"), span("PASSED TO MODEL AS")),
          uiOutput("forecast_mapping")
        )
      )
    ),
    div(
      class = "two-column forecast-layout",
      div(
        class = "content-card",
        div(class = "card-heading", h3("Forecast settings"), p("Define starting conditions and simulation effort.")),
        selectInput("source_area", "Start locations", c("Upload start locations…", "Draw on map…", "SourceAreas.gpkg"), selected = "SourceAreas.gpkg"),
        dateInput("start_date", "Start date", value = "2026-09-01"),
        sliderInput("n_days", "Duration", 14, 365, 90, step = 7, post = " days"),
        sliderInput("n_individuals", "Simulated trajectories", 50, 2000, 500, step = 50),
        numericInput("seed_forecast", "Random seed", 24051987),
        actionButton("run_forecast", if (ready) "Rerun demonstration forecast" else "Run demonstration forecast", class = "primary-action full-width")
      ),
      div(class = "content-card plot-card", div(class = "card-heading", h3("Forecast trajectories"), p("Preview of simulated movement from the selected start locations.")), plotOutput("forecast_map", height = "460px"))
    ),
    div(class = "next-row", actionButton("go_outputs", "Continue to spatial outputs  →", class = "primary-action"))
  )
}

outputs_ui <- function(ready) {
  tagList(
    section_header("05", "Spatial outputs", "Translate simulated trajectories into decision-ready landscape and connectivity summaries.", if (ready) "Outputs ready" else "Ready to derive"),
    div(
      class = "output-toolbar",
      radioButtons("output_type", NULL, c("Visitation heatmap", "Betweenness", "Interpatch connectivity"), selected = "Visitation heatmap", inline = TRUE),
      actionButton("derive_outputs", if (ready) "Refresh outputs" else "Generate outputs", class = "primary-action")
    ),
    div(
      class = "two-column output-layout",
      div(class = "content-card plot-card output-map", div(class = "card-heading heading-with-action", div(h3(textOutput("output_title", inline = TRUE)), p("Illustrative result from 500 trajectories over 90 days.")), status_badge("Dummy data", "demo")), plotOutput("output_map", height = "560px")),
      div(
        class = "output-side",
        div(class = "content-card", h3("Summary"), div(class = "vertical-summary", summary_item("Landscape reached", "61.4%"), summary_item("Patch connections", "17 of 24"), summary_item("Median arrival", "38 days"))),
        div(class = "content-card", h3("Export"), p(class = "muted", "Package maps, tables and the reproducibility manifest."), checkboxInput("export_geotiff", "GeoTIFF surfaces", TRUE), checkboxInput("export_gpkg", "GeoPackage vectors", TRUE), checkboxInput("export_report", "Methods report", TRUE), actionButton("download_package", "Download output package", class = "secondary-action full-width"))
      )
    )
  )
}

shinyApp(ui, server)
