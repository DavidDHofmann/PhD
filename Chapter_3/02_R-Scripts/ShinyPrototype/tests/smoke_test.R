app_env <- new.env()
source("app.R", local = app_env)

shiny::testServer(app_env$server, {
  stopifnot(all(c(
    "Water_Start", "Water_Along", "Water_End", "Moonlight_Nearest"
  ) %in% extracted_terms()))

  session$setInputs(
    water_start = FALSE,
    water_along = TRUE,
    water_end = TRUE,
    distance_water_along = TRUE,
    moonlight_nearest = TRUE,
    temperature_nearest = FALSE
  )

  stopifnot(!("Water_Start" %in% extracted_terms()))
  stopifnot("DistanceToWater_Along" %in% extracted_terms())

  # Render predictor controls, then verify pre-formula preparation settings.
  invisible(output$predictor_settings)
  terms_now <- extracted_terms()
  water_index <- which(terms_now == "Water_Along")
  distance_index <- which(terms_now == "DistanceToWater_Along")
  stopifnot("scaled_Water_Along" %in% transformed_predictors())

  changed_settings <- list("Square root", "yes")
  names(changed_settings) <- c(
    paste0("predictor_form_", distance_index),
    paste0("predictor_scaled_", distance_index)
  )
  do.call(session$setInputs, changed_settings)
  stopifnot("scaled_sqrt_DistanceToWater_Along" %in% transformed_predictors())

  excluded_setting <- list(FALSE)
  names(excluded_setting) <- paste0("predictor_include_", water_index)
  do.call(session$setInputs, excluded_setting)
  stopifnot(!("scaled_Water_Along" %in% transformed_predictors()))

  included_setting <- list(TRUE)
  names(included_setting) <- paste0("predictor_include_", water_index)
  do.call(session$setInputs, included_setting)

  session$setInputs(model_formula_text = paste(
    "case ~ sl + log_sl +",
    "scale(sqrt(Water_Along)) +",
    "DistanceToWater_Along +",
    "Water_Along:DistanceToWater_Along + (1 | step_id)"
  ))

  stopifnot(identical(model_terms(), c("Water_Along", "DistanceToWater_Along")))

  mapping_html <- as.character(output$forecast_mapping)
stopifnot(any(grepl("WaterStatic.tif", mapping_html, fixed = TRUE)))
stopifnot(any(grepl("MODEL ALIAS", mapping_html, fixed = TRUE)))
})

cat("Shiny server test: OK\n")
