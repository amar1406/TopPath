# app.R — Posit Connect friendly wrapper
# This wrapper:
#   • Sets up reticulate for portable Python handling
#   • Sources your existing scripts under top_path/
#   • Launches Shiny using ui/server if defined

options(shiny.maxRequestSize = 2000*1024^2)

# ---- Reticulate: portable defaults ----
Sys.setenv(RETICULATE_AUTOCONDA = "FALSE", RETICULATE_USE_UV = "FALSE")
suppressPackageStartupMessages({
  library(shiny)
  library(reticulate)
})

# Log Python config (for debugging on server logs)
try({ print(py_config()) }, silent = TRUE)

# ---- Source your existing app scripts ----
src_files <- c(
  "top_path/util.R",
  "top_path/top_path_function_2.0.R",
  "top_path/plot_lr_wiring.R",
  "top_path/Draw_MLnet.R",
  "top_path/Run_scMLnet.R",
  "top_path/convert_h5ad_to_cds.R",
  "top_path/shiny_path.R",
  "top_path/shiny_path_3.R",
  "top_path/Shiny_Path_7.R"
)

for (f in src_files) {
  if (file.exists(f)) {
    message("Sourcing: ", f)
    source(f, local = TRUE)
  }
}

# ---- Launch the app ----
if (exists("ui", inherits = TRUE) && exists("server", inherits = TRUE)) {
  shinyApp(ui = ui, server = server)
} else if (exists("app_ui", inherits = TRUE) && exists("app_server", inherits = TRUE)) {
  shinyApp(ui = app_ui, server = app_server)
} else {
  # Fallback UI to help you debug missing ui/server definitions
  ui_fallback <- fluidPage(
    tags$h2("App wrapper is loaded"),
    tags$p("Could not find 'ui' and 'server' objects after sourcing your scripts."),
    verbatimTextOutput("pyinfo")
  )
  server_fallback <- function(input, output, session) {
    output$pyinfo <- renderPrint({
      ret <- try(py_config(), silent = TRUE)
      if (inherits(ret, "try-error")) "reticulate not ready" else ret
    })
  }
  shinyApp(ui = ui_fallback, server = server_fallback)
}
