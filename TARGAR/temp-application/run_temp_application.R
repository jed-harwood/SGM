##############################################
## TEMP APPLICATION: USER-FACING RUN SCRIPT ##
##############################################

.script_file <- tryCatch(
  normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(.script_file)) {
  .file_arg <- commandArgs(FALSE)
  .file_arg <- .file_arg[grepl("^--file=", .file_arg)]
  if (length(.file_arg) > 0) {
    .script_file <- normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = TRUE)
  }
}
.script_dir <- if (!is.na(.script_file)) dirname(.script_file) else getwd()
.app_dir <- .script_dir

source(file.path(.app_dir, "R", "temp_application_backend.R"), chdir = TRUE)

## Pick the year and model family here when running interactively.
## Command-line arguments override these values.
year <- 2011
model <- "both"          # "both", "targar", or "graphicalvar"
save_results <- TRUE

## Optional overrides. Leave these empty for the standard temperature runs.
targar_args <- list()
graphicalvar_args <- list()
net_size_tar <- NULL

cli <- temp_parse_cli_args(
  default_year = year,
  default_model = model,
  default_save_results = save_results
)

if (isTRUE(cli$help)) {
  temp_print_cli_help()
  quit(save = "no", status = 0)
}

year <- cli$year
model <- cli$model
save_results <- cli$save_results

temp_results <- run_temp_analysis(
  year = year,
  model = model,
  save_results = save_results,
  targar_args = targar_args,
  graphicalvar_args = graphicalvar_args,
  net_size_tar = net_size_tar,
  print_summary = TRUE
)
