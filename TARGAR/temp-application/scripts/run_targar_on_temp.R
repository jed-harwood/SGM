####################################
## Run the TAR-GAR training script ##
####################################

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
.app_dir <- normalizePath(file.path(.script_dir, ".."), mustWork = TRUE)

source(file.path(.app_dir, "R", "temp_application_backend.R"), chdir = TRUE)

cli <- temp_parse_cli_args(default_year = 2011, default_model = "targar")

if (isTRUE(cli$help)) {
  temp_print_cli_help()
  quit(save = "no", status = 0)
}

temp_results <- run_temp_analysis(
  year = cli$year,
  model = "targar",
  save_results = cli$save_results,
  print_summary = TRUE
)

results_tar <- temp_results$targar
