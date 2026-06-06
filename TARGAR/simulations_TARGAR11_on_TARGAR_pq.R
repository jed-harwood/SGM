##############################
## Simulation Script:
## TAR-GAR(p = 1, q = 1) fit on TAR-GAR(p, q) data
##
## - data.ar.order: AR order of the data-generating TAR-GAR model
## - data.q: Polynomial order of the data-generating TAR-GAR model
## - fit.ar.order: Fixed at 1
## - fit.q: Fixed at 1
## - Records 0S/eBIC selected metrics and the full eBIC path
##############################


## Clear the environment
rm(list = ls())


## Call the necessary packages
library(SGM)
library(doParallel)
library(foreach)
library(mnormt)
library(quadprog)
library(nloptr)


## Source the simulation wrappers. This lets the script run from the repository
## root or from inside the TARGAR folder.
script.dir = tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
wrapper.file = file.path(script.dir, "simulation-auxiliary-scripts",
                         "simulations_TARGAR_wrappers.R")
if (!file.exists(wrapper.file) &&
    file.exists(file.path("TARGAR", "simulation-auxiliary-scripts",
                          "simulations_TARGAR_wrappers.R"))) {
  script.dir = "TARGAR"
  wrapper.file = file.path(script.dir, "simulation-auxiliary-scripts",
                           "simulations_TARGAR_wrappers.R")
}
source(wrapper.file)


#####################################
#### Parallel Thread Setup
####################################
num.thread = 25
setup_targar_simulation_parallel(num.thread)


#####################################
#### Model Setup
#####################################
## Choose data.ar.order = 1, 2, or 3 for the data-generating model.
## The fitted model is always TAR-GAR(p = 1, q = 1).
data.ar.order = 1
config = targar_on_targar_default_config(dgp_ar_order = data.ar.order)

## User-facing simulation settings. Change these values here for experiments.
d = config$d
n = config$n
data.q = config$q
n.rep = config$n_rep
model = config$model
theta0 = config$theta0
theta1 = config$theta1
edge.prob = 2 / d
num.pass = config$num.pass

fit.ar.order = 1
fit.q = 1

## Output settings.
save.results = TRUE
keep.fits = TRUE
output.dir = file.path(script.dir, "results")

config$d = d
config$n = n
config$q = data.q
config$n_rep = n.rep
config$model = model
config$theta0 = theta0
config$theta1 = theta1
config$edge.prob = edge.prob
config$num.pass = num.pass
config$num.thread = num.thread
config$fit_ar_order = fit.ar.order
config$fit_q = fit.q
config$output_dir = output.dir
config$save_results = save.results
config$keep_fits = keep.fits


#####################################
#### Data Generation
#####################################
sim.setup = prepare_targar_simulation(config)

filename = paste0(config$case, "_d", d, "_n", n, "_model", model,
                  "_rep", n.rep, "_modular.RData")
A.tr = sim.setup$A
net.tr = sim.setup$net.tr
deg = sim.setup$deg
L.tr = sim.setup$L
theta0.tr = config$theta0
theta1.tr = config$theta1
v0.tr = sim.setup$v0
eta.tr = sim.setup$eta
eta.comp.tr = sim.setup$eta.comp
truth = sim.setup$truth
data = sim.setup$data
lambda.v = sim.setup$lambda.v
rho.v = sim.setup$rho.v
net.thre = sim.setup$net.thre

summary(deg)
sum(net.tr) / 2
summary(v0.tr)
eta.tr
c(data.ar.order = config$ar_order, data.q = config$q,
  fit.ar.order = config$fit_ar_order, fit.q = config$fit_q)


#########################
## Fit TAR-GAR(1,1)
#########################
results.TARGAR = fit_targar_replicates(
  setup = sim.setup,
  use_parallel = TRUE
)


#######################################################
## Model Selection and Evaluation Metrics  (TAR-GAR) ##
#######################################################
targar.results = extract_targar_metrics(
  results = results.TARGAR,
  setup = sim.setup
)


##############################
## Results in the Main Text ##
##############################
summary.tables = summarize_targar_metrics(targar.results)

## Selected 0S/eBIC metrics
summary.tables$L.0S.ebic.err
summary.tables$theta.0S.ebic.err
summary.tables$R1.0S.ebic.err
summary.tables$eta.0S.ebic.err
summary.tables$power.0S.ebic
summary.tables$fdr.0S.ebic
summary.tables$F1.0S.ebic
summary.tables$v0.0S.ebic
summary.tables$mean.ebic.0S.selec

## Full eBIC path, averaged across replicates
summary.tables$mean.ebic.0S.path


################
## Save Results
################
if (save.results) {
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  output.file = file.path(output.dir, filename)

  simulation.output = list(
    setup = sim.setup,
    results.TARGAR = if (keep.fits) results.TARGAR else NULL,
    targar.results = targar.results,
    summary.tables = summary.tables
  )

  save(simulation.output, file = output.file)
  message("Saved TAR-GAR-on-TAR-GAR simulation results to: ", output.file)
}
