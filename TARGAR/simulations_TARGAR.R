##############################
## Simulation Script for TAR-GAR:
## - ar.order: AR order, one of 1, 2, or 3
## - q: Polynomial order in the graph Laplacian
## - d: Number of graph nodes
## - n: Sample size
## - edge.prob: Edge probability
## - n.rep: Number of samples to simulate
## - model: Use "LN"; can also use "L" for the unnormalized graph Laplacian
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
## Manually set the simulation, as in simulations_GAR_JCGS.R.
## Change these values directly for experiments.
ar.order = 1
d = 100
n = 500
q = 1
n.rep = 100
## The original matched TAR-GAR(1,1) forecasting script generated n + 100
## rows, then fit the first n rows. This preserves that initialization.
generation.n = if (ar.order == 1 && q == 1) n + 100 else n
model = "LN"
theta0 = 1
theta1 = 2
edge.prob = 2 / d
num.pass = 3
stationary = TRUE

## Original scripts used set.seed(1) for the graph and set.seed(5*p+n)
## for data generation. Here d is the graph dimension p.
graph.seed = 1
data.seed = 5 * d + n

lambda.C = c(1.5, 1, 0.5, 0.25, 0.1)
if (d == 100) {
  C.thre = exp(seq(log(1), log(0.05), length.out = 10))
} else if (d == 250) {
  C.thre = exp(seq(log(1), log(0.075), length.out = 10))
} else {
  C.thre = exp(seq(log(1), log(0.1), length.out = 10))
}
graph.min = 0.5
graph.max = 1
selfloop = FALSE
isolate = FALSE

eps.thre = 1e-6
max.iter = 50000
deg.max.iter = 50000
lap.z.max.iter = 50000
eta.max.iter = 1000

## Build the graph first so lambda.max is available for manual eta vectors.
graph.config = list(
  d = d,
  model = model,
  edge.prob = edge.prob,
  graph_seed = graph.seed,
  graph_min = graph.min,
  graph_max = graph.max,
  selfloop = selfloop,
  isolate = isolate
)
graph.setup = prepare_targar_graph(graph.config)
lambda.max = graph.setup$lambda.max

## True TAR-GAR filter coefficients. Edit this vector directly.
eta = c(0.2, 0.7 / lambda.max)

case = paste0("TARGAR_order", ar.order, "_q", q)

## Output settings.
save.results = TRUE
keep.fits = TRUE
output.dir = file.path(script.dir, "results")

config = list(
  ar_order = ar.order,
  fit_ar_order = ar.order,
  q = q,
  fit_q = q,
  d = d,
  n = n,
  generation_n = generation.n,
  n_rep = n.rep,
  model = model,
  theta0 = theta0,
  theta1 = theta1,
  eta = eta,
  edge.prob = edge.prob,
  graph_seed = graph.seed,
  data_seed = data.seed,
  graph_min = graph.min,
  graph_max = graph.max,
  selfloop = selfloop,
  isolate = isolate,
  lambda.C = lambda.C,
  C.thre = C.thre,
  num.thread = num.thread,
  num.pass = num.pass,
  stationary = stationary,
  eps.thre = eps.thre,
  max_iter = max.iter,
  deg_max_iter = deg.max.iter,
  lap_z_max_iter = lap.z.max.iter,
  eta_max_iter = eta.max.iter,
  case = case,
  output_dir = output.dir,
  save_results = save.results,
  keep_fits = keep.fits
)


#####################################
#### Data Generation
#####################################
sim.setup = prepare_targar_simulation(config, graph_setup = graph.setup)

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
generation.n = sim.setup$generation_n
lambda.max = sim.setup$lambda.max
lambda.v = sim.setup$lambda.v
rho.v = sim.setup$rho.v
net.thre = sim.setup$net.thre

summary(deg)
sum(net.tr) / 2
summary(v0.tr)
eta.tr


#########################
## Fit the Models (TAR-GAR)
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
summary.tables$R.0S.ebic.err
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
  message("Saved TAR-GAR simulation results to: ", output.file)
}
