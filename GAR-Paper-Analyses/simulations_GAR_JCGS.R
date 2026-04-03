##############################
## Simulation Script for GAR:
## - Dimension: p
## - Sample size: n
## - Edge probability: edge.prob
## - rep: Number of samples to simulate
## - model: Use as LN; Can also use L for the (unnormalized) graph Laplacian
## - Includes the GF Measure in the results
##############################


## Clear the environment
rm(list = ls())

## Call the necessary packages
library(SGM)
library(glasso)
library(doParallel)
library(doFuture)
library(future)
library(foreach)
library(mnormt)


source("GenData.R")
source("simulation-auxiliary-scripts/simulations_GAR_wrappers.R")


#####################################
#### Parallel Thread Setup
####################################
num.thread = 25
setup_simulation_parallel(num.thread)


#####################################
#### Model Setup & Data Generation
#####################################
p = 100
n = 100
rep = 100
model = "LN"
theta0 = 1
theta1 = 2

sim.setup = simulation_data_setup(p = p, n = n, rep = rep, model = model, theta0 = theta0, theta1 = theta1)
gar.setup = simulation_gar_setup(p = p, n = n)

filename = sim.setup$filename
A.tr = sim.setup$A.tr
net.tr = sim.setup$net.tr
deg = sim.setup$deg
L = sim.setup$L
LN = sim.setup$LN
L.tr = sim.setup$L.tr
theta0.tr = sim.setup$theta0.tr
theta1.tr = sim.setup$theta1.tr
v0.tr = sim.setup$v0.tr
data = sim.setup$data
Omega.tr = sim.setup$Omega.tr
Sigma.tr = sim.setup$Sigma.tr
A.ggm.tr = sim.setup$A.ggm.tr

lambda.v = gar.setup$lambda.v
rho.v = gar.setup$rho.v
net.thre = gar.setup$net.thre

summary(deg)
summary(diag(A.tr))
summary(diag(LN))
summary(LN[upper.tri(LN)])
sum(net.tr) / 2
summary(v0.tr)
sum(A.ggm.tr) / 2


######################
## Fit the Models (GAR)
######################
results.GAR = simulation_fit_gar(
  data = data,
  n = n,
  rep = rep,
  model = model,
  lambda.v = lambda.v,
  net.thre = net.thre,
  rho.v = rho.v
)


###################################################
## Model Selection and Evaluation Metrics  (GAR) ##
###################################################
gar.results = simulation_process_gar_results(
  results.GAR = results.GAR,
  rep = rep,
  net.tr = net.tr,
  L.tr = L.tr,
  theta0.tr = theta0.tr,
  theta1.tr = theta1.tr,
  v0.tr = v0.tr,
  A.ggm.tr = A.ggm.tr,
  Sigma.tr = Sigma.tr,
  Omega.tr = Omega.tr
)


##########################
# Oracle Estimator for L #
##########################
oracle.results = simulation_fit_oracle(
  data = data,
  rep = rep,
  p = p,
  n = n,
  theta0.tr = theta0.tr,
  v0.tr = v0.tr,
  net.tr = net.tr,
  model = model,
  theta1.tr = theta1.tr,
  L.tr = L.tr
)


#########################
# Fit the models (GLASSO)
#########################
glasso.fit = simulation_fit_glasso(
  data = data,
  rep = rep,
  p = p,
  n = n,
  A.ggm.tr = A.ggm.tr
)


##################################################
# Model Selection for GLASSO models via eBIC    ##
##################################################
glasso.results = simulation_process_glasso_results(
  glasso_fit = glasso.fit,
  data = data,
  rep = rep,
  p = p,
  n = n,
  Sigma.tr = Sigma.tr,
  Omega.tr = Omega.tr
)


##############################
## Results in the Main Text ##
##############################
summary.tables = simulation_results_tables(
  gar_results = gar.results,
  oracle_results = oracle.results,
  glasso_results = glasso.results
)

table1.metric = summary.tables$table1.metric
table2.metric = summary.tables$table2.metric
table3.metric = summary.tables$table3.metric

print(table1.metric)
print(table2.metric)
print(table3.metric)
