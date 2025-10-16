library(mlmRev)
library(coda)
library(jagsUI)
library(data.table)
library(ivd)

saeb <- ivd::saeb

## Calculate school-level SES
school_ses <- saeb[, .(school_ses = mean(student_ses, na.rm = TRUE)), by = school_id]

## Join the school_ses back to the original dataset
saeb <- saeb[school_ses, on = "school_id"]

## Define student level SES as deviation from the school SES
saeb$student_ses <- saeb$student_ses - saeb$school_ses

## Grand mean center school ses
saeb$school_ses <- c(scale(saeb$school_ses, scale = FALSE))

source("R/helpers.R")
source("R/cholesky_code.R")
source("R/summary.R")

location_formula =  math_proficiency ~ 1  + (1|school_id)
scale_formula = ~ 1 + (1 |school_id)
dat <- prepare_data_for_nimble(data = saeb, location_formula = location_formula, scale_formula = scale_formula)

Kr <- ncol(dat$data$Z)
Sr <- ncol(dat$data$Z_scale)
P <- Kr + Sr

workers = 4
ss_prior_p = 0.5
niter = 2000
nburnin = 2000

## Create the bval vector for the spike-and-slab priors on scale effects
bval_scale <- matrix(rep(ss_prior_p, Sr), ncol = 1)

## Assemble the list of data and constants for JAGS
jags_data <- list(
  Y = dat$data$Y,
  X = dat$data$X,
  Z = dat$data$Z,
  X_scale = dat$data$X_scale,
  Z_scale = dat$data$Z_scale,
  N = length(dat$data$Y),
  J = dat$groups,
  K = ncol(dat$data$X),
  S = ncol(dat$data$X_scale),
  P = P,
  Kr = Kr,
  Sr = Sr,
  groupid = dat$group_id,
  bval = bval_scale
)

## Parameters to monitor
params_to_save <- c("beta", "zeta", "sigma_rand", "R", "v_final", "ss", "mu", "tau", "loglik")

## Generate the model code string
jags_model_string <- create_ss_melsm_jags_model(Kr = Kr, Sr = Sr)

# --- FIX for parallel processing: Write model to a temporary file ---
temp_model_file <- tempfile(fileext = ".txt")
writeLines(jags_model_string, temp_model_file)
on.exit(unlink(temp_model_file), add = TRUE) # Ensures cleanup

library(tictoc)
tic()
results <- jagsUI::jags(
                     data = jags_data,
                     parameters.to.save = params_to_save,
                     model.file = temp_model_file,
                     n.chains = workers,
                     n.iter = niter + nburnin,
                     n.burnin = nburnin,
                     parallel = TRUE
                   )
toc()
#results$model


out <- list()
out$samples <- results$samples
out$summary <- results$summary
out$rhat_values <- results$summary[, "Rhat"]
out$n_eff <- results$summary[, "n.eff"]

out$nimble_constants <- jags_data
out$X_location_names <- colnames(dat$data$X)
out$X_scale_names <- colnames(dat$data$X_scale)
out$Z_location_names <- colnames(dat$data$Z)
out$Z_scale_names <- colnames(dat$data$Z_scale)
out$Y <- data.frame("group_id" = dat$group_id, "Y" = dat$data$Y)
out$workers <- workers

out$Z_scale  <- dat$data$Z_scale
out$X_scale  <- dat$data$X_scale

summary.ivd(out)




## Compute logLik:
## Check that Y,  mu and tau are of same length, in case grep picks up other variables
if(length(grep("mu", colnames(combined_chains[[1]]))) != length(grep("tau", colnames(combined_chains[[1]]))) &
   length(grep("mu", colnames(combined_chains[[1]]))) != length(out$Y)) {
  stop("mu and tau are not of same lenght -- check ivd.R")
}

## Collect mu and tau
## Get mu's across chains
mu_combined <- lapply(combined_chains, function(chain) {
  mu_indices <- grep("mu", colnames(chain))
  mu_samples <- chain[, mu_indices, drop = FALSE]
  return(mu_samples)
})

## Get tau's across chains
tau_combined <- lapply(combined_chains, function(chain) {
  tau_indices <- grep("tau", colnames(chain))
  tau_samples <- chain[, tau_indices, drop = FALSE]
  return(tau_samples)
})

N <- length( dat$data$Y)
chains <- length(mu_combined)  # Number of chains
iterations <- nrow(mu_combined[[1]])  # Number of iterations (assuming all chains have same iterations)

## Initialize the array for log-likelihoods: iterations x chains x N
logLik_array <- array(NA, dim = c(iterations, chains, N))

## Loop over chains and iterations to compute log-likelihood
for (chain_idx in 1:chains) {
  for (iter in 1:iterations) {
    ## Extract mu and tau for this iteration and chain, results in vectors of length N
    mu_values <- mu_combined[[chain_idx]][iter, ]
    tau_values <- tau_combined[[chain_idx]][iter, ]

    ## Compute log-likelihood for each observation in Y
    logLik_array[iter, chain_idx, ] <- dnorm(dat$data$Y, mean = mu_values, sd = tau_values, log = TRUE)
  }
}
out$logLik_array <- logLik_array

library(loo)
loo(out$logLik_array)

#### ivd tests -----

out <- ivd(location_formula = math_proficiency ~ 1 + (1|school_id),
           scale_formula =  ~ 1 + (1|school_id),
           data = saeb,
           niter = 3000, nburnin = 5000, WAIC = TRUE, workers = 4)

ivd::summary(out)
