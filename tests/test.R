library(data.table)
library(tictoc)

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
#source("R/cholesky_code.R")
source("R/summary.R")

location_formula =  math_proficiency ~ student_ses * school_ses  + (1|school_id)
scale_formula = ~ student_ses * school_ses + (1 |school_id)
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
