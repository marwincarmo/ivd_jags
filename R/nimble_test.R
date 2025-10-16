library(mlmRev)
library(brms)
library(nimble)
library(coda)

# generate example data -----------------------------------------------

## Sample size
n_students = 80
n_schools = 160

## School data
N <- n_students * n_schools
school_id <- rep(1:n_schools, each = n_students)

# Random effects' diagonal SD matrix
tau <- matrix(
  c(0.3, 0,
    0, 0.1), 2, 2)

## Random effects' correlation matrix
R <- matrix(c(
  1, 0.5,
  0.5, 1), 2, 2)

## Random effects' covariance matrix
Sigma <- tau %*% R %*% t(tau)

## Random effects
u <- MASS::mvrnorm(
             n = n_schools,
             mu = c(0, 0), Sigma
           )


## Location model components
loc_rand_intc <- u[, 1][school_id]
linpred <- 0 + loc_rand_intc

 ## Scale model components
scl_rand_intc <- u[, 2][school_id]
school_sd  <- exp(0 + scl_rand_intc)

y_residuals <- rnorm(N, 0, school_sd)

y <- linpred + y_residuals

school_dat <- data.frame(
  y = y,
  school_id = school_id
)

# A correlation of 0.5 between the random effects is expected

# nimble model --------------------------------------------------------

model_code <- nimbleCode({

  for (i in 1:N){
    y[i] ~ dnorm(mu[i], sd = tau[i])
    mu[i] <- loc_int + u[groupid[i], 1]
    tau[i] <-  exp( scl_int + u[groupid[i], 2] )
  }

  ## Obtain correlated random effects
  for(j in 1:J) {
    ## normal scaling for random effects
    for( k in 1:P ){
      z[k,j] ~ dnorm(0, sd = 1)
    }
    ## Transpose U to get lower cholesky
    ## then compute the hadamard (element-wise) product with the ss vector
    u[j,1:P] <- t( sigma_rand[1:P, 1:P] %*% L[1:P, 1:P]  %*% z[1:P,j])
  }

  # Fixed intercept location
  loc_int ~ dnorm(0, sd = 10)
  # Fixed intercept scale
  scl_int ~ dnorm(0, sd = 10)
  ## Random effects SD
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 3), 0, )
  }

  ## Upper cholesky of random effects correlation
  U[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  L[1:P, 1:P] <- t(U[1:P, 1:P])

  ##
  R[1:P, 1:P] <- t(U[1:P, 1:P] ) %*% U[1:P, 1:P] # using upper cholesky
  #R[1:P, 1:P] <- L[1:P, 1:P] %*% t(L[1:P, 1:P] ) # using lower cholesky

})


constants <- list(N = nrow(school_dat), # total sample size
                   P = 2, # number of random effects
                   J = length(unique(school_dat$school_id)), # number of schools
                   groupid = school_dat$school_id) # school ids

nimble_data <- list(y = school_dat$y)

inits <- list(loc_int = rnorm(1, 5, 10),
              scl_int =  rnorm(1, 1, 3),
              sigma_rand = diag(rlnorm(constants$P, 0, 1)),
              U = diag(1,constants$P)
              )

school_model <- nimbleModel(code = model_code, name = "school_model", constants = constants,
                    data = nimble_data, inits = inits)


mcmc.out <- nimbleMCMC(code = model_code, constants = constants,
                       data = nimble_data, inits = inits,
                       nchains = 4, niter = 8000, nburnin = 2000,
                       monitors = c("loc_int", "scl_int", "sigma_rand", "L", "U", "R"),
                       summary = TRUE, WAIC = TRUE)

## Compute Rhats and n_eff:
# Convert nimbleMCMC output to an mcmc.list object
mcmc_chains <- mcmc.list(lapply(mcmc.out$samples, mcmc))

# Compute R-hat values
rhat_values <- gelman.diag(mcmc_chains, multivariate = FALSE)$psrf[, 1]

# Compute effective sample size (ESS)
ess_values <- effectiveSize(mcmc_chains)

summary_stats <- summary(mcmc_chains)

summary_table <- data.frame(
  Estimate = summary_stats$statistics[, "Mean"],
  Est.Error = summary_stats$statistics[, "SD"],
  Low95CI = summary_stats$quantiles[, "2.5%"],
  Upp95CI = summary_stats$quantiles[, "97.5%"],
  Rhat = rhat_values,
  ESS = ess_values
)

round(summary_table[!is.na(summary_table$Rhat), ],2)

# brms model --------------------------------------------------------------
#
mod0 <- brm( bf( y ~ 1 + ( 1|c|school_id),
                 sigma ~ 1 + ( 1 |c| school_id ) ),
             cores = 4, iter = 6000, warmup = 2000,
             data = school_dat,
            backend = "cmdstanr")
summary(mod0)

## Nimble cholesky (from L)
## > round(summary_table[!is.na(summary_table$Rhat), ],2)
##                  Estimate Est.Error Low95CI Upp95CI Rhat    ESS
## L[2, 1]              0.60      0.08    0.43    0.74 1.01 584.88
## L[2, 2]              0.79      0.06    0.67    0.90 1.01 619.12
## R[2, 1]              0.60      0.08    0.43    0.74 1.01 584.88
## R[1, 2]              0.60      0.08    0.43    0.74 1.01 584.88
## R[2, 2]              1.00      0.00    1.00    1.00 1.00   0.00
## U[1, 2]              0.60      0.08    0.43    0.74 1.01 584.88
## U[2, 2]              0.79      0.06    0.67    0.90 1.01 619.12
## loc_int             -0.02      0.05   -0.16    0.05 1.08 325.71
## scl_int             -0.28      0.73   -2.74    0.02 1.30 508.02
## sigma_rand[1, 1]     0.31      0.02    0.27    0.35 1.01 468.75
## sigma_rand[2, 2]     0.41      0.77    0.11    2.95 1.28 717.99

## Nimble cholesky (from U)
## > round(summary_table[!is.na(summary_table$Rhat), ],2)
##                  Estimate Est.Error Low95CI Upp95CI Rhat     ESS
## L[2, 1]              0.45      0.10    0.25    0.63 1.00 1297.14
## L[2, 2]              0.89      0.05    0.78    0.97 1.01 1299.26
## R[2, 1]              0.45      0.10    0.25    0.63 1.00 1297.14
## R[1, 2]              0.45      0.10    0.25    0.63 1.00 1297.14
## R[2, 2]              1.00      0.00    1.00    1.00 1.00    0.00
## U[1, 2]              0.45      0.10    0.25    0.63 1.00 1297.14
## U[2, 2]              0.89      0.05    0.78    0.97 1.01 1299.26
## loc_int             -0.03      0.02   -0.08    0.02 1.01  632.41
## scl_int             -0.01      0.01   -0.03    0.01 1.00 1492.79
## sigma_rand[1, 1]     0.28      0.02    0.25    0.32 1.00  668.79
## sigma_rand[2, 2]     0.09      0.01    0.08    0.11 1.00 1739.16

## brms
## ~school_id (Number of levels: 160)
##                                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sd(Intercept)                      0.28      0.02     0.25     0.32 1.00     3917     6679
## sd(sigma_Intercept)                0.09      0.01     0.08     0.11 1.00     8284     9900
## cor(Intercept,sigma_Intercept)     0.45      0.10     0.25     0.63 1.00     6814    10504

## Regression Coefficients:
##                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## Intercept          -0.03      0.02    -0.08     0.01 1.00     2812     5689
## sigma_Intercept    -0.01      0.01    -0.03     0.01 1.00     7455    11130

library(ivd)

out <- ivd(location_formula = y ~  1 + ( 1 | school_id),
           scale_formula =  ~ 1  + (1 | school_id),
           data = school_dat,
           niter = 6000, nburnin = 2000, WAIC = TRUE,
           workers = 4, n_eff = 'local')
