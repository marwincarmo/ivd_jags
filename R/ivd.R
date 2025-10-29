library(nimble)
library(data.table)
library(future.apply)
library(coda)

load("data/saeb.rda")
source("R/helpers.R")

school_ses <- saeb[, .(school_ses = mean(student_ses, na.rm = TRUE)), by = school_id]

## Join the school_ses back to the original dataset
saeb <- saeb[school_ses, on = "school_id"]

## Define student level SES as deviation from the school SES
saeb$student_ses <- saeb$student_ses - saeb$school_ses

## Grand mean center school ses
saeb$school_ses <- c(scale(saeb$school_ses, scale = FALSE))

## Define model
location_formula = math_proficiency ~  1 + ( 1 | school_id)
scale_formula =  ~ 1  + (1 | school_id)
data = saeb
niter = 1000
nburnin = 500
WAIC = TRUE
workers = 4
n_eff = 'local'
ss_prior_p = 0.5

if(is.null(nburnin)) {
  nburnin <- niter
}
niter <- niter + nburnin
dat <- prepare_data_for_nimble(data = data, location_formula = location_formula, scale_formula = scale_formula)
data <- dat[[1]]
groups <- dat$groups
group_id <- dat$group_id

  ## Nimble part:
  ## Nimble constants
constants <- list(N = length(data$Y),
                  J = groups,
                  K = ncol(data$X),  ## number of fixed location effects
                  Kr = ncol(data$Z), ## number of random location effects
                  S = ncol(data$X_scale),  ## number of fixed scale effects
                  Sr = ncol(data$Z_scale),  ## number of random scale effects
                  P = ncol(data$Z) + ncol(data$Z_scale),  ## number of random effects
                  groupid = group_id,
                  bval = matrix(c(rep(1,  ncol(data$Z)), rep(ss_prior_p, ncol(data$Z_scale)) ), ncol = 1)) ## Prior probability for dbern
  ## Nimble inits
inits <- list(beta = rnorm(constants$K, 5, 10), ## TODO: Check inits
              zeta =  rnorm(constants$S, 1, 3),
              sigma_rand = diag(rlnorm(constants$P, 0, 1)),
              L = diag(1,constants$P),
              zscore = rnorm(constants$P))

modelCode <- nimbleCode({
  ## Likelihood components:
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], sd = tau[i]) ## explicitly ask for SD not precision
    ## Check if K (number of fixed location effects) an S (number of fixed scale effecs)
    ## are greater than 1, if not, use simplified computation to avoid indexing issues in nimble
    ## Location
    ## Check if we have more than just an intercept:
    if(K>1) {
      if(Kr>1) {
        mu[i] <- sum(beta[1:K] * X[i, 1:K]) + sum( u[groupid[i], 1:Kr] * Z[i, 1:Kr] )
      } else {
        mu[i] <- sum(beta[1:K] * X[i, 1:K]) + u[groupid[i], 1]
      }
    } else {
      mu[i] <- beta[1] + u[groupid[i], 1] * Z[i, 1]
    }

    ## Scale
    ## Check if we have more than just an fixed intercept:
    if(S>1) {
      if(Sr>1) {
        tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + sum(u[groupid[i], (Kr+1):(Kr+Sr)] * Z_scale[i, 1:Sr]) )
      } else {
        tau[i] <- exp( sum(zeta[1:S] * X_scale[i, 1:S]) + u[groupid[i], (Kr+1)] )
      }
    } else {
      ## This assumes that if there is only one fixed intercept in scale, there is also exactly one random intercept in scale,
      ## and no other effects
      tau[i] <- exp( zeta[1] + u[groupid[i], (Kr+1)] )
    }
  }
  ## Obtain correlated random effects
  for(j in 1:J) {
    ## Bernoulli for Spike and Slab
    for(p in 1:P){
      ss[p,j] ~ dbern(bval[p,1]) ## bval is a constant
    }
    ## normal scaling for random effects
    for( k in 1:P ){
      z[k,j] ~ dnorm(0, sd = 1)
    }
    ## Transpose L to get lower cholesky
    ## then compute the hadamard (element-wise) product with the ss vector
    u[j,1:P] <- t( sigma_rand[1:P, 1:P] %*% L[1:P, 1:P]  %*% z[1:P,j] * ss[1:P,j] )
  }
  ## Priors:
  ## Fixed effects: Location
  for (k in 1:K) {
    beta[k] ~ dnorm(0, sd = 1000)
  }
  ## Fixed effects: Scale
  for (s in 1:S) {
    zeta[s] ~ dnorm(0, sd = 3)
  }
  ## Random effects SD
  for(p in 1:P){
    sigma_rand[p,p] ~ T(dt(0, 1, 1), 0, )  # Half-cauchy with 1,1
  }

  ## Correlations between random effects
  for(p in 1:P){
    zscore[p] ~ dnorm(0, sd = 1)
    rho[p] <- tanh(zscore[p])
  }

  ## Lower cholesky of random effects correlation
  L[1:P, 1:P] <- invvec_to_corr(V = rho[1:P], P = P)

                                        #L[1:P, 1:P] ~ dlkj_corr_cholesky(eta = 1, p = P)
  ##
  R[1:P, 1:P] <- L[1:P, 1:P]  %*% t(L[1:P, 1:P])
})


future::plan(multisession, workers = workers)

results <- future_lapply(1:workers, function(x) {
  compiled_model <- build_ivd_model(
    code = modelCode,
    constants = constants,
    dummy_data = data,
    dummy_inits = inits,
    useWAIC = WAIC
  )
  run_MCMC_compiled_model(
    compiled = compiled_model,
    seed = x,
    new_data = data,
    new_inits = inits,
    niter = niter,
    nburnin = nburnin,
    useWAIC = WAIC)},

  future.seed = TRUE,
  future.packages = c("nimble"),
  future.globals = list(
    modelCode = modelCode,
    constants = constants,
    data = data,
    inits = inits,
    niter = niter,
    nburnin = nburnin,
    WAIC = WAIC,
    build_ivd_model = build_ivd_model,
    run_MCMC_compiled_model = run_MCMC_compiled_model,
    invvec_to_corr = invvec_to_corr))

## Prepare object to be returned
out <- list()

## When WAIC=TRUE, 'results' is a list of lists. When WAIC=FALSE, it's a list of matrices.
## This block normalizes it so 'mcmc_chains' is always a list of mcmc objects.
if (WAIC) {
  ## Store WAIC results separately
  out$WAIC_results <- lapply(results, `[[`, "WAIC")
  ## Extract just the samples for the mcmc list
  mcmc_chains <- lapply(results, function(x) coda::as.mcmc(x$samples))
} else {
  mcmc_chains <- lapply(results, coda::as.mcmc)
}
mcmc_chains0 <- lapply(results0, as.mcmc)
combined_chains0 <- coda::mcmc.list(mcmc_chains0)

## Compute logLik:
## Check that Y,  mu and tau are of same length, in case grep picks up other variables
## This check now works correctly because combined_chains[[1]] is always an mcmc object.
if(length(grep("mu", colnames(combined_chains[[1]]))) != length(grep("tau", colnames(combined_chains[[1]]))) &
   length(grep("mu", colnames(combined_chains[[1]]))) != length(data$Y)) {
  stop("mu and tau are not of same lenght -- check ivd.R")
}

## Collect mu and tau for logLik calculation
mu_combined <- lapply(combined_chains, function(chain) {
  mu_indices <- grep("mu", colnames(chain))
  chain[, mu_indices, drop = FALSE]
})

tau_combined <- lapply(combined_chains, function(chain) {
  tau_indices <- grep("tau", colnames(chain))
  chain[, tau_indices, drop = FALSE]
})

N <- length(data$Y)
chains <- length(mu_combined)
iterations <- nrow(mu_combined[[1]])

## Initialize and fill the array for log-likelihoods
if (WAIC) {
  logLik_array <- array(NA, dim = c(iterations, chains, N))
  for (chain_idx in 1:chains) {
    for (iter in 1:iterations) {
      mu_values <- mu_combined[[chain_idx]][iter, ]
      tau_values <- tau_combined[[chain_idx]][iter, ]
      logLik_array[iter, chain_idx, ] <- dnorm(data$Y, mean = mu_values, sd = tau_values, log = TRUE)
    }
  }
  out$logLik_array <- logLik_array
}

## Compute Rhats and n_eff:
## x <- combined_chains # This is now consistently the correct object
parameters <- ncol(combined_chains[[1]])
samples_array <- array(NA, dim = c(iterations, chains, parameters))



## Use the monitor function from rstan to obtain Rhat (coda's gelman.rhat does not work reliably)
print("Compiling results...")


## Split Rhat and split n_eff:
## Vehtari et al doi:10.1214/20-BA1221 available at
## http://www.stat.columbia.edu/~gelman/research/published/rhat.pdf

## Split each chain into two halves for split-Rhat calculation
split_chains <- list()
for (i in 1:chains) {
  split_chains[[(i * 2) - 1]] <- combined_chains[[i]][1:(iterations / 2), ]
  split_chains[[i * 2]] <- combined_chains[[i]][(iterations / 2 + 1):iterations, ]
}

m <- length(split_chains) # Number of split chains
n <- nrow(split_chains[[1]]) # Length of each split chain

## Calculate B and W iteratively
chain_means <- sapply(split_chains, colMeans)
grand_mean <- rowMeans(chain_means)

B <- n * apply(chain_means, 1, var) #n * rowSums((t(chain_means) - grand_mean)^2) / (m - 1)

chain_variances <- sapply(split_chains, function(chain) apply(chain, 2, var))
W <- rowMeans(chain_variances)

## Calculate R-hat
vtp <- (n - 1) * W / n + B / n
Rhat <- sqrt(vtp / W)

## out$rhat_values <- Rhat
## if( any(out$rhat_values[!is.na(out$rhat_values)] > 1.1) ) warning("Some R-hat values are greater than 1.10 -- increase warmup and/or sampling iterations." )

## ## Initialize a new array with double the chains, half the iterations
## split_samples <- array(NA, dim = c(iterations / 2, chains * 2, parameters))

## ## Fill the 3D array with the data from the list
## for (i in seq_along(combined_chains)) {
##   samples_array[, i, ] <- combined_chains[[i]]
## }

## ## Split each chain into two halves
## for (c in 1:chains) {
##   ## First half of the iterations for the first split chain
##   split_samples[, (c * 2) - 1, ] <- samples_array[1:(iterations / 2), c, ]
##   ## Second half of the iterations for the second split chain
##   split_samples[, c * 2, ] <- samples_array[(iterations / 2 + 1):iterations, c, ]
## }

## ## m - number of chains after splitting
## m2 <- dim(split_samples )[2]
## ## n be the length of th chain
## n2 <- dim(split_samples )[1]
## s2 <- m2*n2

## ## B ingredients
## tdm <- apply(split_samples, 3, function(slice) {
##   apply(slice, 2, mean)
## })
## tdd <- colMeans(tdm)

## ## Eq. 3.1
## result <- tdm - matrix(tdd, nrow = nrow(tdm), ncol = ncol(tdm), byrow = TRUE)
## B2 <- apply(result, 2, function(x) sum(x^2)*n/(m-1))

## ## Eq. 3.2
## W <- apply(split_samples, 3, function(param_samples) {
##   chain_variances <- apply(param_samples, 2, function(chain) {
##     chain_mean <- mean(chain)
##     sum((chain - chain_mean)^2) / (n - 1)  # s2m: Variance for each chain
##   })
##   mean(chain_variances)  # Average variance across all split chains
## })

## ## Eq. 3.3
## vtp <- (n-1)*W/n + B/n

## ## Compute R-hat
## Rhat2 <- sqrt(vtp / W)

if (n_eff == 'local') {
  ## Compute split-chain n_eff
  ## ACF is computed using FFT as per Vehtari et al.
  ## Compute for multiple chains, following eq 10:

  rho_t_list <- lapply(1:parameters, function(p) {

    ## This returns a list of vectors of potentially different lengths.
    s2m_rtm_list <- lapply(split_chains, function(chain) {
      var_chain <- var(chain[, p])
      if (is.na(var_chain) || var_chain == 0) {
        return(numeric(0)) # Return empty vector for constant chains
      }
      acf_values <- .autocorrelation_fft(chain[, p])
      position <- min(which(acf_values[-length(acf_values)] + acf_values[-1] < 0), na.rm = TRUE)

      rho <- if (is.finite(position)) {
               acf_values[2:(position + 1)]
             } else {
               numeric(0) # No valid position found
             }
      return(var_chain * rho)
    })

    ## Find the maximum length of the calculated ACF vectors
    max_len <- max(sapply(s2m_rtm_list, length))

    ## If all chains for this parameter were constant, return a vector of zeros
    if (max_len == 0) {
      return(rep(0, n - 1))
    }

    ## Pad the vectors with NA and bind them into a matrix
    s2m_rtm_matrix <- do.call(cbind, lapply(s2m_rtm_list, function(x) {
      c(x, rep(NA, max_len - length(x)))
    }))

    ## Now rowMeans will work correctly
    avg_s2m_rtm <- rowMeans(s2m_rtm_matrix, na.rm = TRUE)

    ## The final rho_t needs to be of length n-1. Pad with zeros.
    full_avg_s2m_rtm <- c(avg_s2m_rtm, rep(0, (n - 1) - length(avg_s2m_rtm)))

    rho_t <- 1 - (W[p] - full_avg_s2m_rtm) / vtp[p]
    return(rho_t)
  })

  ## Calculate effective sample size
  n_eff <- m * n / (1 + 2 * sapply(rho_t_list, sum, na.rm = TRUE))
  ## n_eff cannot be larger than the total number of samples
  n_eff <- pmin(n_eff, n * m)
  n_eff <- round(n_eff)

#########----------------------------------------------------------

  mn_s2m_ptm <- apply(split_samples, 3, function(param_samples) {

    chain_variances <- apply(param_samples, 2, function(chain) {
      chain_mean <- mean(chain)
      sum((chain - chain_mean)^2) / (n - 1)  # s2m: Variance for each chain
    })

    chain_rho <- apply(param_samples, 2, function(samp_per_chain) {
      acf_values <- .autocorrelation_fft( samp_per_chain )
      ## Truncate according to Geyer (1992)
      position <-  min( seq(2:length(acf_values))[acf_values[-length(acf_values)] + acf_values[-1] < 0] )
      ## position contains NA for constants, needs to be addressed here:

      if( !is.na(position) ) {
        ## Pad with NA's so that all vectors are of same length. Saves me storing the position object
        ## pad with NA so that mean() can be calculated over differing rho's per chains
        rho <- append(acf_values[1:position+1], rep(NA, length(acf_values)-position), after = position)
      } else {
        rho <- rep(NA, n)
      }
    })

    s2m_rtm <- lapply(seq_along(chain_variances), function(i) {
      chain_variances[i] * chain_rho[,i]
    })

    ## average across chains
    ## Convert list to a matrix
    matrix_form <- do.call(cbind, s2m_rtm)
    avg_s2m_rtm <- rowMeans(matrix_form, na.rm = TRUE)
  })

  ## Eq 10: W - mn_s2m_ptm
  numerator <- matrix(W, nrow = nrow(mn_s2m_ptm), ncol = length(W), byrow = TRUE) - mn_s2m_ptm
  rho_t <- 1 - numerator / matrix(vtp, nrow = nrow(numerator), ncol = length(vtp), byrow = TRUE)

  n_eff <- round( n*m / ( 1+2*colSums(rho_t, na.rm = TRUE) ) )

}  else if(n_eff == 'stan') {
  ## Based on rstan, takes forever...
  #monitor_results <- rstan::monitor(samples_array, print = FALSE)
  #n_eff <- monitor_results$n_eff

  samples_array <- array(NA, dim = c(iterations, chains, parameters))
    for (i in seq_along(combined_chains)) {
      samples_array[, i, ] <- combined_chains[[i]]
    }
    monitor_results <- rstan::monitor(samples_array, print = FALSE)
    n_eff <- monitor_results$n_eff
}

## Extract and print R-hat values
out$rhat_values <- Rhat
if( any(out$rhat_values[!is.na(out$rhat_values)] > 1.1) ) warning("Some R-hat values are greater than 1.10 -- increase warmup and/or sampling iterations." )

## Effective sample size
out$n_eff <- n_eff

## Save the rest to the out object
out$samples <- combined_chains
out$nimble_constants <- constants
out$X_location_names <- colnames(data$X) # save fixed effects names for summary table renaming
out$X_scale <- data$X_scale
out$Z_location_names <- colnames(data$Z) # save random effects names for summary table renaming
out$Z_scale <- data$Z_scale
out$Y <- data.frame("group_id" = group_id, "Y" = data$Y)
out$workers <- workers
