#' Create a complete and dynamic JAGS model string for the SS-MELSM
#'
#' This function encapsulates all the logic for programmatically generating
#' a JAGS model for a Spike-and-Slab Mixed-Effects Location Scale Model.
#' It handles any number of random location and scale effects.
#'
#' @param Kr Integer. The number of random effects in the location model.
#' @param Sr Integer. The number of random effects in the scale model.
#'
#' @return A character string containing the complete JAGS model code.

 create_ss_melsm_jags_model <- function(Kr, Sr) {

    P <- Kr + Sr
    if (P == 0) { stop("Model must have at least one random effect.") }

    # --- Part 1: Define the corrected model template ---
    # Global definitions are now correctly placed after the group loop.
    jags_model_template <- "
    model {
      # --- Likelihood Section ---
      for(i in 1:N) {
        Y[i] ~ dnorm(mu[i], prec[i])
        prec[i] <- 1 / pow(tau[i], 2)
        ##DYNAMIC_LIKELIHOOD_LOCATION##
        ##DYNAMIC_LIKELIHOOD_SCALE##

        # Calculate and store log-likelihood for WAIC
        loglik[i] <- logdensity.norm(Y[i], mu[i], prec[i])
      }

      # --- Group-Level Section (j = group) ---
      for(j in 1:J) {
        ##DYNAMIC_GROUP_PRIORS##
      }

      # --- Model-Level (Global) Priors & Definitions ---
      ##DYNAMIC_GLOBAL_PRIORS##
    }"

    # --- Part 2: Generate the code snippets ---
    loc_lik_code <- if (Kr > 0) sprintf("mu[i] <- inprod(beta[], X[i, 1:K]) + inprod(v_final[groupid[i], 1:%d], Z[i, 1:%d])", Kr, Kr) else "mu[i] <- inprod(beta[], X[i, 1:K])"
    scale_lik_code <- if (Sr > 0) sprintf("log(tau[i]) <- inprod(zeta[], X_scale[i, 1:S]) + inprod(v_final[groupid[i], %d:%d], Z_scale[i, 1:%d])", Kr + 1, P, Sr) else "log(tau[i]) <- inprod(zeta[], X_scale[i, 1:S])"

    # Snippets for GLOBAL priors
    global_priors <- c()
    global_priors <- c(global_priors, "# Priors for fixed effects and variance components")
    global_priors <- c(global_priors, "for (k in 1:K) { beta[k] ~ dnorm(0, 1.0E-6) }")
    global_priors <- c(global_priors, "for (s in 1:S) { zeta[s] ~ dnorm(0, 1/3^2) }")
    global_priors <- c(global_priors, "for (p in 1:P) { sigma_rand[p] ~ dt(0, 1, 1) T(0,) }")

    if (P > 1) {
      n_cors <- P * (P - 1) / 2
      global_priors <- c(global_priors, "\n\t# Priors for correlations (vectorized)")
      global_priors <- c(global_priors, sprintf("\tfor (p in 1:%d) {", n_cors),
                         "\t\tzscore[p] ~ dnorm(0, 1)", "\t\trho_vec[p] <- tanh(zscore[p])", "\t}")
      global_priors <- c(global_priors, "\n\t# Map correlation vector to matrix using direct indexing")
      global_priors <- c(global_priors, sprintf("\tfor (j in 2:%d) {", P),
                         "\t\tfor (i in 1:(j-1)) {",
                         "\t\t\tR[i,j] <- rho_vec[ (j-1)*(j-2)/2 + i ]",
                         "\t\t\tR[j,i] <- R[i,j]", "\t\t}", "\t}")
    }
    global_priors <- c(global_priors, "\tfor(p in 1:P) { R[p,p] <- 1 }")

    global_priors <- c(global_priors, "\n\t# Cholesky decomposition of the correlation matrix")
    global_priors <- c(global_priors, "\tL[1,1] <- 1")
    if (P > 1) {
      for (i in 2:P) { global_priors <- c(global_priors, sprintf("\tL[%d,1] <- R[%d,1]", i, i)) }
      for (j in 2:P) {
        sum_sq_str <- paste0("pow(L[", j, ", ", 1:(j-1), "], 2)", collapse = " + ")
        global_priors <- c(global_priors, sprintf("\tL[%d,%d] <- sqrt(1 - (%s))", j, j, sum_sq_str))
        if (j < P) {
          for (i in (j + 1):P) {
            sum_prod_str <- paste0("L[", i, ", ", 1:(j-1), "] * L[", j, ", ", 1:(j-1), "]", collapse = " + ")
            global_priors <- c(global_priors, sprintf("\tL[%d,%d] <- (1/L[%d,%d]) * (R[%d,%d] - (%s))", i, j, j, j, i, j, sum_prod_str))
          }
        }
      }
      for (i in 1:(P - 1)) {
        for (j in (i + 1):P) {
          global_priors <- c(global_priors, sprintf("\tL[%d,%d] <- 0", i, j))
        }
      }
    }

    # Snippets for GROUP-LEVEL priors
    group_priors <- c()
    group_priors <- c(group_priors, "# Standard normal deviates")
    group_priors <- c(group_priors, "for(p in 1:P){ z[j,p] ~ dnorm(0,1) }")
    group_priors <- c(group_priors, "\n\t\t# Generate correlated random effects v = tau * L * z")
    for (k in 1:P) {
      terms <- paste0("L[", k, ", ", 1:k, "] * z[j, ", 1:k, "]", collapse = " + ")
      group_priors <- c(group_priors, sprintf("\t\tv[j,%d] <- sigma_rand[%d] * (%s)", k, k, terms))
    }
    if (Kr > 0) {
      group_priors <- c(group_priors, "\n\t\t# Assign location random effects (always on)",
        sprintf("\t\tfor (k in 1:%d) {", Kr), "\t\t\tv_final[j, k] <- v[j, k]", "\t\t}")
    }
    if (Sr > 0) {
      group_priors <- c(group_priors, "\n\t\t# Assign scale random effects (with spike-and-slab via loop)",
        sprintf("\t\tfor (k in 1:%d) {", Sr), "\t\t\tss[j, k] ~ dbern(bval[k, 1])",
        sprintf("\t\t\tv_final[j, %d + k] <- v[j, %d + k] * ss[j, k]", Kr, Kr), "\t\t}")
    }

    # --- Part 3: Assemble the final model string ---
    final_model <- jags_model_template
    final_model <- sub("##DYNAMIC_LIKELIHOOD_LOCATION##", loc_lik_code, final_model)
    final_model <- sub("##DYNAMIC_LIKELIHOOD_SCALE##", scale_lik_code, final_model)
    final_model <- sub("##DYNAMIC_GLOBAL_PRIORS##", paste(global_priors, collapse = "\n\t"), final_model)
    final_model <- sub("##DYNAMIC_GROUP_PRIORS##", paste(group_priors, collapse = "\n\t\t"), final_model)

    return(final_model)
  }
