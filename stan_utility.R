# Check transitions that ended with a divergence
check_div <- function(fit, quiet=FALSE) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)

  if (!quiet) print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                    n, N, 100 * n / N))
  if (n > 0) {
    if (!quiet) print('  Try running with larger adapt_delta to remove the divergences')
    if (quiet) return(FALSE)
  } else {
    if (quiet) return(TRUE)
  }
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10, quiet=FALSE) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)

  if (!quiet)
    print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                            n, N, max_depth, 100 * n / N))
  
  if (n > 0) {
    if (!quiet) print('  Run again with max_depth set to a larger value to avoid saturation')
    if (quiet) return(FALSE)
  } else {
    if (quiet) return(TRUE)
  }
}

# Checks the energy fraction of missing information (E-FMI)
check_energy <- function(fit, quiet=FALSE) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      if (!quiet) print(sprintf('Chain %s: E-FMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('E-FMI indicated no pathological behavior')
    if (quiet) return(TRUE)
  } else {
    if (!quiet) print('  E-FMI below 0.2 indicates you may need to reparameterize your model')
    if (quiet) return(FALSE)
  }
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit, quiet=FALSE) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]

  iter <- dim(extract(fit)[[1]])[[1]]

  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      if (!quiet) print(sprintf('n_eff / iter for parameter %s is %s!',
                        rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('n_eff / iter looks reasonable for all parameters')
    if (quiet) return(TRUE)
  }
  else {
    if (!quiet) print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
    if (quiet) return(FALSE)
  }
}

# Checks the potential scale reduction factors
check_rhat <- function(fit, quiet=FALSE) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]

  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      if (!quiet) print(sprintf('Rhat for parameter %s is %s!',
                        rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('Rhat looks reasonable for all parameters')
    if (quiet) return(TRUE)
  } else {
    if (!quiet) print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
    if (quiet) return(FALSE)
  }
}

check_all_diagnostics <- function(fit, quiet=FALSE) {
  if (!quiet) {
    check_n_eff(fit)
    check_rhat(fit)
    check_div(fit)
    check_treedepth(fit)
    check_energy(fit)
  } else {
    warning_code <- 0
    
    if (!check_n_eff(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 0))
    if (!check_rhat(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 1))
    if (!check_div(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 2))
    if (!check_treedepth(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 3))
    if (!check_energy(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 4))

    return(warning_code)
  }
}

parse_warning_code <- function(warning_code) {
  if (bitwAnd(warning_code, bitwShiftL(1, 0)))
    print("n_eff / iteration warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 1)))
    print("rhat warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 2)))
    print("divergence warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 3)))
    print("treedepth warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 4)))
    print("energy warning")
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))

  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent

  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]

  return(list(div_params, nondiv_params))
}

### BEGIN sbc.R
library(rstan)
library(foreach)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() + 2)

ensure_model = function(model) { 
  if (typeof(model) == "character")
    model = stan_model(model)
  model
}

get_stan_draws <- function(gen_model, gen_data, num_draws = 4000) {
  gendata <- sampling(ensure_model(gen_model),
                      data=gen_data,
                      seed=23489432,
                      warmup = 0,
                      chains = 1,
                      iter = num_draws,
                      control=list(max_treedepth = 12),
                      algorithm='Fixed_param')
  as.data.frame(gendata)
}

add_draws_to_constants <- function(draws, data_vars, data_rows, constants) {
  getVarNames <-function(name) {
    names(draws)[grep(paste("^", name, sep=""), names(draws))]
  }
  
  for (j in seq_along(data_vars)) {
    dv <- data_vars[j]
    r <- data_rows[j]
    if (r > 1) {
      constants[[dv]] <- matrix(unlist(draws[getVarNames(dv)], use.names=F), nrow=r)
    } else {
      constants[[dv]] <- unlist(draws[getVarNames(dv)], use.names=F)
    }
  }
  constants
}

thin.df = function(df, skip, desired_quantiles) {
  sapply(colnames(df), function(param) {
    thinned_param = df[seq(1, nrow(df), round(skip[param])),param]
    thinned_param[1:desired_quantiles]
  })
}

thin.fit = function(fit, target_data, desired_quantiles) {
  neffs = summary(fit)$summary[,"n_eff"]
  df = as.data.frame(fit)
  skip = ceiling(nrow(df) / neffs) # assume we will have the required number of neffs after
  needed = max(skip * desired_quantiles * 2) # that 2 is because they're half warmup
  if (needed > nrow(df)) {
    write(paste0("Redoing! needed ", needed, "\n"), file="redos.txt", append=T) 
    df = get_posterior_thetas(fit@stanmodel, target_data, thin=F, num_posterior_draws=needed)
  }
  thin.df(df, skip, desired_quantiles)
}

get_posterior_thetas <- function(target_model, target_data, thin, num_posterior_draws,
                                 desired_quantiles) {
  library(rstan)
  target_model = ensure_model(target_model)
  target_fit = sampling(target_model, data=target_data,
                        chains=1, iter=num_posterior_draws*2, seed=12345
                        , control=list(max_treedepth=20)
                        #, control=list(adapt_delta=0.98)
                        #, control=list(adapt_delta=0.99, max_treedepth=14)
  );
  diag = check_all_diagnostics(target_fit)
  
  if (thin) {
    thin.fit(target_fit, target_data, desired_quantiles)
  } else {
    as.data.frame(target_fit)
  }
}

get_quantiles <- function(param_names, aggregate, comparison) {
  function(theta0s, posterior_thetas) {
    sapply(param_names, function(param) {
      theta0 = theta0s[, param];
      samples_theta0 = posterior_thetas[, param];
      aggregate(comparison(samples_theta0, theta0))
    })
  }
}

get_param_names <- function(gdf, data_vars) {
  getVarNames <-function(name) {
    return(names(gdf)[grep(paste("^", name, sep=""), names(gdf))])
  }
  mapCat <- function(f, coll) {unlist(lapply(coll, f))}
  
  data_names <- mapCat(getVarNames, data_vars)
  nonParamNames <- union(data_names, c("lp__"))
  setdiff(names(gdf), nonParamNames)
}

sbc.ranks <- function(gen_model_file, gen_data, data_vars, data_rows, 
                             target_model_file, num_replicates=1000, num_posterior_draws=4000,
                             thin=F, desired_quantiles=50) {
  gdf <- get_stan_draws(gen_model_file, gen_data, num_draws=num_replicates)
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1], outfile="sbc.log")
  registerDoParallel(cl)
  tryCatch({
    if (num_replicates > nrow(gdf)) {
      stop(paste0("Only ", nrow(gdf), " draws for ", num_replicates, " replicates!"))
    }
    
    target_model = ensure_model(target_model_file)
    paramNames <- get_param_names(gdf, data_vars)
    get_quants <- get_quantiles(paramNames, sum, function(samples, p0) {samples > p0})
    check_all_diagnostics = check_all_diagnostics
    #replications <- sapply(1:num_replicates, function(i) {               
    replications <-  foreach(i=1:num_replicates, #.combine=rbind,
                             .export = c("get_posterior_thetas", "ensure_model", "thin.fit",
                                         "add_draws_to_constants", "thin.df"
                             )) %dopar% {
                               library(rstan)
                               target_data <- add_draws_to_constants(gdf[i,], data_vars, data_rows, gen_data)
                               results <- get_posterior_thetas(target_model, target_data, thin=thin,
                                                               num_posterior_draws=num_posterior_draws,
                                                               desired_quantiles = desired_quantiles)
                               get_quants(gdf[i,paramNames], results)
                             }
    #)
    
    replications = do.call(rbind, replications)
  }, finally= { stopCluster(cl) })
  list(theta0s=gdf[,paramNames], ranks=replications)
}
