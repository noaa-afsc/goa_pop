pop_mod <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  M <- exp(log_M)
  a50C <- exp(log_a50C)
  a50S <- exp(log_a50S)
  q <- exp(log_q)
  F50 <- exp(log_F50)
  F40 <- exp(log_F40)
  F35 <- exp(log_F35)

  spawn_fract <- (spawn_mo - 1) / 12
  spawn_adj <- exp(-M)^(spawn_fract)

  A <- nrow(age_error)
  A1 <- length(ages)
  T <- length(years)
  Ts <- sum(srv_ind)
  Tfa <- sum(fish_age_ind)
  Tsa <- sum(srv_age_ind)
  Tfs <- sum(fish_size_ind)
  L <- length(length_bins)
  g <- 0.00001

  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T)
  initNat <- rep(0, A)
  catch_pred <- rep(0, T)
  srv_pred <- rep(0, Ts)
  srv_var <- rep(0, Ts)
  fish_age_pred <- matrix(0, A1, Tfa)
  srv_age_pred <- matrix(0, A1, Tsa)
  fish_size_pred <- matrix(0, L, Tfs)
  spawn_bio <- tot_bio <- rep(0, T)
  N_spr <- sb_spr <- matrix(1, A, 4)

  # priors -----------------
  if (like_type == "admb") {
    nll_M <- (log(M) - log(mean_M))^2 / (2 * cv_M^2)
    nll_q <- (log(q) - log(mean_q))^2 / (2 * cv_q^2)
    nll_sigmaR <- (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2)
  } else {
    nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
    nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
    nll_sigmaR <- -RTMB::dnorm(
      log(sigmaR / mean_sigmaR),
      0,
      cv_sigmaR,
      log = TRUE
    )
  }

  # function alt ----
  ddirmult <- function(obs, pred, iss, ln_theta, log = TRUE) {
    y_obs <- iss * obs
    dirichlet_parm <- exp(ln_theta) * iss
    logres <- lgamma(iss + 1) - sum(lgamma(y_obs + 1))
    logres <- logres + lgamma(dirichlet_parm) - lgamma(iss + dirichlet_parm)
    logres <- logres +
      sum(lgamma(y_obs + dirichlet_parm * pred) - lgamma(dirichlet_parm * pred))
    if (log) return(logres) else return(exp(logres))
  }

  # selectivity ----
  to_one <- function(x) {
    x / max(x)
  }

  sel_logistic <- function(age, a50, delta, adj = 0) {
    x <- age + adj
    1 / (1 + exp(-log(19) * (x - a50) / delta))
  }

  sel_gamma <- function(age, b50, delta, adj = 0) {
    x <- age + adj
    denom <- 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
    ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
  }

  slx_block <- matrix(0, A, 4)
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0)
  if (block2 == 'gamma') {
    slx_block[, 2] <- RTMButils::sel_gamma(1:A, a50C[2], deltaC[2], adj = 0)
    slx_block[, 3] <- RTMButils::sel_gamma(1:A, a50C[3], deltaC[3], adj = 0)
    slx_block[, 4] <- RTMButils::sel_gamma(1:A, a50C[4], deltaC[4], adj = 0)
  } else {
    slx_block[, 3] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0)
    slx_block[, 2] <- (slx_block[, 1] + slx_block[, 3]) * 0.5
    slx_block[, 4] <- to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj = 0))
    slx_block[, 3] <- to_one(slx_block[, 3])
  }

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]]
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0)

  # mortality ----
  Ft <- exp(log_mean_F + log_Ft)
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t]
    Zat[, t] <- Fat[, t] + M
  }
  Sat <- exp(-Zat)
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  ## Nat ----
  for (t in 1:T) {
    if (bias_switch == 1) {
      bias_adj <- bias_ramp[t] * ((sigmaR^2) / 2)
      Nat[1, t] <- exp(log_mean_R - bias_adj + log_Rt[t])
    } else {
      Nat[1, t] <- exp(log_mean_R + log_Rt[t])
    }
  }
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1])
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M))

  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1]
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1]
  }

  # recruitment likelihood
  switch(
    like_type,
    "admb" = {
      like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
        (2 * sigmaR^2) +
        length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
        wt_rec_var
    },
    "rtmb" = {
      like_rec_main <- -sum(RTMB::dnorm(log_Rt, 0, sigmaR, log = TRUE))
      like_rec_init <- -sum(RTMB::dnorm(init_log_Rt, 0, sigmaR, log = TRUE))
      like_rec <- (like_rec_main + like_rec_init) * wt_rec_var
    }
  )

  # Calculate recruits and biomasses
  recruits <- Nat[1, ]
  spawn_bio <- colSums(Nat * wt_mature)
  tot_bio <- colSums(Nat * waa)

  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa)
  switch(
    like_type,
    "admb" = {
      ssqcatch <- sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2)
    },
    "rtmb" = {
      sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
      like_catch <- -sum(RTMB::dnorm(
        log(catch_obs + g),
        log(catch_pred + g),
        sigma_catch,
        log = TRUE
      )) *
        catch_wt
    }
  )

  ## survey biomass ----
  isrv <- 1
  srv_like <- 0.0
  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q

      switch(
        like_type,
        "admb" = {
          srv_like <- srv_like +
            (log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
              (2 * (srv_sd[isrv] / srv_obs[isrv])^2)
        },
        "rtmb" = {
          log_sd <- sqrt(log(1 + srv_cv[isrv]^2))

          # Toggle bias correction on or off
          if (do_bias_correct) {
            mu <- log(srv_pred[isrv]) - 0.5 * log_sd^2
          } else {
            mu <- log(srv_pred[isrv])
          }

          srv_like <- srv_like -
            RTMB::dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
        }
      )
      isrv <- isrv + 1
    }
  }
  like_srv <- srv_like * srv_wt

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1
  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Calculate base predictions (shared between both likelihood types)
      fish_age_pred[, icomp] <- colSums((Cat[, t] / sum(Cat[, t])) * age_error)

      switch(
        like_type,
        "admb" = {
          offset <- offset -
            fish_age_iss[icomp] *
              sum((fish_age_obs[, icomp] + g) * log(fish_age_obs[, icomp] + g))

          fish_age_lk <- fish_age_lk -
            sum(
              fish_age_iss[icomp] *
                (fish_age_obs[, icomp] + g) *
                log(fish_age_pred[, icomp] + g)
            )
        },
        "rtmb" = {
          # 1. Add robustness constant 'g' and re-normalize so probabilities sum to 1
          pred_prob <- fish_age_pred[, icomp] + g
          pred_prob <- pred_prob / sum(pred_prob)

          # 2. Convert observed proportions to "counts" using the input sample size
          obs_counts <- fish_age_obs[, icomp] * fish_age_iss[icomp]

          # 3. Calculate exact multinomial likelihood
          fish_age_lk <- fish_age_lk -
            RTMB::dmultinom(x = obs_counts, prob = pred_prob, log = TRUE)
        }
      )
      icomp <- icomp + 1
    }
  }
  # Offset shifts ADMB likelihood so a perfect fit = 0.0 (stays 0.0 for RTMB)
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1
  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Calculate base predictions (shared between both likelihood types)
      srv_age_pred[, icomp] <- colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      )

      switch(
        like_type,
        "admb" = {
          offset_sa <- offset_sa -
            srv_age_iss[icomp] *
              sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))

          srv_age_lk <- srv_age_lk -
            srv_age_iss[icomp] *
              sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
        },
        "rtmb" = {
          # 1. Add robustness constant 'g' and re-normalize
          pred_prob <- srv_age_pred[, icomp] + g
          pred_prob <- pred_prob / sum(pred_prob)

          # 2. Convert observed proportions to "counts" using effective sample size
          obs_counts <- srv_age_obs[, icomp] * srv_age_iss[icomp]

          # 3. Calculate exact multinomial likelihood
          srv_age_lk <- srv_age_lk -
            RTMB::dmultinom(x = obs_counts, prob = pred_prob, log = TRUE)
        }
      )
      icomp <- icomp + 1
    }
  }
  # Offset shifts ADMB likelihood so a perfect fit = 0.0 (stays 0.0 for RTMB)
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt

  ## fishery size comp ----
  fish_size_lk <- 0.0
  offset_fs <- 0.0
  icomp <- 1
  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Calculate base predictions (shared between both likelihood types)
      fish_size_pred[, icomp] <- colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      )

      switch(
        like_type,
        "admb" = {
          offset_fs <- offset_fs -
            fish_size_iss[icomp] *
              sum(
                (fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g)
              )

          fish_size_lk <- fish_size_lk -
            fish_size_iss[icomp] *
              sum(
                (fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g)
              )
        },
        "rtmb" = {
          # 1. Add robustness constant 'g' and re-normalize
          pred_prob <- fish_size_pred[, icomp] + g
          pred_prob <- pred_prob / sum(pred_prob)

          # 2. Convert observed proportions to "counts" using effective sample size
          obs_counts <- fish_size_obs[, icomp] * fish_size_iss[icomp]

          # 3. Calculate exact multinomial likelihood
          fish_size_lk <- fish_size_lk -
            RTMB::dmultinom(x = obs_counts, prob = pred_prob, log = TRUE)
        }
      )
      icomp <- icomp + 1
    }
  }
  # Offset shifts ADMB likelihood so a perfect fit = 0.0 (stays 0.0 for RTMB)
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt

  # SPR ------------------------

  valid_idx <- which(
    years >= (1977 + ages[1]) & years <= (max(years) - ages[1])
  )
  n_rec <- length(valid_idx)
  yrs_rec <- years[valid_idx]
  pred_rec <- mean(Nat[1, valid_idx])
  stdev_rec <- sqrt(
    sum((log_Rt[valid_idx] - mean(log_Rt[valid_idx]))^2) / (n_rec - 1)
  )

  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  SB0 <- sum(sb_spr[, 1])
  SBF50 <- sum(sb_spr[, 2])
  SBF40 <- sum(sb_spr[, 3])
  SBF35 <- sum(sb_spr[, 4])

  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # nll ----
  if (like_type == "admb") {
    nll <- ssqcatch
  } else {
    nll <- like_catch
  }
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  if (like_type == "admb") {
    RTMB::REPORT(ssqcatch)
  } else {
    RTMB::REPORT(like_catch)
  }
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)

  return(nll)
}
base <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- (log(M) - log(mean_M))^2 / (2 * cv_M^2) # Prior on natural mortality
  nll_q <- (log(q) - log(mean_q))^2 / (2 * cv_q^2) # Prior on survey catchability
  nll_sigmaR <- (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

  # selectivity ----

  to_one <- function(x) {
    x / max(x)
  }

  sel_logistic <- function(age, a50, delta, adj = 0) {
    x <- age + adj
    sel <- 1 / (1 + exp(-log(19) * (x - a50) / delta))
    # sel / max(sel)
    sel
  }

  sel_gamma <- function(age, b50, delta, adj = 0) {
    x <- age + adj
    denom <- 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
    sel <- ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
    # sel / max(sel)
    sel
  }
  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 3] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 3: double logistic selectivity
  slx_block[, 2] <- (slx_block[, 1] + slx_block[, 3]) * 0.5 # Block 2: average of blocks 1 and 3
  slx_block[, 4] <- to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj = 0)) # Block 4: double logistic selectivity
  slx_block[, 3] <- to_one(slx_block[, 3]) # Normalize block 3 selectivity - must be done after block 2 to match ADMB

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  ssqcatch <- sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)

  ## survey biomass ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      srv_like <- srv_like +
        sum(
          (log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
            (2 * (srv_sd[isrv] / srv_obs[isrv])^2)
        )
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * age_error
      ))
      # Offset for multinomial likelihood
      offset <- offset -
        fish_age_iss[icomp] *
          sum(
            (fish_age_obs[, icomp] + g) *
              log(fish_age_obs[, icomp] + g)
          )
      # Multinomial likelihood for age composition
      fish_age_lk <- fish_age_lk -
        sum(
          fish_age_iss[icomp] *
            (fish_age_obs[, icomp] + g) *
            log(fish_age_pred[, icomp] + g)
        )
      icomp <- icomp + 1
    }
  }
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[, icomp] <- as.numeric(colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      ))
      # Offset for multinomial likelihood
      offset_sa <- offset_sa -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk <- srv_age_lk -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

  ## fishery size comp ----
  icomp <- 1
  fish_size_lk <- 0.0
  offset_fs <- 0.0

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      ))
      # Offset for multinomial likelihood
      offset_fs <- offset_fs -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk <- fish_size_lk -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # likelihood/penalties --------------------
  like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
    (2 * sigmaR^2) +
    length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
    wt_rec_var
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  # nll ----
  nll <- ssqcatch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

base_gamma <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- (log(M) - log(mean_M))^2 / (2 * cv_M^2) # Prior on natural mortality
  nll_q <- (log(q) - log(mean_q))^2 / (2 * cv_q^2) # Prior on survey catchability
  nll_sigmaR <- (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

  # selectivity ----

  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 2] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 2: gamma selectivity
  slx_block[, 3] <- sel_gamma(1:A, a50C[3], deltaC[3], adj = 0) # Block 3: gamma
  slx_block[, 4] <- sel_gamma(1:A, a50C[4], deltaC[4], adj = 0) # Block 4: gamma selectivity

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  ssqcatch <- sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)

  ## survey biomass ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      srv_like <- srv_like +
        sum(
          (log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
            (2 * (srv_sd[isrv] / srv_obs[isrv])^2)
        )
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * age_error
      ))
      # Offset for multinomial likelihood
      offset <- offset -
        fish_age_iss[icomp] *
          sum(
            (fish_age_obs[, icomp] + g) *
              log(fish_age_obs[, icomp] + g)
          )
      # Multinomial likelihood for age composition
      fish_age_lk <- fish_age_lk -
        sum(
          fish_age_iss[icomp] *
            (fish_age_obs[, icomp] + g) *
            log(fish_age_pred[, icomp] + g)
        )
      icomp <- icomp + 1
    }
  }
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[, icomp] <- as.numeric(colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      ))
      # Offset for multinomial likelihood
      offset_sa <- offset_sa -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk <- srv_age_lk -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

  ## fishery size comp ----
  icomp <- 1
  fish_size_lk <- 0.0
  offset_fs <- 0.0

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      ))
      # Offset for multinomial likelihood
      offset_fs <- offset_fs -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk <- fish_size_lk -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # likelihood/penalties --------------------
  like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
    (2 * sigmaR^2) +
    length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
    wt_rec_var
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  # nll ----
  nll <- ssqcatch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

srv_like <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
  nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
  nll_sigmaR <- -RTMB::dnorm(
    log(sigmaR / mean_sigmaR),
    0,
    cv_sigmaR,
    log = TRUE
  )

  # selectivity ----
  to_one <- function(x) {
    x / max(x)
  }

  sel_logistic <- function(age, a50, delta, adj = 0) {
    x <- age + adj
    sel <- 1 / (1 + exp(-log(19) * (x - a50) / delta))
    # sel / max(sel)
    sel
  }

  sel_gamma <- function(age, b50, delta, adj = 0) {
    x <- age + adj
    denom <- 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
    sel <- ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
    # sel / max(sel)
    sel
  }
  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 3] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 3: double logistic selectivity
  slx_block[, 2] <- (slx_block[, 1] + slx_block[, 3]) * 0.5 # Block 2: average of blocks 1 and 3
  slx_block[, 4] <- to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj = 0)) # Block 4: double logistic selectivity
  slx_block[, 3] <- to_one(slx_block[, 3]) # Normalize block 3 selectivity - must be done after block 2 to match ADMB

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(RTMB::dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  )) *
    catch_wt
  ## survey biomass - bias corrected ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      log_sd <- sqrt(log(1 + srv_cv[isrv]^2))
      mu <- log(srv_pred[isrv]) - 0.5 * log_sd^2
      srv_like <- srv_like -
        RTMB::dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * age_error
      ))
      # Offset for multinomial likelihood
      offset <- offset -
        fish_age_iss[icomp] *
          sum(
            (fish_age_obs[, icomp] + g) *
              log(fish_age_obs[, icomp] + g)
          )
      # Multinomial likelihood for age composition
      fish_age_lk <- fish_age_lk -
        sum(
          fish_age_iss[icomp] *
            (fish_age_obs[, icomp] + g) *
            log(fish_age_pred[, icomp] + g)
        )
      icomp <- icomp + 1
    }
  }
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[, icomp] <- as.numeric(colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      ))
      # Offset for multinomial likelihood
      offset_sa <- offset_sa -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk <- srv_age_lk -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

  ## fishery size comp ----
  icomp <- 1
  fish_size_lk <- 0.0
  offset_fs <- 0.0

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      ))
      # Offset for multinomial likelihood
      offset_fs <- offset_fs -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk <- fish_size_lk -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # likelihood/penalties --------------------
  like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
    (2 * sigmaR^2) +
    length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
    wt_rec_var
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  # nll ----
  nll <- like_catch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

srv_like_gamma <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
  nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
  nll_sigmaR <- -RTMB::dnorm(
    log(sigmaR / mean_sigmaR),
    0,
    cv_sigmaR,
    log = TRUE
  )

  # selectivity ----

  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 2] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 3: double logistic selectivity
  slx_block[, 3] <- sel_gamma(1:A, a50C[3], deltaC[3], adj = 0) # Block 2: average of blocks 1 and 3
  slx_block[, 4] <- sel_gamma(1:A, a50C[4], deltaC[4], adj = 0) # Block 4: double logistic selectivity

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  )) *
    catch_wt
  ## survey biomass - bias corrected ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      # survey likelihood (lognormal)
      CV <- srv_sd[isrv] / srv_obs[isrv]
      log_sd <- sqrt(log(1 + CV^2))
      mu <- log(srv_pred[isrv]) - 0.5 * log_sd^2
      srv_like <- srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * age_error
      ))
      # Offset for multinomial likelihood
      offset <- offset -
        fish_age_iss[icomp] *
          sum(
            (fish_age_obs[, icomp] + g) *
              log(fish_age_obs[, icomp] + g)
          )
      # Multinomial likelihood for age composition
      fish_age_lk <- fish_age_lk -
        sum(
          fish_age_iss[icomp] *
            (fish_age_obs[, icomp] + g) *
            log(fish_age_pred[, icomp] + g)
        )
      icomp <- icomp + 1
    }
  }
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[, icomp] <- as.numeric(colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      ))
      # Offset for multinomial likelihood
      offset_sa <- offset_sa -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk <- srv_age_lk -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

  ## fishery size comp ----
  icomp <- 1
  fish_size_lk <- 0.0
  offset_fs <- 0.0

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      ))
      # Offset for multinomial likelihood
      offset_fs <- offset_fs -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk <- fish_size_lk -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # likelihood/penalties --------------------
  like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
    (2 * sigmaR^2) +
    length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
    wt_rec_var
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  # nll ----
  nll <- like_catch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

srv_like2 <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning
  wt_mature <- waa * maa * sex_ratio

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
  nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
  nll_sigmaR <- -RTMB::dnorm(
    log(sigmaR / mean_sigmaR),
    0,
    cv_sigmaR,
    log = TRUE
  )

  # function alt ----
  ddirmult <- function(obs, pred, iss, ln_theta, log = TRUE) {
    # expected counts and Dirichlet parameter
    y_obs <- iss * obs
    dirichlet_parm <- exp(ln_theta) * iss
    # base integration constants
    logres <- lgamma(iss + 1) - sum(lgamma(y_obs + 1))
    # theta scaling
    logres <- logres + lgamma(dirichlet_parm) - lgamma(iss + dirichlet_parm)
    logres <- logres +
      sum(lgamma(y_obs + dirichlet_parm * pred) - lgamma(dirichlet_parm * pred))

    if (log) {
      return(logres)
    } else {
      return(exp(logres))
    }
  }
  # selectivity ----
  to_one <- function(x) {
    x / max(x)
  }

  sel_logistic <- function(age, a50, delta, adj = 0) {
    x <- age + adj
    sel <- 1 / (1 + exp(-log(19) * (x - a50) / delta))
    # sel / max(sel)
    sel
  }

  sel_gamma <- function(age, b50, delta, adj = 0) {
    x <- age + adj
    denom <- 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
    sel <- ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
    # sel / max(sel)
    sel
  }
  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 3] <- sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 3: double logistic selectivity
  slx_block[, 2] <- (slx_block[, 1] + slx_block[, 3]) * 0.5 # Block 2: average of blocks 1 and 3
  slx_block[, 4] <- to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj = 0)) # Block 4: double logistic selectivity
  slx_block[, 3] <- to_one(slx_block[, 3]) # Normalize block 3 selectivity - must be done after block 2 to match ADMB

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    # year-specific bias adjustment if using random effects bias_switch = 0
    if (bias_switch == 1) {
      bias_adj <- bias_ramp[t] * ((sigmaR^2) / 2)
      Nat[1, t] <- exp(log_mean_R - bias_adj + log_Rt[t]) # recruitment in year t
    } else {
      Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # recruitment in year t
    }
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # recriuitment likelihood
  like_rec_main <- -sum(RTMB::dnorm(log_Rt, 0, sigmaR, log = TRUE))
  like_rec_init <- -sum(RTMB::dnorm(init_log_Rt, 0, sigmaR, log = TRUE))

  if (bias_switch == 1) {
    like_rec_main <- like_rec_main - sum((1 - 0.5 * bias_ramp) * log(sigmaR))
  }
  like_rec <- (like_rec_main + like_rec_init) * wt_rec_var

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(RTMB::dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  )) *
    catch_wt

  ## survey biomass ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      log_sd <- sqrt(log(1 + srv_cv[isrv]^2))
      mu <- log(srv_pred[isrv] + g) - 0.5 * log_sd^2
      srv_like <- srv_like -
        RTMB::dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  # fishery age comp ----
  pred <- (t(Cat) / colSums(Cat)) %*% age_error
  fish_age_pred <- t(pred)[, fish_age_ind == 1] + g
  fish_age_pred <- t(t(fish_age_pred) / colSums(fish_age_pred))
  fish_age_lk <- 0.0

  switch(
    comp_type,
    "mult" = {
      for (i in 1:ncol(fish_age_pred)) {
        obs_count <- fish_age_obs[, i] * fish_age_iss[i]

        fish_age_lk <- fish_age_lk -
          RTMB::dmultinom(x = obs_count, prob = fish_age_pred[, i], log = TRUE)
      }
    },
    "dm" = {
      for (i in 1:ncol(fish_age_pred)) {
        fish_age_lk <- fish_age_lk -
          ddirmult(
            obs = fish_age_obs[, i],
            pred = fish_age_pred[, i],
            iss = fish_age_iss[i],
            ln_theta = log_theta_fac,
            log = TRUE
          )
      }
    }
  )

  like_fish_age <- fish_age_lk * fish_age_wt

  ## survey age comp ----
  pred_srv <- (t(Nat * slx_srv) / colSums(Nat * slx_srv)) %*% age_error
  srv_age_pred <- t(pred_srv)[, srv_age_ind == 1] + g
  srv_age_pred <- t(t(srv_age_pred) / colSums(srv_age_pred))
  srv_age_lk <- 0.0
  switch(
    comp_type,
    "mult" = {
      for (i in 1:ncol(srv_age_pred)) {
        obs_count <- srv_age_obs[, i] * srv_age_iss[i]
        srv_age_lk <- srv_age_lk -
          RTMB::dmultinom(x = obs_count, prob = srv_age_pred[, i], log = TRUE)
      }
    },

    "dm" = {
      for (i in 1:ncol(srv_age_pred)) {
        srv_age_lk <- srv_age_lk -
          ddirmult(
            obs = srv_age_obs[, i],
            pred = srv_age_pred[, i],
            iss = srv_age_iss[i],
            ln_theta = log_theta_sac,
            log = TRUE
          )
      }
    }
  )

  like_srv_age <- srv_age_lk * srv_age_wt

  ## fishery size comp ----
  size_years <- which(fish_size_ind == 1)
  fish_size_lk <- 0.0
  for (i in 1:Tfs) {
    t <- size_years[i]
    cat_prop <- Cat[, t] / sum(Cat[, t])
    pred <- as.vector(t(cat_prop) %*% saa_array[,, fish_saa_ind[t]]) + g
    fish_size_pred[, i] <- pred / sum(pred)
  }

  switch(
    comp_type,
    "mult" = {
      for (i in 1:Tfs) {
        obs_count <- fish_size_obs[, i] * fish_size_iss[i]
        fish_size_lk <- fish_size_lk -
          RTMB::dmultinom(x = obs_count, prob = fish_size_pred[, i], log = TRUE)
      }
    },

    "dm" = {
      for (i in 1:Tfs) {
        fish_size_lk <- fish_size_lk -
          ddirmult(
            obs = fish_size_obs[, i],
            pred = fish_size_pred[, i],
            iss = fish_size_iss[i],
            ln_theta = log_theta_fsc,
            log = TRUE
          )
      }
    }
  )
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # nll ----
  nll <- like_catch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  if (comp_type == "dm") {
    RTMB::REPORT(log_theta_fac)
    RTMB::REPORT(log_theta_sac)
    RTMB::REPORT(log_theta_fsc)
  }

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

srv_like2_gamma <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning
  wt_mature <- waa * maa * sex_ratio

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
  nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
  nll_sigmaR <- -RTMB::dnorm(
    log(sigmaR / mean_sigmaR),
    0,
    cv_sigmaR,
    log = TRUE
  )

  # function alt ----
  ddirmult <- function(obs, pred, iss, ln_theta, log = TRUE) {
    # expected counts and Dirichlet parameter
    y_obs <- iss * obs
    dirichlet_parm <- exp(ln_theta) * iss
    # base integration constants
    logres <- lgamma(iss + 1) - sum(lgamma(y_obs + 1))
    # theta scaling
    logres <- logres + lgamma(dirichlet_parm) - lgamma(iss + dirichlet_parm)
    logres <- logres +
      sum(lgamma(y_obs + dirichlet_parm * pred) - lgamma(dirichlet_parm * pred))

    if (log) {
      return(logres)
    } else {
      return(exp(logres))
    }
  }
  # selectivity ----

  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- RTMButils::sel_logistic(1:A, a50C[1], deltaC[1], adj = 0) # Block 1: logistic selectivity
  slx_block[, 2] <- RTMButils::sel_gamma(1:A, a50C[2], deltaC[2], adj = 0) # Block 3: double logistic selectivity
  slx_block[, 3] <- RTMButils::sel_gamma(1:A, a50C[3], deltaC[3], adj = 0) # Block 2: average of blocks 1 and 3
  slx_block[, 4] <- RTMButils::sel_gamma(1:A, a50C[4], deltaC[4], adj = 0) # Block 4: double logistic selectivity

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- RTMButils::sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    # year-specific bias adjustment if using random effects bias_switch = 0
    if (bias_switch == 1) {
      bias_adj <- bias_ramp[t] * ((sigmaR^2) / 2)
      Nat[1, t] <- exp(log_mean_R - bias_adj + log_Rt[t]) # recruitment in year t
    } else {
      Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # recruitment in year t
    }
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # recriuitment likelihood
  like_rec_main <- -sum(RTMB::dnorm(log_Rt, 0, sigmaR, log = TRUE))
  like_rec_init <- -sum(RTMB::dnorm(init_log_Rt, 0, sigmaR, log = TRUE))

  if (bias_switch == 1) {
    like_rec_main <- like_rec_main - sum((1 - 0.5 * bias_ramp) * log(sigmaR))
  }
  like_rec <- (like_rec_main + like_rec_init) * wt_rec_var

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  )) *
    catch_wt

  ## survey biomass ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      log_sd <- sqrt(log(1 + srv_cv[isrv]^2))
      mu <- log(srv_pred[isrv] + g) - 0.5 * log_sd^2
      srv_like <- srv_like -
        RTMB::dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  # fishery age comp ----
  pred <- (t(Cat) / colSums(Cat)) %*% age_error
  fish_age_pred <- t(pred)[, fish_age_ind == 1] + g
  fish_age_pred <- t(t(fish_age_pred) / colSums(fish_age_pred))
  fish_age_lk <- 0.0

  switch(
    comp_type,
    "mult" = {
      for (i in 1:ncol(fish_age_pred)) {
        obs_count <- fish_age_obs[, i] * fish_age_iss[i]
        fish_age_lk <- fish_age_lk -
          RTMB::dmultinom(x = obs_count, prob = fish_age_pred[, i], log = TRUE)
      }
    },
    "dm" = {
      for (i in 1:ncol(fish_age_pred)) {
        fish_age_lk <- fish_age_lk -
          ddirmult(
            obs = fish_age_obs[, i],
            pred = fish_age_pred[, i],
            iss = fish_age_iss[i],
            ln_theta = log_theta_fac,
            log = TRUE
          )
      }
    }
  )

  like_fish_age <- fish_age_lk * fish_age_wt

  ## survey age comp ----
  pred_srv <- (t(Nat * slx_srv) / colSums(Nat * slx_srv)) %*% age_error
  srv_age_pred <- t(pred_srv)[, srv_age_ind == 1] + g
  srv_age_pred <- t(t(srv_age_pred) / colSums(srv_age_pred))
  srv_age_lk <- 0.0
  switch(
    comp_type,
    "mult" = {
      for (i in 1:ncol(srv_age_pred)) {
        obs_count <- srv_age_obs[, i] * srv_age_iss[i]
        srv_age_lk <- srv_age_lk -
          RTMB::dmultinom(x = obs_count, prob = srv_age_pred[, i], log = TRUE)
      }
    },

    "dm" = {
      for (i in 1:ncol(srv_age_pred)) {
        srv_age_lk <- srv_age_lk -
          ddirmult(
            obs = srv_age_obs[, i],
            pred = srv_age_pred[, i],
            iss = srv_age_iss[i],
            ln_theta = log_theta_sac,
            log = TRUE
          )
      }
    }
  )

  like_srv_age <- srv_age_lk * srv_age_wt

  ## fishery size comp ----
  size_years <- which(fish_size_ind == 1)
  fish_size_lk <- 0.0
  for (i in 1:Tfs) {
    t <- size_years[i]
    cat_prop <- Cat[, t] / sum(Cat[, t])
    pred <- as.vector(t(cat_prop) %*% saa_array[,, fish_saa_ind[t]]) + g
    fish_size_pred[, i] <- pred / sum(pred)
  }

  switch(
    comp_type,
    "mult" = {
      for (i in 1:Tfs) {
        obs_count <- fish_size_obs[, i] * fish_size_iss[i]
        fish_size_lk <- fish_size_lk -
          RTMB::dmultinom(x = obs_count, prob = fish_size_pred[, i], log = TRUE)
      }
    },

    "dm" = {
      for (i in 1:Tfs) {
        fish_size_lk <- fish_size_lk -
          ddirmult(
            obs = fish_size_obs[, i],
            pred = fish_size_pred[, i],
            iss = fish_size_iss[i],
            ln_theta = log_theta_fsc,
            log = TRUE
          )
      }
    }
  )
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # nll ----
  nll <- like_catch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  if (comp_type == "dm") {
    RTMB::REPORT(log_theta_fac)
    RTMB::REPORT(log_theta_sac)
    RTMB::REPORT(log_theta_fsc)
  }

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}
# ADMB model 20.1 bridged to RTMB
# no input parameter or data changes
# base <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")

#   RTMB::getAll(pars, data)

#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass

#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning

#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero

#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit

#   # priors -----------------
#   # Priors on key parameters (negative log-likelihood contributions)
#   nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

# # selectivity ----

# 	to_one <- function(x) {
#   	x / max(x)
# 	}

# 	sel_logistic <- function(age, a50, delta, adj=0) {
#   	x = age + adj
#   	sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
#   	# sel / max(sel)
#   	sel
# 	}

#   sel_gamma <- function(age, b50, delta, adj=0) {
#   	x = age + adj
#   	denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
#   	sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
#   	# sel / max(sel)
#   	sel
# 	}
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,3] = sel_gamma(1:A, a50C[2], deltaC[2], adj=0) # Block 3: double logistic selectivity
#   slx_block[,2] = (slx_block[,1] + slx_block[,3]) * 0.5          # Block 2: average of blocks 1 and 3
#   slx_block[,4] = to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
#   slx_block[,3] = to_one(slx_block[,3])                          # Normalize block 3 selectivity - must be done after block 2 to match ADMB

#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }

#   slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)

#   # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year

#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)

#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }

#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year

#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)

#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)

#   ## survey biomass ----
#   isrv = 1
#   srv_like = 0.0

#   for(t in 1:T) {
#     if(srv_ind[t]==1) {
#       srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#       # Survey likelihood (lognormal, using observed and predicted survey biomass)
#       srv_like = srv_like + sum((log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
#                                   (2 * (srv_sd[isrv] / srv_obs[isrv])^2))
#       isrv = isrv + 1
#     }
#   }

#   like_srv = srv_like * srv_wt # Weighted survey likelihood

#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] *
#         sum((fish_age_obs[,icomp] + g) *
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] *
#                                         (fish_age_obs[,icomp] + g) *
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0

#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD

#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F35 * slx_fish[A,T])))

#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }

#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35

#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2

#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec

#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)

#   # nll ----
#   nll = ssqcatch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen

#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)

#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)

#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)

#   RTMB::REPORT(ssqcatch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }

# # RTMB model with:
# # survey likelihood bias correction
# # no input parameter or data changes
# srv_like <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")

#   RTMB::getAll(pars, data)

#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass

#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning

#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero

#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit

#   # priors -----------------
#   # # Priors on key parameters (negative log-likelihood contributions)
#   # nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   # nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   # nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

#   # priors ----
#    nll_M = -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
#    nll_q= -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
#    nll_sigmaR= -RTMB::dnorm(log(sigmaR / mean_sigmaR), 0, cv_sigmaR, log = TRUE)

#   to_one <- function(x) {
#   	x / max(x)
# 	}

# 	sel_logistic <- function(age, a50, delta, adj=0) {
#   	x = age + adj
#   	sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
#   	# sel / max(sel)
#   	sel
# 	}

#   sel_gamma <- function(age, b50, delta, adj=0) {
#   	x = age + adj
#   	denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
#   	sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
#   	# sel / max(sel)
#   	sel
# 	}
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,3] = sel_gamma(1:A, a50C[2], deltaC[2], adj=0) # Block 3: double logistic selectivity
#   slx_block[,2] = (slx_block[,1] + slx_block[,3]) * 0.5          # Block 2: average of blocks 1 and 3
#   slx_block[,4] = to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
#   slx_block[,3] = to_one(slx_block[,3])                          # Normalize block 3 selectivity - must be done after block 2 to match ADMB

#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }

#   slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)

#   # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year

#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)

#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }

#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year

#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)

#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
#   sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
#   like_catch <- -sum(dnorm(log(catch_obs + g), log(catch_pred + g), sigma_catch, log = TRUE)) * catch_wt
#   ## survey biomass - bias corrected ----
#   isrv = 1
#   srv_like = 0.0

#   for(t in 1:T) {
#     if(srv_ind[t]==1) {
#       srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#       # survey likelihood (lognormal)
#       CV = srv_sd[isrv] / srv_obs[isrv]
#       log_sd = sqrt(log(1 + CV^2))
#       mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
#       srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
#       isrv = isrv + 1
#     }
#   }

#   like_srv = srv_like * srv_wt # Weighted survey likelihood

#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] *
#         sum((fish_age_obs[,icomp] + g) *
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] *
#                                         (fish_age_obs[,icomp] + g) *
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0

#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD

#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F35 * slx_fish[A,T])))

#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }

#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35

#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2

#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec

#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)

#   # nll ----
#   nll = like_catch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen

#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)

#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)

#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)

#   RTMB::REPORT(like_catch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }

# # RTMB model with:
# # survey likelihood bias correction
# # add gamma parameters, no data changes
# srv_like_gamma <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")

#   RTMB::getAll(pars, data)

#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass

#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning

#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero

#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit

#   # priors -----------------
#   # Priors on key parameters (negative log-likelihood contributions)
#   nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

#   to_one <- function(x) {
#   	x / max(x)
# 	}

# 	sel_logistic <- function(age, a50, delta, adj=0) {
#   	x = age + adj
#   	sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
#   	# sel / max(sel)
#   	sel
# 	}

#   sel_gamma <- function(age, b50, delta, adj=0) {
#   	x = age + adj
#   	denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
#   	sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
#   	# sel / max(sel)
#   	sel
# 	}
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,2] = to_one(sel_gamma(1:A, a50C[2], deltaC[2], adj=0))
#   slx_block[,3] = to_one(sel_gamma(1:A, a50C[3], deltaC[3], adj=0))
#   slx_block[,4] = to_one(sel_gamma(1:A, a50C[4], deltaC[4], adj=0))

#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }

#     slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)

#     # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year

#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)

#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }

#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year

#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)

#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
#   sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
#   like_catch <- -sum(dnorm(log(catch_obs + g), log(catch_pred + g), sigma_catch, log = TRUE)) * catch_wt

#   ## survey biomass - bias corrected ----
#   srv_like = 0.0
#   srv_pred_all = colSums(Nat * slx_srv * waa) * q
#   srv_pred = srv_pred_all[srv_ind == 1]
#   log_sd = sqrt(log(1 + srv_cv^2))
#   mu = log(srv_pred + g) - 0.5 * log_sd^2
#   srv_like = -sum(RTMB::dnorm(log(srv_obs), mu, log_sd, log = TRUE))
#   like_srv = srv_like * srv_wt

#   # isrv = 1
#   # srv_like = 0.0

#   # for(t in 1:T) {
#   #   if(srv_ind[t]==1) {
#   #     srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#   #     # survey likelihood (lognormal)
#   #     CV = srv_sd[isrv] / srv_obs[isrv]
#   #     log_sd = sqrt(log(1 + CV^2))
#   #     mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
#   #     srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
#   #     isrv = isrv + 1
#   #   }
#   # }

#   # like_srv = srv_like * srv_wt # Weighted survey likelihood

#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] *
#         sum((fish_age_obs[,icomp] + g) *
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] *
#                                         (fish_age_obs[,icomp] + g) *
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0

#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD

#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F35 * slx_fish[A,T])))

#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }

#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35

#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2

#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec

#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)

#   # nll ----
#   nll = like_catch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen

#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)

#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)

#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)

#   RTMB::REPORT(like_catch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }

# # RTMB model with:
# # survey likelihood bias correction
# # no input parameter or data changes
# srv_like_double_log <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")

#   RTMB::getAll(pars, data)

#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass

#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning

#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero

#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit

#   # priors -----------------
#   # Priors on key parameters (negative log-likelihood contributions)
#   nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

#   # selectivity ----

# 	sel_logistic <- function(age, a50, delta, adj=0) {
#   	x = age + adj
#   	sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
#   	sel / max(sel)
# 	}

#   sel_gamma <- function(age, b50, delta, adj=0) {
#   	x = age + adj
#   	denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
#   	sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
#   	sel / max(sel)

# 	}
#   sel_double_logistic <- function(age, a50_a, a50_d, delta, delta2, adj=0) {
#     a50_d = a50_a + a50_d
#     x = age + adj
#     asc = 1 / (1 + exp(-delta * (x - a50_a)))
#     desc = 1 / (1 + exp(delta2 * (x - a50_d)))
#     sel = asc * desc
#     sel / max(sel)
#   }
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,2] = sel_double_logistic(1:A, a50C[4], a50C[5], deltaC[4], deltaC[5], adj=0) # Block 2: gamma selectivity
#   slx_block[,3] = sel_gamma(1:A, a50C[2], deltaC[2], adj=0) # Block 3: gamma selectivity
#   slx_block[,4] = sel_gamma(1:A, a50C[3], deltaC[3], adj=0) # Block 4: gamma selectivity

#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }

#   slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)

#   # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year

#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)

#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }

#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year

#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)

#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
#   sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
#   like_catch <- -sum(dnorm(log(catch_obs + g), log(catch_pred + g), sigma_catch, log = TRUE)) * catch_wt
#   ## survey biomass - bias corrected ----
#   isrv = 1
#   srv_like = 0.0

#   for(t in 1:T) {
#     if(srv_ind[t]==1) {
#       srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#       # survey likelihood (lognormal)
#       CV = srv_sd[isrv] / srv_obs[isrv]
#       log_sd = sqrt(log(1 + CV^2))
#       mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
#       srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
#       isrv = isrv + 1
#     }
#   }

#   like_srv = srv_like * srv_wt # Weighted survey likelihood

#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] *
#         sum((fish_age_obs[,icomp] + g) *
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] *
#                                         (fish_age_obs[,icomp] + g) *
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0

#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD

#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F35 * slx_fish[A,T])))

#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }

#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35

#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2

#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec

#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)

#   # nll ----
#   nll = like_catch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen

#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)

#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)

#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)

#   RTMB::REPORT(like_catch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }

# # RTMB model with:
# # survey likelihood bias correction
# # all time blocks logistic selectivity
# all_logistic <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")

#   RTMB::getAll(pars, data)

#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass

#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning

#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero

#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit

#   # priors -----------------
#   # Priors on key parameters (negative log-likelihood contributions)
#   nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

#   # selectivity ----

# 	sel_logistic <- function(age, a50, delta, adj=0) {
#   	x = age + adj
#   	sel = 1 / (1 + exp(-log(19) * (x - a50) / delta))
#   	sel / max(sel)
# 	}

#   sel_gamma <- function(age, b50, delta, adj=0) {
#   	x = age + adj
#   	denom = 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
#   	sel = ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
#   	sel / max(sel)

# 	}
#   sel_double_logistic <- function(age, a50_a, a50_d, delta, delta2, adj=0) {
#     a50_d = a50_a + a50_d
#     x = age + adj
#     asc = 1 / (1 + exp(-delta * (x - a50_a)))
#     desc = 1 / (1 + exp(delta2 * (x - a50_d)))
#     sel = asc * desc
#     sel / max(sel)
#   }
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,2] = sel_logistic(1:A, a50C[2], deltaC[2], adj=0) # Block 2: gamma selectivity
#   slx_block[,3] = sel_logistic(1:A, a50C[3], deltaC[3], adj=0) # Block 3: gamma selectivity
#   slx_block[,4] = sel_logistic(1:A, a50C[4], deltaC[4], adj=0) # Block 4: gamma selectivity

#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }

#   slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)

#   # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year

#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)

#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }

#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year

#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)

#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
#   sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
#   like_catch <- -sum(dnorm(log(catch_obs + g), log(catch_pred + g), sigma_catch, log = TRUE)) * catch_wt
#   ## survey biomass - bias corrected ----
#   isrv = 1
#   srv_like = 0.0

#   for(t in 1:T) {
#     if(srv_ind[t]==1) {
#       srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#       # survey likelihood (lognormal)
#       CV = srv_sd[isrv] / srv_obs[isrv]
#       log_sd = sqrt(log(1 + CV^2))
#       mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
#       srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
#       isrv = isrv + 1
#     }
#   }

#   like_srv = srv_like * srv_wt # Weighted survey likelihood

#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] *
#         sum((fish_age_obs[,icomp] + g) *
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] *
#                                         (fish_age_obs[,icomp] + g) *
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1

#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0

#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD

#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) /
#     (1 - exp(-(M + F35 * slx_fish[A,T])))

#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }

#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35

#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2

#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec

#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)

#   # nll ----
#   nll = like_catch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen

#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)

#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)

#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)

#   RTMB::REPORT(like_catch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }

all_double_log <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M <- exp(log_M) # Natural mortality
  a50C <- exp(log_a50C) # Age at 50% selectivity (fishery)
  a50S <- exp(log_a50S) # Age at 50% selectivity (survey)
  q <- exp(log_q) # Survey catchability
  F50 <- exp(log_F50) # Fishing mortality at 50% spawning biomass
  F40 <- exp(log_F40) # Fishing mortality at 40% spawning biomass
  F35 <- exp(log_F35) # Fishing mortality at 35% spawning biomass

  # Spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # Fraction of year before spawning
  spawn_adj <- exp(-M)^(spawn_fract) # Mortality adjustment for spawning

  # Index values and dimensions
  A <- nrow(age_error) # Number of ages in model
  A1 <- length(ages) # Number of ages in comps
  T <- sum(catch_ind) # Number of fishery years
  Ts <- sum(srv_ind) # Number of survey years
  Tfa <- sum(fish_age_ind) # Number of fishery age comp years
  Tsa <- sum(srv_age_ind) # Number of survey age comp years
  Tfs <- sum(fish_size_ind) # Number of fishery size comp years
  L <- length(length_bins) # Number of length bins
  g <- 0.00001 # Small number to avoid division by zero

  # Containers for model outputs
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat <- rep(0, A) # Initial numbers-at-age
  catch_pred <- rep(0, T) # Predicted catch
  srv_pred <- rep(0, Ts) # Predicted survey index
  srv_var <- rep(0, Ts) # Survey variance
  fish_age_pred <- matrix(0, A1, Tfa) # Predicted fishery age comps
  srv_age_pred <- matrix(0, A1, Tsa) # Predicted survey age comps
  fish_size_pred <- matrix(0, L, Tfs) # Predicted fishery size comps
  spawn_bio <- tot_bio <- rep(0, T) # Spawning and total biomass
  N_spr <- sb_spr <- matrix(1, A, 4) # Numbers and spawning biomass per recruit

  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M <- (log(M) - log(mean_M))^2 / (2 * cv_M^2) # Prior on natural mortality
  nll_q <- (log(q) - log(mean_q))^2 / (2 * cv_q^2) # Prior on survey catchability
  nll_sigmaR <- (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability

  # selectivity ----

  sel_logistic <- function(age, a50, delta, adj = 0) {
    x <- age + adj
    sel <- 1 / (1 + exp(-log(19) * (x - a50) / delta))
    sel / max(sel)
  }

  sel_gamma <- function(age, b50, delta, adj = 0) {
    x <- age + adj
    denom <- 0.5 * (sqrt(b50^2 + 4 * delta^2) - b50)
    sel <- ((x / b50)^(b50 / denom)) * exp((b50 - x) / denom)
    sel / max(sel)
  }
  sel_double_logistic <- function(age, a50_a, a50_d, delta, delta2, adj = 0) {
    a50_d <- a50_a + a50_d
    x <- age + adj
    asc <- 1 / (1 + exp(-delta * (x - a50_a)))
    desc <- 1 / (1 + exp(delta2 * (x - a50_d)))
    sel <- asc * desc
    sel / max(sel)
  }
  slx_block <- matrix(0, A, 4) # Selectivity blocks for fishery
  slx_block[, 1] <- sel_double_logistic(
    1:A,
    a50C[1],
    a50C[2],
    deltaC[1],
    deltaC[2],
    adj = 0
  ) # Block 1: logistic selectivity
  slx_block[, 2] <- sel_double_logistic(
    1:A,
    a50C[3],
    a50C[4],
    deltaC[3],
    deltaC[4],
    adj = 0
  ) # Block 2: gamma selectivity
  slx_block[, 3] <- sel_double_logistic(
    1:A,
    a50C[5],
    a50C[6],
    deltaC[5],
    deltaC[6],
    adj = 0
  ) # Block 3: gamma selectivity
  slx_block[, 4] <- sel_double_logistic(
    1:A,
    a50C[7],
    a50C[8],
    deltaC[7],
    deltaC[8],
    adj = 0
  ) # Block 4: gamma selectivity

  for (t in 1:T) {
    slx_fish[, t] <- slx_block[, fish_block_ind[t]] # Assign selectivity by year
  }

  slx_srv <- sel_logistic(1:A, a50S, deltaS, adj = 0) # Survey selectivity (logistic)

  # mortality ----
  # Calculate fishing mortality for each year
  Ft <- exp(log_mean_F + log_Ft) # Annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # Fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # Total mortality at age and year
  }
  Sat <- exp(-Zat) # Survivorship at age and year

  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for (t in 1:T) {
    Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- exp(log_mean_R - (a - 1) * M + init_log_Rt[a - 1]) # Initial numbers for ages 2 to A-1
  }
  Nat[A, 1] <- exp(log_mean_R - (A - 1) * M) / (1 - exp(-M)) # Plus group (oldest age class)

  # Remaining columns: survivors from previous year
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # Survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # Plus group accumulates survivors
  }

  # Calculate recruits and biomasses
  recruits <- Nat[1, ] # Recruitment time series
  spawn_bio <- colSums(Nat * wt_mature) # Spawning biomass by year
  tot_bio <- colSums(Nat * waa) # Total biomass by year

  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj <- Sat[, T]^(spawn_fract)
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  ## catch ----
  # Calculate predicted catch at age and year
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa) # Predicted catch biomass
  # ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  )) *
    catch_wt
  ## survey biomass - bias corrected ----
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q # Predicted survey index
      # survey likelihood (lognormal)
      CV <- srv_sd[isrv] / srv_obs[isrv]
      log_sd <- sqrt(log(1 + CV^2))
      mu <- log(srv_pred[isrv]) - 0.5 * log_sd^2
      srv_like <- srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }

  like_srv <- srv_like * srv_wt # Weighted survey likelihood

  ## fishery age comp ----
  fish_age_lk <- 0.0
  offset <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * age_error
      ))
      # Offset for multinomial likelihood
      offset <- offset -
        fish_age_iss[icomp] *
          sum(
            (fish_age_obs[, icomp] + g) *
              log(fish_age_obs[, icomp] + g)
          )
      # Multinomial likelihood for age composition
      fish_age_lk <- fish_age_lk -
        sum(
          fish_age_iss[icomp] *
            (fish_age_obs[, icomp] + g) *
            log(fish_age_pred[, icomp] + g)
        )
      icomp <- icomp + 1
    }
  }
  fish_age_lk <- fish_age_lk - offset
  like_fish_age <- fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood

  ## survey age comp ----
  srv_age_lk <- 0.0
  offset_sa <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[, icomp] <- as.numeric(colSums(
        (Nat[, t] * slx_srv) / sum(Nat[, t] * slx_srv) * age_error
      ))
      # Offset for multinomial likelihood
      offset_sa <- offset_sa -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_obs[, icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk <- srv_age_lk -
        srv_age_iss[icomp] *
          sum((srv_age_obs[, icomp] + g) * log(srv_age_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  srv_age_lk <- srv_age_lk - offset_sa
  like_srv_age <- srv_age_lk * srv_age_wt # Weighted survey age comp likelihood

  ## fishery size comp ----
  icomp <- 1
  fish_size_lk <- 0.0
  offset_fs <- 0.0

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[, icomp] <- as.numeric(colSums(
        (Cat[, t] / sum(Cat[, t])) * saa_array[,, fish_saa_ind[t]]
      ))
      # Offset for multinomial likelihood
      offset_fs <- offset_fs -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_obs[, icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk <- fish_size_lk -
        fish_size_iss[icomp] *
          sum((fish_size_obs[, icomp] + g) * log(fish_size_pred[, icomp] + g))
      icomp <- icomp + 1
    }
  }
  fish_size_lk <- fish_size_lk - offset_fs
  like_fish_size <- fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood

  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt, pred_rec = Nat[1, ], year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df <- df[years >= (1977 + ages[1]) & years <= (max(years) - ages[1]), ]
  n_rec <- nrow(df)
  yrs_rec <- df$year
  pred_rec <- mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec <- sqrt(
    sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)
  ) # Recruitment SD

  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for (a in 2:A) {
    N_spr[a, 1] <- N_spr[a - 1, 1] * exp(-M)
    N_spr[a, 2] <- N_spr[a - 1, 2] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_spr[a, 3] <- N_spr[a - 1, 3] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_spr[a, 4] <- N_spr[a - 1, 4] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A, 1] <- N_spr[A - 1, 1] * exp(-M) / (1 - exp(-M))
  N_spr[A, 2] <- N_spr[A - 1, 2] *
    exp(-(M + F50 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_spr[A, 3] <- N_spr[A - 1, 3] *
    exp(-(M + F40 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_spr[A, 4] <- N_spr[A - 1, 4] *
    exp(-(M + F35 * slx_fish[A - 1, T])) /
    (1 - exp(-(M + F35 * slx_fish[A, T])))

  # Calculate spawning biomass per recruit for reference points
  for (a in 1:A) {
    sb_spr[a, 1] <- N_spr[a, 1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a, 2] <- N_spr[a, 2] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F50 * slx_fish[a, T]))
    sb_spr[a, 3] <- N_spr[a, 3] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F40 * slx_fish[a, T]))
    sb_spr[a, 4] <- N_spr[a, 4] *
      wt_mature[a] *
      exp(-spawn_fract * (M + F35 * slx_fish[a, T]))
  }

  # Calculate reference point spawning biomasses
  SB0 <- sum(sb_spr[, 1]) # Unfished spawning biomass per recruit
  SBF50 <- sum(sb_spr[, 2]) # Spawning biomass per recruit at F50
  SBF40 <- sum(sb_spr[, 3]) # Spawning biomass per recruit at F40
  SBF35 <- sum(sb_spr[, 4]) # Spawning biomass per recruit at F35

  # SPR penalties to enforce reference point constraints
  sprpen <- 100. * (SBF50 / SB0 - 0.5)^2
  sprpen <- sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen <- sprpen + 100. * (SBF35 / SB0 - 0.35)^2

  # Scale reference points by mean recruitment
  B0 <- SB0 * pred_rec
  B40 <- SBF40 * pred_rec
  B35 <- SBF35 * pred_rec

  # likelihood/penalties --------------------
  like_rec <- (sum(c(log_Rt, init_log_Rt)^2) /
    (2 * sigmaR^2) +
    length(c(log_Rt, init_log_Rt)) * log(sigmaR)) *
    wt_rec_var
  f_regularity <- wt_fmort_reg * sum(log_Ft^2)

  # nll ----
  nll <- like_catch
  nll <- nll + like_srv
  nll <- nll + like_fish_age
  nll <- nll + like_srv_age
  nll <- nll + like_fish_size
  nll <- nll + like_rec
  nll <- nll + f_regularity
  nll <- nll + nll_M
  nll <- nll + nll_q
  nll <- nll + nll_sigmaR
  nll <- nll + sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}
urm <- function(pars, data) {
  # load ----
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  RTMB::getAll(pars, data)

  # setup ----
  # transform log parameters to natural scale
  M <- exp(log_M) # natural mortality
  q <- exp(log_q) # survey catchability
  slx_pars <- exp(log_slx_pars) # selectivity
  F50 <- exp(log_F50)
  F40 <- exp(log_F40)
  F35 <- exp(log_F35)
  R_init <- exp(log_mean_R_init)
  F_init <- exp(log_F_init)

  # spawning adjustments
  spawn_fract <- (spawn_mo - 1) / 12 # fraction of year before spawning

  # dimensions ----
  A <- nrow(age_error) # ages in model
  A1 <- length(ages) # ages in comps
  T <- length(years) # fishery years
  Ts <- sum(srv_ind) # survey years
  Tfa <- sum(fish_age_ind) # fishery age comp years
  Tsa <- sum(srv_age_ind) # survey age comp years
  Tfs <- sum(fish_size_ind) # fishery size comp years
  L <- length(length_bins) # lengths in comps
  n_curves <- nrow(slx_pars) # number of selectivity curves
  g <- 0.00001 # small number to avoid division by zero
  Z_init <- cum_Z_init <- rep(0, A) # total mortality pre-model equilibrium

  # containers ----
  Bat <- Cat <- Nat <- Fat <- Zat <- Sat <- slx_fish <- matrix(0, A, T) # biomass, catch, numbers, F, Z, S, selectivity at age
  catch_pred <- spawn_bio <- tot_bio <- rep(0, T)
  srv_pred <- srv_var <- rep(0, Ts)
  fish_age_pred <- matrix(0, A1, Tfa)
  srv_age_pred <- matrix(0, A1, Tsa)
  fish_size_pred <- matrix(0, L, Tfs)
  slx_curves <- matrix(0, A, n_curves)

  # priors ----
  nll_M <- -RTMB::dnorm(log(M), log(mean_M), cv_M, log = TRUE)
  nll_q <- -RTMB::dnorm(log(q), log(mean_q), cv_q, log = TRUE)
  nll_sigmaR <- -RTMB::dnorm(
    log(sigmaR / mean_sigmaR),
    0,
    cv_sigmaR,
    log = TRUE
  )

  wt_mature <- waa * maa * sex_ratio

  # selectivity ----
  for (s in 1:n_curves) {
    slx_curves[, s] <- RTMButils::get_slx(
      ages = 1:A,
      type = slx_type[s],
      pars = slx_pars[s, ],
      adj = 1
    )
  }

  # assign to the fishery over time
  for (t in 1:T) {
    slx_fish[, t] <- slx_curves[, fish_block_ind[t]]
  }

  # survey selectivity
  slx_srv <- slx_curves[, srv_slx_ind]

  # mortality ----
  Ft <- exp(log_mean_F + log_Ft) # annual fishing mortality on natural scale
  for (t in 1:T) {
    Fat[, t] <- Ft[t] * slx_fish[, t] # fishing mortality at age and year
    Zat[, t] <- Fat[, t] + M # total mortality at age and year
  }
  Sat <- exp(-Zat) # survivorship at age and year

  # numbers-at-age ----

  for (t in 1:T) {
    # year-specific bias adjustment if using random effects bias_switch = 0
    if (bias_switch == 1) {
      bias_adj <- bias_ramp[t] * ((sigmaR^2) / 2)
      Nat[1, t] <- exp(log_mean_R - bias_adj + log_Rt[t]) # recruitment in year t
    } else {
      Nat[1, t] <- exp(log_mean_R + log_Rt[t]) # recruitment in year t
    }
  }

  # initial total mortality at age (historical pre-model equilibrium)
  for (a in 1:A) {
    Z_init[a] <- M + F_init * slx_fish[a, 1]
  }

  # cumulative mortality experienced by each cohort prior to Year 1
  for (a in 2:A) {
    cum_Z_init[a] <- cum_Z_init[a - 1] + Z_init[a - 1]
  }

  # column 1: initial numbers-at-age for each cohort cells 1 to A-1
  for (a in 2:(A - 1)) {
    Nat[a, 1] <- R_init * exp(-cum_Z_init[a] + init_log_Rt[a - 1])
  }

  # specific deviation for initial plus group - deterministic
  Nat[A, 1] <- R_init * exp(-cum_Z_init[A]) / (1 - exp(-Z_init[A]))

  # forward
  for (t in 2:T) {
    for (a in 2:A) {
      Nat[a, t] <- Nat[a - 1, t - 1] * Sat[a - 1, t - 1] # survivors from previous age and year
    }
    Nat[A, t] <- Nat[A, t] + Nat[A, t - 1] * Sat[A, t - 1] # plus group
  }

  # recriuitment likelihood
  like_rec_main <- -sum(RTMB::dnorm(log_Rt, 0, sigmaR, log = TRUE))
  if (bias_switch == 1) {
    like_rec_main <- like_rec_main - sum((1 - 0.5 * bias_ramp) * log(sigmaR))
  }
  like_rec_init <- -sum(RTMB::dnorm(init_log_Rt, 0, sigmaR_init, log = TRUE))
  like_rec <- (like_rec_main + like_rec_init) * wt_rec_var

  recruits <- Nat[1, ]
  spawn_bio <- colSums(Nat * wt_mature)
  tot_bio <- colSums(Nat * waa)

  # mortality adjustment for last year spawning adjustment
  spawn_adj <- Sat[, T]^spawn_fract
  spawn_bio[T] <- sum(Nat[, T] * spawn_adj * wt_mature)

  # catch ----
  Cat <- Fat / Zat * Nat * (1 - Sat)
  catch_pred <- colSums(Cat * waa)
  sigma_catch <- sqrt(log(catch_cv^2 + 1.0))
  like_catch <- -sum(RTMB::dnorm(
    log(catch_obs + g),
    log(catch_pred + g),
    sigma_catch,
    log = TRUE
  ))

  # survey biomass ----
  # w/log-normal bias correction
  isrv <- 1
  srv_like <- 0.0

  for (t in 1:T) {
    if (srv_ind[t] == 1) {
      srv_pred[isrv] <- sum(Nat[, t] * slx_srv * waa) * q
      log_sd <- sqrt(log(1 + srv_cv[isrv]^2))
      mu <- log(srv_pred[isrv] + g) - 0.5 * log_sd^2
      srv_like <- srv_like -
        RTMB::dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv <- isrv + 1
    }
  }
  like_srv <- srv_like * srv_wt

  # fishery age comp ----
  fish_age_lk <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_age_ind[t] == 1) {
      # predicted age composition (with ageing error)
      fish_age_pred[, icomp] <- as.vector(colSums(
        (Cat[, t] / (sum(Cat[, t]))) * age_error
      ))
      obs_count <- fish_age_obs[, icomp] * fish_age_iss[icomp]
      fish_age_lk <- fish_age_lk -
        RTMB::dmultinom(
          x = obs_count,
          prob = fish_age_pred[, icomp],
          log = TRUE
        )
      icomp <- icomp + 1
    }
  }
  like_fish_age <- fish_age_lk * fish_age_wt # weighted fishery age comp likelihood

  # survey age comp ----
  srv_age_lk <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (srv_age_ind[t] == 1) {
      # predicted age composition (with ageing error)
      srv_age_pred[, icomp] <- as.vector(colSums(
        (Nat[, t] * slx_srv) / (sum(Nat[, t] * slx_srv)) * age_error
      ))
      obs_count <- srv_age_obs[, icomp] * srv_age_iss[icomp]
      srv_age_lk <- srv_age_lk -
        RTMB::dmultinom(x = obs_count, prob = srv_age_pred[, icomp], log = TRUE)
      icomp <- icomp + 1
    }
  }
  like_srv_age <- srv_age_lk * srv_age_wt # weighted survey age comp likelihood

  # fishery size comp ----
  fish_size_lk <- 0.0
  icomp <- 1

  for (t in 1:T) {
    if (fish_size_ind[t] == 1) {
      fish_size_pred[, icomp] <- as.vector(colSums(
        (Cat[, t] / (sum(Cat[, t]))) * saa_array[,, fish_saa_ind[t]]
      ))
      obs_count <- fish_size_obs[, icomp] * fish_size_iss[icomp]
      fish_size_lk <- fish_size_lk -
        RTMB::dmultinom(
          x = obs_count,
          prob = fish_size_pred[, icomp],
          log = TRUE
        )
      icomp <- icomp + 1
    }
  }
  like_fish_size <- fish_size_lk * fish_size_wt # weighted fishery size comp likelihood

  # SPR ----
  Nspr_0 <- rep(1, A)
  for (a in 2:A) {
    Nspr_0[a] <- Nspr_0[a - 1] * exp(-M)
  }
  Nspr_0[A] <- Nspr_0[A] / (1 - exp(-M))
  SB0 <- sum(Nspr_0 * wt_mature * exp(-spawn_fract * M)) # unfished spawning biomass per recruit

  # biomass reference points
  N_50 <- N_40 <- N_35 <- rep(1, A)
  for (a in 2:A) {
    N_50[a] <- N_50[a - 1] * exp(-(M + F50 * slx_fish[a - 1, T]))
    N_40[a] <- N_40[a - 1] * exp(-(M + F40 * slx_fish[a - 1, T]))
    N_35[a] <- N_35[a - 1] * exp(-(M + F35 * slx_fish[a - 1, T]))
  }
  N_50[A] <- N_50[A] / (1 - exp(-(M + F50 * slx_fish[A, T])))
  N_40[A] <- N_40[A] / (1 - exp(-(M + F40 * slx_fish[A, T])))
  N_35[A] <- N_35[A] / (1 - exp(-(M + F35 * slx_fish[A, T])))

  sbpr_50 <- sum(
    N_50 * wt_mature * exp(-spawn_fract * (M + F50 * slx_fish[, T]))
  )
  sbpr_40 <- sum(
    N_40 * wt_mature * exp(-spawn_fract * (M + F40 * slx_fish[, T]))
  )
  sbpr_35 <- sum(
    N_35 * wt_mature * exp(-spawn_fract * (M + F35 * slx_fish[, T]))
  )

  sprpen <- 1000.0 *
    ((sbpr_50 / SB0 - 0.50)^2 +
      (sbpr_40 / SB0 - 0.40)^2 +
      (sbpr_35 / SB0 - 0.35)^2)

  yrs_rec <- years[years >= (1977 + ages[1]) & years <= (max(years) - ages[1])]
  n_rec <- length(yrs_rec)
  pred_rec <- mean(Nat[1, years %in% yrs_rec])
  stdev_rec <- sqrt(
    sum((log_Rt[years %in% yrs_rec] - mean(log_Rt[years %in% yrs_rec]))^2) /
      (n_rec - 1)
  )

  B0 <- SB0 * pred_rec
  B50 <- B0 * 0.50
  B40 <- B0 * 0.40
  B35 <- B0 * 0.35

  # F penalty ----
  f_regularity <- -sum(RTMB::dnorm(
    x = log_Ft,
    mean = 0,
    sd = sigmaF,
    log = TRUE
  )) *
    wt_fmort_reg

  # joint negative log-likelihood
  nll <- like_catch +
    like_srv +
    like_fish_age +
    like_srv_age +
    like_fish_size +
    like_rec +
    f_regularity +
    nll_M +
    nll_q +
    nll_sigmaR +
    sprpen

  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(slx_pars)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(sigmaF)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(init_log_Rt)
  RTMB::ADREPORT(init_log_Rt)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_curves)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)

  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)

  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)

  # note: ADREPORT reference points to get uncertainty estimates (delta method)
  RTMB::ADREPORT(B0)
  RTMB::ADREPORT(B50)
  RTMB::ADREPORT(B40)
  RTMB::ADREPORT(B35)
  RTMB::ADREPORT(F35)
  RTMB::ADREPORT(F40)
  RTMB::ADREPORT(F50)

  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)

  RTMB::REPORT(like_catch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec_init)
  RTMB::REPORT(like_rec_main)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll)

  return(nll)
}
