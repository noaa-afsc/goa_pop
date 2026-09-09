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
