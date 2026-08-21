# 2026 GOA Pacific ocean perch assessment
# 2026-08
# ben.williams@noaa.gov

# structural changes
# m2020.1-2023 (aka m23.0) bridged to RTMB (m25.0)
# m25.0 - update with data through 2026
# m25.1 - change likelihood to use full form for catch and survey biomass
# m25.2 - m25.1 w/2nd fishery selectivity time block to gamma from average of 1st and 3rd

# input changes
# a - use surveyISS for survey_age_iss
# b - use model-based survey biomass
# c - use a & b
# m25.XXr - Francis reweight

# load ----
library(RTMButils)
library(afscassess)
library(afscdata)
library(patchwork)
library(tidyverse)
theme_set(afscassess::theme_report())

source(here::here(2026, 'r', "models.R"))

# globals ----
year <- 2026
rec_age <- 2
plus_age <- 25
lengths <- 16:45

# data ----
# 2023 data & pars from ADMB model
data0 <- readRDS(here::here(2025, 'rtmb_bridge', "data.rds"))
pars0 <- readRDS(here::here(2025, 'rtmb_bridge', "pars.rds"))

# update data to 2026
catch <- vroom::vroom(here::here(year, "data", "output", "fish_catch.csv"))
yld_rat <- vroom::vroom(here::here(year, "data", "output", "yld_rat.csv"))
srv <- vroom::vroom(here::here(
  year,
  "data",
  "raw",
  "goa_region_bts_biomass_data.csv"
)) %>%
  rename(biomass = biomass_mt) %>%
  mutate(
    cv = sqrt(biomass_var) / biomass,
    lse = sqrt(log(1 + cv^2)),
    lci = exp(-1.96 * lse) * biomass,
    uci = exp(1.96 * lse) * biomass,
    type = "db",
    id = "2025",
    sd = sqrt(biomass_var)
  ) %>%
  select(year, biomass, lci, uci, type, cv, sd, id)
fac <- vroom::vroom(here::here(year, "data", "output", "fish_age_comp.csv"))
sac <- vroom::vroom(here::here(year, "data", "output", "goa_bts_age_comp.csv"))
# fsc = vroom::vroom(here::here(year, "data", "output", "fish_length_comp.csv")) # this doesn't get updated
saa <- vroom::vroom(here::here(year, "data", "output", "saa.csv"))
waa <- vroom::vroom(here::here(year, "data", "output", "waa.csv"))$wbar
maa <- data0$maa
ae <- vroom::vroom(here::here(year, "data", "output", "ae_model.csv"))
size_age <- unname(saa[, -1])
saa_array <- data0$saa_array
saa_array[,, 2] <- as.matrix(size_age)


data <- list(
  ages = rec_age:plus_age,
  years = catch$year,
  length_bins = lengths,
  waa = waa,
  maa = maa,
  wt_mature = waa * maa * 0.5,
  spawn_mo = 5,
  fish_block_ind = c(data0$fish_block_ind, 4, 4, 4),
  slx_type = c("1", "1"),
  srv_slx_ind = 2,
  catch_obs = catch$catch,
  catch_ind = rep(1, nrow(catch)),
  catch_wt = rep(50, nrow(catch)),
  # catch_wt = 1,
  srv_yrs = srv$year,
  srv_ind = ifelse(catch$year %in% srv$year, 1, 0),
  srv_obs = srv$biomass,
  srv_sd = srv$sd,
  srv_cv = srv$cv,
  srv_wt = 1,
  fish_age_yrs = fac$year,
  fish_age_ind = ifelse(catch$year %in% fac$year, 1, 0),
  fish_age_obs = unname(t(as.matrix(fac[, -(1:4)]))),
  fish_age_iss = sqrt(fac$n_s),
  fish_age_wt = 1,
  srv_age_yrs = sac$year,
  srv_age_ind = ifelse(catch$year %in% sac$year, 1, 0),
  srv_age_iss = c(data0$srv_age_iss, 30.82, 26.6),
  srv_age_obs = unname(t(as.matrix(sac[, -1]))),
  srv_age_wt = 1,
  fish_size_yrs = data0$fish_size_yrs,
  fish_size_ind = ifelse(catch$year %in% data0$fish_size_yrs, 1, 0),
  fish_size_obs = unname(data0$fish_size_obs),
  fish_size_iss = data0$fish_size_iss,
  fish_size_wt = 1,
  age_error = unname(as.matrix(ae)),
  fish_saa_ind = as.integer(c(data0$fish_saa_ind, 2, 2, 2)), # 2 fishery saa indices
  saa_array = saa_array,
  wt_fmort_reg = 0.1,
  wt_rec_var = 1,
  mean_M = 0.0614,
  cv_M = 0.1,
  mean_q = 1.15,
  cv_q = 0.447213595,
  mean_sigmaR = 1.7,
  cv_sigmaR = 0.2,
  yield_ratio = yld_rat$yld,
  sex_ratio = 0.5,
  like_type = "admb",
  block2 = "avg",
  do_bias_correct = FALSE,
  bias_switch = 0 # 1 = use, 0 = don't
)

# original ADMB starting pars
pars <- list(
  log_M = log(0.0614),
  log_a50C = log(c(6, 10, 10)),
  deltaC = c(1.5, 4.5, 4.5),
  log_a50S = log(7.3),
  deltaS = 3.8,
  log_q = log(1.15),
  log_mean_R = 3.0,
  init_log_Rt = rep(0.0, 26),
  log_Rt = rep(0.0, sum(data$catch_ind)),
  log_mean_F = -2.5,
  log_Ft = rep(0.0, sum(data$catch_ind)),
  log_F35 = -2.1,
  log_F40 = -2.3,
  log_F50 = -2.7,
  sigmaR = 0.7
)

# limits ----
lower <- c(
  log_M = log(0.03), # Minimum natural mortality
  log_a50C = rep(log(0.5), 3), # Minimum age at 50% selectivity
  deltaC = rep(0.1, 3), # Minimum selectivity slope
  log_a50S = log(0.5), # Minimum survey selectivity age
  deltaS = 0.1, # Minimum survey selectivity slope
  log_q = log(0.01), # Minimum catchability
  log_mean_R = 2, # Minimum mean recruitment
  init_log_Rt = rep(-10, 26), # Minimum initial recruitment devs
  log_Rt = rep(-10, sum(data$catch_ind)), # Minimum recruitment devs
  log_mean_F = -5, # Minimum mean fishing mortality
  log_Ft = rep(-10, sum(data$catch_ind)), # Minimum F devs
  log_F35 = -5, # Minimum F35
  log_F40 = -5, # Minimum F40
  log_F50 = -5, # Minimum F50
  sigmaR = 0.1 # Minimum recruitment SD
)
upper <- c(
  log_M = log(0.27), # Maximum natural mortality
  log_a50C = rep(log(20), 3), # Maximum age at 50% selectivity
  deltaC = rep(20, 3), # Maximum selectivity slope
  log_a50S = log(20), # Maximum survey selectivity age
  deltaS = 20, # Maximum survey selectivity slope
  log_q = log(10), # Maximum catchability
  log_mean_R = 10, # Maximum mean recruitment
  init_log_Rt = rep(10, 26), # Maximum initial recruitment devs
  log_Rt = rep(10, sum(data$catch_ind)), # Maximum recruitment devs
  log_mean_F = 1, # Maximum mean fishing mortality
  log_Ft = rep(5, sum(data$catch_ind)), # Maximum F devs
  log_F35 = 1, # Maximum F35
  log_F40 = 1, # Maximum F40
  log_F50 = 1, # Maximum F50
  sigmaR = 1.9 # Maximum recruitment SD
)

# survey data ----
# 2023 vast
v23 <- read.csv(here::here(2023, "dev", "mb_vs_db", "vast_2023.csv")) %>%
  filter(Estimate > 0) %>%
  mutate(
    est = Estimate / 1000,
    ef = exp(1.96 * Std..Error.for.ln.Estimate.)
  ) %>%
  mutate(lwr = est / ef, upr = est * ef) %>%
  select(year = Time, est, lwr, upr) %>%
  mutate(id = 'mb-2023', type = "model-based")
# 2025 sdmTMB
v25 <- readRDS(here::here(2025, "data", "user_input", "vast_2025.rds")) %>%
  mutate(
    id = 'mb-2025',
    est = est / 1000,
    lwr = lwr / 1000,
    upr = upr / 1000,
    type = "model-based"
  )
# 2023 design-based
b23 <- read.csv(here::here(
  2023,
  'data',
  'output',
  'goa_total_bts_biomass.csv'
)) %>%
  mutate(id = 'db-2023', type = "design-based") %>%
  select(year, est = biomass, lwr = lci, upr = uci, id, type)

# data setup ----
# m.1 = change the likelihoods to RTMB
data1 <- data
data1$catch_wt <- 1
data1$catch_cv <- rep(0.10, length(data1$years))
data1$like_type <- "rtmb"

# m.2 = change to gamma slx 2nd time block
data2 <- data1
data2$block2 <- "gamma"

parsg <- pars
parsg$log_a50C <- c(pars$log_a50C, 0.693)
parsg$deltaC <- c(pars$deltaC, 4.5)

# set bounds
id <- max(grep("^log_a50C", names(lower)))
lowerg <- append(lower, c(log_a50C4 = log(0.1)), after = id)
upperg <- append(upper, c(log_a50C4 = log(20)), after = id)
id <- max(grep("^deltaC3", names(lower)))
lowerg <- append(lowerg, c(deltaC4 = 0.1), after = id)
upperg <- append(upperg, c(deltaC4 = 20), after = id)


# base model ----
# updates to run in pop_mod
data0$like_type <- "admb"
data0$block2 <- "avg"
data0$bias_switch <- 0

# run model, rewight for fun
m23 <- run_model(pop_mod, data0, pars0)
m23r <- run_model_reweight(m23)$model

fit_check(m23)
fit_check(m23r)

m23$proj
m23r$proj

saveRDS(m23, here::here(year, "alt", "results", "m23.rds"))
saveRDS(m23r, here::here(year, "alt", "results", "m23r.rds"))

# update to 2026 data
m25.0 <- run_model(
  pop_mod,
  data = data,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.0)
m25.1 <- run_model(
  pop_mod,
  data = data1,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.1)
m25.2 <- run_model(
  pop_mod,
  data = data2,
  pars = parsg,
  lower = lowerg,
  upper = upperg
)
fit_check(m25.2)

m25.0$proj
m25.1$proj
m25.2$proj

saveRDS(m25.0, here::here(year, "alt", "results", "m25.0.rds"))
saveRDS(m25.1, here::here(year, "alt", "results", "m25.1.rds"))
saveRDS(m25.2, here::here(year, "alt", "results", "m25.2.rds"))

m25.0r <- run_model_reweight(m25.0, max_one = FALSE)$model
m25.1r <- run_model_reweight(m25.1, max_one = FALSE)$model
m25.2r <- run_model_reweight(m25.2, max_one = FALSE)$model

fit_check(m25.0r)
fit_check(m25.1r)
fit_check(m25.2r)

m25.0r$proj
m25.1r$proj
m25.2r$proj

saveRDS(m25.0r, here::here(year, "alt", "results", "m25.0r.rds"))
saveRDS(m25.1r, here::here(year, "alt", "results", "m25.1r.rds"))
saveRDS(m25.2r, here::here(year, "alt", "results", "m25.2r.rds"))


# Input sample size ----
iss <- c(
  121,
  180,
  114,
  83,
  189,
  205,
  243,
  56,
  130,
  215,
  135,
  218,
  249,
  270,
  164,
  146
) # from surveyISS R package
dataa <- data
data1a <- data1
data2a <- data2
dataa$srv_age_iss <- data1a$srv_age_iss <- data2a$srv_age_iss <- iss

m25.0a <- run_model(
  pop_mod,
  data = dataa,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.0a)
m25.1a <- run_model(
  pop_mod,
  data = data1a,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.1a)
m25.2a <- run_model(
  pop_mod,
  data = data2a,
  pars = parsg,
  lower = lowerg,
  upper = upperg
)
fit_check(m25.2a)

saveRDS(m25.0a, here::here(year, "alt", "results", "m25.0a.rds"))
saveRDS(m25.1a, here::here(year, "alt", "results", "m25.1a.rds"))
saveRDS(m25.2a, here::here(year, "alt", "results", "m25.2a.rds"))


m25.0ar <- run_model_reweight(m25.0a, max_one = FALSE)$model
fit_check(m25.0ar)
m25.1ar <- run_model_reweight(m25.1a, max_one = FALSE)$model
fit_check(m25.1ar)
m25.2ar <- run_model_reweight(m25.2a, max_one = FALSE)$model
fit_check(m25.2ar)

saveRDS(m25.0ar, here::here(year, "alt", "results", "m25.0ar.rds"))
saveRDS(m25.1ar, here::here(year, "alt", "results", "m25.1ar.rds"))
saveRDS(m25.2ar, here::here(year, "alt", "results", "m25.2ar.rds"))

# model-based survey ----
datab <- data
data1b <- data1
data2b <- data2
data2b$srv_obs <- data1b$srv_obs <- datab$srv_obs <- v25$est
data2b$srv_sd <- data1b$srv_sd <- datab$srv_sd <- v25$est *
  sqrt(exp(v25$se^2) - 1) # base model uses sd
data2b$srv_cv <- data1b$srv_cv <- datab$srv_sd / datab$srv_obs

m25.0b <- run_model(
  pop_mod,
  data = datab,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.0b)
m25.1b <- run_model(
  pop_mod,
  data = data1b,
  pars = pars,
  lower = lower,
  upper = upper
)
fit_check(m25.1b)
m25.2b <- run_model(
  pop_mod,
  data = data2b,
  pars = parsg,
  lower = lowerg,
  upper = upperg
)
fit_check(m25.2b)

saveRDS(m25.0b, here::here(year, "alt", "results", "m25.0b.rds"))
saveRDS(m25.1b, here::here(year, "alt", "results", "m25.1b.rds"))
saveRDS(m25.2b, here::here(year, "alt", "results", "m25.2b.rds"))

m25.0br <- run_model_reweight(m25.0b, max_one = FALSE)$model
fit_check(m25.0br)
m25.1br <- run_model_reweight(m25.1b, max_one = FALSE)$model
fit_check(m25.1br)
m25.2br <- run_model_reweight(m25.2b, max_one = FALSE)$model
fit_check(m25.2br)

saveRDS(m25.0br, here::here(year, "alt", "results", "m25.0br.rds"))
saveRDS(m25.1br, here::here(year, "alt", "results", "m25.1br.rds"))
saveRDS(m25.2br, here::here(year, "alt", "results", "m25.2br.rds"))

# iss & mb survey ----
datac <- dataa
data1c <- data1a
data2c <- data2a

data2c$srv_obs <- data1c$srv_obs <- datac$srv_obs <- v25$est
data2c$srv_sd <- data1c$srv_sd <- datac$srv_sd <- v25$est *
  sqrt(exp(v25$se^2) - 1) # base model uses sd
data2c$srv_cv <- data1c$srv_cv <- datac$srv_sd / datac$srv_obs

m25.0c <- run_model(
  pop_mod,
  data = datac,
  pars = pars,
  lower = lower,
  upper = upper
)
m25.1c <- run_model(
  pop_mod,
  data = data1c,
  pars = pars,
  lower = lower,
  upper = upper
)
m25.2c <- run_model(
  pop_mod,
  data = data2c,
  pars = parsg,
  lower = lowerg,
  upper = upperg
)

fit_check(m25.0c)
fit_check(m25.1c)
fit_check(m25.2c)

saveRDS(m25.0c, here::here(year, "alt", "results", "m25.0c.rds"))
saveRDS(m25.1c, here::here(year, "alt", "results", "m25.1c.rds"))
saveRDS(m25.2c, here::here(year, "alt", "results", "m25.2c.rds"))

m25.0cr <- run_model_reweight(m25.0c, max_one = FALSE)$model
fit_check(m25.0cr)
m25.1cr <- run_model_reweight(m25.1c, max_one = FALSE)$model
fit_check(m25.1cr)
m25.2cr <- run_model_reweight(m25.2c, max_one = FALSE)$model
fit_check(m25.2cr)

saveRDS(m25.0cr, here::here(year, "alt", "results", "m25.0cr.rds"))
saveRDS(m25.1cr, here::here(year, "alt", "results", "m25.1cr.rds"))
saveRDS(m25.2cr, here::here(year, "alt", "results", "m25.2cr.rds"))


data.frame(year = data$years, block = data$fish_block_ind)

get_pars(m25.0b, "m25.0b") |>
  left_join(get_pars(m25.1b, "m25.1b")) |>
  left_join(get_pars(m25.2b, "m25.2b"))


plot_ssb(m25.0, m25.0b)
plot_ssb(m25.0b, m25.1b)
plot_ssb(m25.1b, m25.2b)


plot_ssb <- function(m1, m2) {
  tibble(year = m1$rpt$years, ssb = m1$rpt$spawn_bio, id = "m1") |>
    bind_rows(tibble(year = m2$rpt$years, ssb = m2$rpt$spawn_bio, id = "m2")) |>
    ggplot(aes(year, ssb, color = id)) +
    geom_line()
}
plot_ssb(m23, m25.0)
plot_ssb(m25.0, m25.0r)
plot_ssb(m25.0, m25.1)
plot_ssb(m25.0r, m25.1r)
plot_ssb(m25.1r, m25.2)
plot_ssb(m25.1r, m25.2r)

get_pars(m25.0, "m25") %>%
  left_join(get_pars(m25.1, "m25.1"))

get_likes(m25.0, "m25") %>%
  left_join(get_likes(m25.1, "m25.1"))

oo <- resids(m25.0)
ooo <- resids(m25.0, var = "srv_age")
aaa <- resids(m25.0r, var = "srv_age")
oo$agg
ooo$agg

ooo$ss
aaa$agg
aaa$ss
get_pars(m25.2, "m25.2") %>%
  left_join(get_pars(m25.0, "m25.0")) %>%
  left_join(get_pars(m25.1, "m25.1")) %>%
  left_join(get_pars(m23, "m23")) %>%
  print(n = nrow(.))


get_pars(m25.0c, "m25.0c") |>
  left_join(get_pars(m25.1c, "m25.1c")) |>
  left_join(get_pars(m25.2c, "m25.2c"))

plot_ssb(m25.0, m25.2c)

get_likes(m25.1c, "m25.1c", addl = "ssq_catch") |>
  left_join(get_likes(m25.0c, "m25.0c")) |>
  left_join(get_likes(m25.2c, "m25.2c"))
