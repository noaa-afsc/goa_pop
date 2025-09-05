# September 2025 plan team examinations

# load ----
library(RTMButils)
# devtools::unload('RTMButils')
library(tidyverse)
# library(Matrix)
# library(here)
# library(scico)
library(patchwork)
theme_set(afscassess::theme_report())

# source(here::here(2025, 'r', "utils.R"))
source(here::here(2025, 'r', "models.R"))

# data ----
# setup data and pars for all models 
# 2023 data & pars from ADMB model
data = readRDS(here::here(2025, 'rtmb_bridge', "data.rds"))
# pars = readRDS(here::here(2025, 'rtmb_bridge', "pars.rds"))



# original ADMB starting pars
pars = list(log_M = log(0.0614),
            log_a50C = log(c(6, 2.5, 2.5)),
            deltaC = c(1.5,4.5, 4.5),
            log_a50S = log(7.3),
            deltaS = 3.8,
            log_q = log(1.15),
            log_mean_R = 3.0,
            init_log_Rt = rep(0, 26),
            log_Rt = rep(0, sum(data$catch_ind)),
            log_mean_F = 0,
            log_Ft = rep(0, sum(data$catch_ind)),
            log_F35 = 0,
            log_F40 = 0,
            log_F50 = 0,
            sigmaR = 1.7)

# add an a50C and deltaC to estimate 2nd time block selectivity
pars1 = list(log_M = log(0.0614),
            log_a50C = log(c(6, 2.5, 2.5, 2.5)),
            deltaC = c(1.5,4.5, 4.5, 4.5),
            log_a50S = log(7.3),
            deltaS = 3.8,
            log_q = log(1.15),
            log_mean_R = 3.0,
            init_log_Rt = rep(0, 26),
            log_Rt = rep(0, sum(data$catch_ind)),
            log_mean_F = 0,
            log_Ft = rep(0, sum(data$catch_ind)),
            log_F35 = 0,
            log_F40 = 0,
            log_F50 = 0,
            sigmaR = 1.7)

# update catch likelihood data 
data1 = data
data1$catch_wt = 1
data1$catch_cv = rep(0.10, length(data1$years))

# GAP re-stratified data
# data1 = readRDS(here::here(2025, "research", 'gap', 'dat.rds'))

# management models -----
m25.0 = run_model(base, data, pars) # base model run
m25.1 = run_model(srv_like, data1, pars) # update survey biomass likelihood
m25.2 = run_model(slx_scale, data1, pars1) # update selectivity time block
# m25.3 = run_model(ctch_like, data1, pars1) # update selectivity time block

fit_check(m25.0)
fit_check(m25.1)
fit_check(m25.2)

## francis reweight ----
m25.0_rwt = run_model_reweight(model=base, data=data, pars=pars,  iters = 10)
data$fish_age_wt = 2.867722;   data$srv_age_wt = 2.011436; data$fish_size_wt = 0.4984959
m25.0a = run_model(base, data, pars) # base model run

m25.1_rwt = run_model_reweight(model=srv_like, data=data1, pars=pars,  iters = 10)
data1$fish_age_wt = 2.848452; data1$srv_age_wt = 1.998594; data1$fish_size_wt = 0.4985743
m25.1a = run_model(srv_like, data1, pars) # update survey biomass likelihood

m25.2_rwt = run_model_reweight(model=slx_scale, data=data1, pars=pars1,  iters = 10)
data1$fish_age_wt = 2.756654 ; data1$srv_age_wt = 1.990389; data1$fish_size_wt = 0.5251724
m25.2a = run_model(slx_scale, data1, pars1) # update selectivity time block

# reset weights
data$fish_age_wt = data$srv_age_wt = data$fish_size_wt = data1$fish_age_wt = data1$srv_age_wt = data1$fish_size_wt = 1

## mgmt results ----
### tables ----
left_join(get_likes(m25.0, model = 'm25.0'),
          get_likes(m25.0a, model = 'm25.0a')) %>% 
  left_join(get_likes(m25.1, model = 'm25.1', addl = "catch_like") %>% 
              dplyr::mutate(item = ifelse(item=="like_catch", "ssqcatch", item))) %>%
  left_join(get_likes(m25.1a, model = 'm25.1a', addl = "catch_like") %>% 
              dplyr::mutate(item = ifelse(item=="like_catch", "ssqcatch", item))) %>%
  left_join(get_likes(m25.2, model = 'm25.2', addl = "catch_like") %>% 
              dplyr::mutate(item = ifelse(item=="like_catch", "ssqcatch", item))) %>%
  left_join(get_likes(m25.2a, model = 'm25.2a', addl = "catch_like") %>% 
              dplyr::mutate(item = ifelse(item=="like_catch", "ssqcatch", item))) %>%
  mutate(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total")) %>% 
  relocate(Likelihood, m25.0, m25.1, m25.2, m25.0a, m25.1a, m25.2a) %>% 
  select(-item) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "like_tbl_slx.csv"), delim=",")
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::width(j = 1, width = 1.5) 

left_join(get_pars(m25.2, model = 'm25.2'),
          get_pars(m25.0a, model = 'm25.0a')) %>% 
  left_join(get_pars(m25.1, model = 'm25.1')) %>% 
  left_join(get_pars(m25.1a, model = 'm25.1a')) %>%
  left_join(get_pars(m25.0, model = 'm25.0')) %>% 
  left_join(get_pars(m25.2a, model = 'm25.2a')) %>%
  mutate(Item = c("M",  'a50-1', 'a50-2', 'a50-3', 'a50-4', 'delta-1', 'delta-2', 
                  'delta-3', 'delta-4', 'a50 survey', 'delta survey', "q", "sigma R", "Log mean recruitment", 
                  "Log mean F", "2024 Total biomass", "2024 Spawning biomass", "2024 OFL", 
                  "2024 F OFL", " 2024 ABC", "2024 F ABC")) %>% 
  relocate(Item, m25.0, m25.1, m25.2, m25.0a, m25.1a) %>%
  select(-item) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "par_tbl_slx.csv"), delim=",")
  flextable::flextable() %>% 
  flextable::autofit() %>%
  flextable::width(j = 1, width = 2)  %>% 
  flextable::colformat_double(
    i = c(16:18,20),
    big.mark = ",", 
    digits = 0, 
    na_str = "N/A"
  ) %>% 
  flextable::colformat_double(
    i = c(1:15,19,21),
    big.mark = ",", 
    digits = 4, 
    na_str = "N/A"
  ) 

### figs ----
  
#### catch pred ----
m25.0$rpt$catch_pred %>% 
  as.data.frame() %>% 
  mutate(year = m25.0$rpt$years,
         id = "m25.0") %>% 
  bind_rows(
    m25.0a$rpt$catch_pred %>% 
      as.data.frame() %>% 
      mutate(year = m25.0$rpt$years,
             id = "m25.0a")) %>% 
  bind_rows(
    m25.1$rpt$catch_pred %>% 
      as.data.frame() %>% 
      mutate(year = m25.0$rpt$years,
             id = "m25.1")) %>% 
  bind_rows(
    m25.1a$rpt$catch_pred %>% 
      as.data.frame() %>% 
      mutate(year = m25.0$rpt$years,
             id = "m25.1a")) %>% 
  bind_rows(
    m25.2$rpt$catch_pred %>% 
      as.data.frame() %>% 
      mutate(year = m25.0$rpt$years,
             id = "m25.2")) %>% 
  bind_rows(
    m25.2a$rpt$catch_pred %>% 
      as.data.frame() %>% 
      mutate(year = m25.0$rpt$years,
             id = "m25.2a")
  ) %>%
  pivot_longer(-c(year, id)) %>% 
  # mutate(Model = factor(id, levels = c("m25.0", "m25.a", "m25.b", "m25.1", "m25.2", "m25.3"))) %>% 
  # mutate(block = as.numeric(gsub("V", "", name))) %>%
  # filter(!year %in% c(1:35)) %>% 
  ggplot(aes(year, value, color = id)) + 
  geom_line() +
  scico::scale_color_scico_d("Model", palette = 'batlow', end = 0.8) +
  expand_limits(y = 0) +
  theme(legend.position = c(0.8, 0.42)) +
  ylab("Catch (t)") +
  xlab("Year")
ggsave(here::here(2025, 'sep_pt', 'figs', 'catch.png'), units = "in", width=6.5, height=6.5)

#### spawn bio ----
  m25.0$rpt$spawn_bio %>% 
    as.data.frame() %>% 
    mutate(year = m25.0$rpt$years,
           id = "m25.0") %>% 
    bind_rows(
      m25.0a$rpt$spawn_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.0a")) %>% 
    bind_rows(
      m25.1$rpt$spawn_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.1")) %>% 
    bind_rows(
      m25.1a$rpt$spawn_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.1a")) %>% 
    bind_rows(
      m25.2$rpt$spawn_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.2")) %>% 
    bind_rows(
      m25.2a$rpt$spawn_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.2a")
    ) %>%
    pivot_longer(-c(year, id)) %>% 
    # mutate(Model = factor(id, levels = c("m25.0", "m25.a", "m25.b", "m25.1", "m25.2", "m25.3"))) %>% 
    # mutate(block = as.numeric(gsub("V", "", name))) %>%
    # filter(!year %in% c(1:35)) %>% 
    ggplot(aes(year, value, color = id)) + 
    geom_line() +
    scico::scale_color_scico_d("Model", palette = 'batlow', end = 0.8) +
    expand_limits(y = 0) +
    theme(legend.position = c(0.8, 0.22)) +
    ylab("Spawning Biomass (t)") +
    xlab("Year")
  ggsave(here::here(2025, 'sep_pt', 'figs', 'spawn_bio.png'), units = "in", width=6.5, height=6.5)
  
#### tot bio ----
  m25.0$rpt$tot_bio %>% 
    as.data.frame() %>% 
    mutate(year = m25.0$rpt$years,
           id = "m25.0") %>% 
    bind_rows(
      m25.0a$rpt$tot_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.0a")) %>% 
    bind_rows(
      m25.1$rpt$tot_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.1")) %>% 
    bind_rows(
      m25.1a$rpt$tot_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.1a")) %>% 
    bind_rows(
      m25.2$rpt$tot_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.2")) %>% 
    bind_rows(
      m25.2a$rpt$tot_bio %>% 
        as.data.frame() %>% 
        mutate(year = m25.0$rpt$years,
               id = "m25.2a")) %>% 
    pivot_longer(-c(year, id)) %>% 
    # filter(id %in% c("m25", "m25a", "m25b")) %>% 
    # mutate(block = as.numeric(gsub("V", "", name))) %>%
    # filter(!year %in% c(1:35)) %>% 
    ggplot(aes(year, value, color = id)) + 
    geom_line() +
    scico::scale_color_scico_d("Model", palette = 'batlow') +
    expand_limits(y = 0) +
    theme(legend.position = c(0.8, 0.22)) +
    xlab("Year") + 
    ylab("Total biomass (t)")
  ggsave(here::here(2025, 'sep_pt', 'figs', 'tot_bio.png'), units = "in", width=6.5, height=6.5)
  
  
#### survey slx ----  
  m25.0$rpt$slx_srv %>% 
    as.data.frame() %>% 
    mutate(age = 2:29,
           id = "m25.0") %>% 
    bind_rows(
      m25.0a$rpt$slx_srv %>% 
        as.data.frame() %>% 
        mutate(age = 2:29,
               id = "m25.0a")) %>% 
    bind_rows(
      m25.1$rpt$slx_srv %>%
        as.data.frame() %>%
        mutate(age = 2:29,
               id = "m25.1")) %>%
    bind_rows(
      m25.1a$rpt$slx_srv %>%
        as.data.frame() %>%
        mutate(age = 2:29,
               id = "m25.1a")) %>%
    bind_rows(
      m25.2$rpt$slx_srv %>%
        as.data.frame() %>%
        mutate(age = 2:29,
               id = "m25.2")) %>%
    bind_rows(
      m25.2a$rpt$slx_srv %>%
        as.data.frame() %>%
        mutate(age = 2:29,
               id = "m25.2a")
    ) %>%
    pivot_longer(-c(age, id)) %>% 
    # mutate(Model = factor(id, levels = c("m25.0", "m25.a", "m25.b", "m25.1", "m25.2", "m25.3"))) %>% 
    # mutate(year = as.numeric(gsub("V", "", name))) %>%
    # filter(!year %in% c(1:35)) %>% 
    ggplot(aes(age, value, color = id)) + 
    geom_line() +
    scico::scale_color_scico_d("Model", palette = 'batlow') +
    expand_limits(x = 0, y = 0) +
    theme(legend.position = c(0.8, 0.2)) +
    xlab("Age") +
    ylab("Survey Selectivity")
  
  ggsave(here::here(2025, 'sep_pt', 'figs', 'srv_slx.png'), units = "in", width=6.5, height=6.5)  

#### fishery slx ----
  m25.0$rpt$slx_block %>% 
    as.data.frame() %>% 
    rename(B1 = V1, B2 = V2, B3 = V3, B4 = V4) %>% 
    mutate(age = 2:29,
           model = "m25.0") %>% 
    bind_rows(
      m25.0a$rpt$slx_block %>% 
        as.data.frame() %>% 
        rename(B1 = V1, B2 = V2, B3 = V3, B4 = V4) %>% 
        mutate(age = 2:29,
               model = "m25.0a")
    ) %>%
    bind_rows(
      m25.2$rpt$slx_block %>% 
        as.data.frame() %>% 
        rename(B1 = V1, B2 = V2, B3 = V3, B4 = V4) %>% 
        mutate(age = 2:29,
               model = "m25.2")
    ) %>%
    bind_rows(
      m25.2a$rpt$slx_block %>% 
        as.data.frame() %>% 
        rename(B1 = V1, B2 = V2, B3 = V3, B4 = V4) %>% 
        mutate(age = 2:29,
               model = "m25.2a")
    ) %>% 
    pivot_longer(-c(age, model)) %>% 
    # mutate(Model = factor(model, levels = c("m25.0", "m25.a", "m25.b", "m25.1", "m25.2", "m25.3"))) %>% 
    ggplot(aes(age, value, color = name)) + 
    geom_line() +
    scico::scale_color_scico_d("Block", palette = 'roma') +
    facet_wrap(~model) +
    xlab("Age") +
    ylab("Fishery selectivity")
  
  ggsave(here::here(2025, 'sep_pt', 'figs', 'fish_slx_block.png'), units = "in", width=6.5, height=6.5)  
  
  
#### comp residuals ----
  out = resids(obs = data$fish_size_obs, pred = m25.2a$rpt$fish_size_pred, yrs = data$fish_size_yrs, iss=data$fish_size_iss, ind = data$length_bins, label = 'Length (cm)')
  out = resids(obs = data$fish_age_obs, 
               pred = m25.2a$rpt$fish_age_pred, 
               yrs = data$fish_age_yrs, 
               iss=data$fish_age_iss, 
               ind = data$ages, label = 'Age')
  out = resids(obs = data1$srv_age_obs, pred = m25.2a$rpt$srv_age_pred, yrs = data$srv_age_yrs, iss=data1$srv_age_iss, ind = data$ages, label = 'Age')
  
 (out$pearson + out$osa) + plot_layout(guides = "collect")
  ggsave(here::here(2025, 'sep_pt', 'figs', 'm25b_rwt_fish_size_osa.png'), units = "in", width=6.5, height=6.5)  
  

  out$agg / (out$qq + out$ss) + plot_layout(guides = "collect")
  ggsave(here::here(2025, 'sep_pt', 'figs', 'm25b_rwt_fish_size_agg.png'), units = "in", width=6.5, height=6.5)  
  
  
  out$annual
  ggsave(here::here(2025, 'sep_pt', 'figs', 'm25b_rwt_fish_size_annual.png'), units = "in", width=6.5, height=6.5)  
  
#### retro ----
  
  
# GAP restrat ----
  # not presented
  
m25c = run_model(slx_scale, data1, pars1) # base model run
fit_check(m25c)
model_test(m25b, m25c)
### francis reweight ----
m25c_rwt = run_model_reweight(model=slx_scale, data=data1, pars=pars1,  iters = 10)

data2 = data
data2$fish_size_wt = 0
m25d = run_model(slx_scale, data2, pars1) # base model run
fit_check(m25d)
model_test(m25b, m25c)

#### tables ----
left_join(get_likes(m25b, model = 'm25b'),
          get_likes(m25c, model = 'm25c')) %>% 
  left_join(get_likes(m25b_rwt, model = '25b-rwt')) %>%
  left_join(get_likes(m25c_rwt, model = '25c-rwt')) %>%
  mutate(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total")) %>% 
  relocate(Likelihood) %>% 
  select(-item) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "gap_like_tbl.csv"), delim=",")
flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::width(j = 1, width = 1.5) 

left_join(get_pars(m25b, model = 'm25b'),
          get_pars(m25c, model = 'm25c')) %>% 
  left_join(get_pars(m25b_rwt, model = '25b-rwt')) %>%
  left_join(get_pars(m25c_rwt, model = '25c-rwt')) %>%
  mutate(Item = c("M",  'a50-1', 'a50-2', 'a50-3', 'a50-4', 'delta-1', 'delta-2', 
                  'delta-3', 'delta-4', 'a50 survey', 'delta survey', "q", "sigma R", "Log mean recruitment", 
                  "Log mean F", "2024 Total biomass", "2024 Spawning biomass", "2024 OFL", 
                  "2024 F OFL", " 2024 ABC", "2024 F ABC")) %>% 
  relocate(Item) %>% 
  select(-item) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "gap_par_tbl.csv"), delim=",")
flextable::flextable() %>% 
  flextable::autofit() %>%
  flextable::width(j = 1, width = 2)  %>% 
  flextable::colformat_double(
    i = c(16:18,20),
    big.mark = ",", 
    digits = 0, 
    na_str = "N/A"
  ) %>% 
  flextable::colformat_double(
    i = c(1:15,19,21),
    big.mark = ",", 
    digits = 4, 
    na_str = "N/A"
  ) 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# research models ----
# not currently being presented
# ###  data ----
# vast = read_csv("2023/dev/mb_vs_db/vast_2023.csv") %>% 
#   select(year = Time, biomass_mt = Estimate, se = `Std. Error for Estimate`) %>% 
#   filter(biomass_mt>0) %>% 
#   mutate(biomass_mt  = biomass_mt / 1000,
#          se = se / 1000)
#   
# # Vast survey biomass
# data2 = data
# data2$srv_obs = vast$biomass_mt
# data2$srv_sd = vast$se  
# 
# # time varying fishery selectivity
# data3 = data1
# 
# m25c = run_model(slx_scale, data, pars) # base model run
# m25a = run_model(srv_like, data, pars) # update survey biomass likelihood
# m25b = run_model(slx_scale, data, pars1) # update selectivity time block
# 
# fit_check(m25)
# fit_check(m25a)
# fit_check(m25b)
# ## letter or number
# model_test(m25, m25a)
# model_test(m25a, m25b)
# 
# ## francis reweight ----
# m25_rwt = run_model_reweight(model=base, data=data, pars=pars,  iters = 10)
# m25a_rwt = run_model_reweight(model=srv_like, data=data, pars=pars,  iters = 10)
# m25b_rwt = run_model_reweight(model=slx_scale, data=data, pars=pars1,  iters = 10)
# 
# ## mgmt results ----
# ### tables ----
# left_join(get_likes(m25, model = 'm25'),
#           get_likes(m25a, model = 'm25a')) %>% 
#   left_join(get_likes(m25b, model = '25b')) %>%
#   left_join(get_likes(m25_rwt, model = '25-rwt')) %>%
#   left_join(get_likes(m25a_rwt, model = '25a-rwt')) %>%
#   left_join(get_likes(m25b_rwt, model = '25b-rwt')) %>%
#   mutate(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total")) %>% 
#   relocate(Likelihood) %>% 
#   select(-item) %>% 
#   vroom::vroom_write(here::here(2025, "sep_pt", "tables", "like_tbl_slx.csv"), delim=",")
# flextable::flextable() %>%
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 1.5) 
# 
# left_join(get_pars(m25b, model = 'm25b'),
#           get_pars(m25a, model = 'm25a')) %>% 
#   left_join(get_pars(m25, model = 'm25')) %>% 
#   left_join(get_pars(m25_rwt, model = 'm25-rwt')) %>%
#   left_join(get_pars(m25a_rwt, model = 'm25a-rwt')) %>% 
#   left_join(get_pars(m25b_rwt, model = 'm25b-rwt')) %>%
#   mutate(Item = c("M",  'a50-1', 'a50-2', 'a50-3', 'a50-4', 'delta-1', 'delta-2', 
#                   'delta-3', 'delta-4', 'a50 survey', 'delta survey', "q", "sigma R", "Log mean recruitment", 
#                   "Log mean F", "2024 Total biomass", "2024 Spawning biomass", "2024 OFL", 
#                   "2024 F OFL", " 2024 ABC", "2024 F ABC")) %>% 
#   relocate(Item, m25, m25a, m25b) %>% 
#   select(-item) %>% 
#   vroom::vroom_write(here::here(2025, "sep_pt", "tables", "par_tbl_slx.csv"), delim=",")
# flextable::flextable() %>% 
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 2)  %>% 
#   flextable::colformat_double(
#     i = c(14:16,18),
#     big.mark = ",", 
#     digits = 0, 
#     na_str = "N/A"
#   ) %>% 
#   flextable::colformat_double(
#     i = c(1:13,17,19),
#     big.mark = ",", 
#     digits = 4, 
#     na_str = "N/A"
#   ) 
# 
# 
# # merge block 3 & 4
# data2$fish_block_ind[data2$fish_block_ind == 4] <- 3
# n_blk = sum(data2$fish_block_ind==3)
# 
# # data2$cv_M = 0.45
# # adjust pars, add rw pars
# pars2 = pars
# pars2$log_a50C = pars2$log_a50C[-3] # remove last block since is now tv
# # pars2$deltaC = pars2$deltaC[-3] # remove last block since is now tv
# pars2$log_a50C_rw = rep(1, n_blk)
# # pars2$deltaC_rw = rep(5, n_blk)
# pars2$sigma_rw = 0.05
# 
# 
# # time-varying fishery slx and time-varying q
# data3 = data2
# data3$mean_M
# data3$cv_M
# # data3$mean_q = 1
# pars3 = pars2
# pars3$log_q = rep(1, sum(data$srv_ind))
# pars3$sigma_q = 0.05
# 
# # time-varying survey slx
# pars4 = pars
# pars4$log_a50S = rep(pars$log_a50S, sum(data1$srv_ind))
# pars4$sigma_slx_srv = 0.1
# 
# # free M
# data5 = data
# data5$mean_M = 0.0614
# data5$cv_M = 0.3173802  
# m25c = run_model(slx_scale, data1, pars1) # gap re-stratification
# m25c_rwt = run_model_reweight(model=slx_scale, data=data1, pars=pars1,  iters = 10)
# m25.1 = run_model(slx_scale, data2, pars1) # VAST survey biomass
# m25.1_rwt = run_model_reweight(model=slx_scale, data=data2, pars=pars1,  iters = 10)
# 
# fit_check(m25c)
# fit_check(m25c_rwt)
# fit_check(m25.1)
# fit_check(m25.1_rwt)
# 
# model_test(m25b, m25c)
# model_test(m25c, m25.1)
# 
# #### tables ----
# left_join(get_likes(m25c, model = 'm25c'),
#           get_likes(m25c_rwt, model = '25c-rwt')) %>% 
#   left_join(get_likes(m25.2, model = 'm25.2')) %>%
#   left_join(get_likes(m25.2_rwt, model = 'm25.2-rwt')) %>%
#   left_join(get_likes(m25.1, model = 'm25.1')) %>%
#   left_join(get_likes(m25.1_rwt, model = 'm25.1-rwt')) %>%
#   mutate(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total")) %>% 
#   relocate(Likelihood) %>% 
#   select(-item) %>% 
#   # vroom::vroom_write(here::here(2025, "sep_pt", "tables", "like_tbl_slx.csv"), delim=",")
#   flextable::flextable() %>%
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 1.5) 
# 
# 
# left_join(get_pars(m25b, model = 'm25b'),
#           get_pars(m25b_rwt, model = 'm25b-rwt')) %>% 
#   left_join(get_pars(m25c, model = 'm25c')) %>% 
#   left_join(get_pars(m25c_rwt, model = 'm25c-rwt')) %>%
#   left_join(get_pars(m25.2, model = 'm25.2')) %>%
#   left_join(get_pars(m25.2_rwt, model = 'm25.2-rwt')) %>%
#   # left_join(get_pars(m25.9$rpt, m25.9$proj, model = 'm25.9')) %>%
#   # left_join(get_pars(m25.10$rpt, m25.10$proj, model = 'm25.10')) %>%
#   mutate(Item = c("M",  'a50-1', 'a50-2', 'a50-3','a50-4', 'delta-1', 'delta-2', 
#                   'delta-3','delta-4', 'a50 survey', 'delta survey', "q", "sigma R", "Log mean recruitment", 
#                   "Log mean F", "2024 Total biomass", "2024 Spawning biomass", "2024 OFL", 
#                   "2024 F OFL", " 2024 ABC", "2024 F ABC")) %>% 
#   relocate(Item) %>% 
#   select(-item) %>% 
#   # vroom::vroom_write(here::here(2025, "sep_pt", "tables", "par_tbl_slx.csv"), delim=",")
#   flextable::flextable() %>% 
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 2)  %>% 
#   flextable::colformat_double(
#     i = c(16:18,20),
#     big.mark = ",", 
#     digits = 0, 
#     na_str = "N/A"
#   ) %>% 
#   flextable::colformat_double(
#     i = c(1:15,19,21),
#     big.mark = ",", 
#     digits = 4, 
#     na_str = "N/A"
#   ) 
# 
# plot(data$srv_yrs, data$srv_obs, pch=19, col='gray', ylim = c(0, 2100000))
# # lines(data6$srv_yrs, data6$srv_obs, pch=19, type = "p")
# lines(data$srv_yrs, m25$rpt$srv_pred, col=1)
# lines(data$srv_yrs, m25a$rpt$srv_pred, col=2)
# lines(data$srv_yrs, m25b$rpt$srv_pred, col=3)
# lines(data$srv_yrs, m25c$rpt$srv_pred, col=4)
# lines(data$srv_yrs, rwt$rpt$srv_pred, col=5)
# lines(data$srv_yrs, m25.10$rpt$srv_pred, col=6)
# 
# plot(data1$years, m25$rpt$spawn_bio, ylim = c(0, 600000), type="l")
# lines(data$years, m25a$rpt$spawn_bio, col=2)
# lines(data$years, m25b$rpt$spawn_bio, col=3)
# lines(data$years, m25c$rpt$spawn_bio, col=4)
# lines(data$years, rwt$rpt$spawn_bio, col=5)
# lines(data$years, m25.10$rpt$spawn_bio, col=6)
# 
# plot(data1$years, m25$rpt$tot_bio, ylim = c(0, 1800000), type="l")
# lines(data$years, m25a$rpt$tot_bio, col=2)
# lines(data$years, m25b$rpt$tot_bio, col=3)
# lines(data$years, m25c$rpt$tot_bio, col=4)
# lines(data$years, rwt$rpt$tot_bio, col=5)
# lines(data$years, m25.10$rpt$tot_bio, col=6)
# 
# plot(data1$years, m25d$rpt$catch, type="l")
# lines(data1$years, rwt$rpt$catch, type="l", col=5)
# 
# plot(data1$years, m25$rpt$recruits, ylim = c(0, 480), type="l")
# lines(data$years, m25a$rpt$recruits, col=2)
# lines(data$years, m25b$rpt$recruits, col=3)
# lines(data$years, m25c$rpt$recruits, col=4)
# lines(data$years, rwt$rpt$recruits, col=4)
# lines(data$years, m25.10$rpt$recruits, col=6)
# 
# 
# out = resids(obs = data$fish_size_obs, pred = rwt$rpt$fish_size_pred, yrs = data$fish_size_yrs, iss=data$fish_size_iss, ind = data$length_bins, label = 'Length (cm)')
# out = resids(obs = data$fish_age_obs, pred = rwt$rpt$fish_age_pred, yrs = data$fish_age_yrs, iss=data$fish_age_iss, ind = data$ages, label = 'Age')
# out = resids(obs = data1$srv_age_obs, pred = rwt$rpt$srv_age_pred, yrs = data$srv_age_yrs, iss=data1$srv_age_iss, ind = data$ages, label = 'Age')
# 
# 
# m25.8$rpt$q
# m25.8$rpt$sigma_q
# library(patchwork)
# 
# (out$pearson + out$osa) /
#   (out$agg + out$ss) + plot_layout(guides = "collect")
# 
# out$annual
# 
# rwt$rpt$slx_srv %>% 
#   as.data.frame() %>% 
#   mutate(age = 2:29) %>% 
#   pivot_longer(-age) %>% 
#   mutate(year = as.numeric(gsub("V", "", name))) %>%
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(age, value, color = year, group = name)) + 
#   geom_line() +
#   scico::scale_color_scico(palette = 'roma')
# 
# ggsave(here::here(2025, 'sep_pt', 'figs', 'srv_slx.png'), units = "in", width=6.5, height=6.5)
# 
# m25c_rwt$rpt$slx_block %>% 
#   as.data.frame() %>% 
#   # select(V36:63) %>% 
#   mutate(age = 2:29) %>% 
#   pivot_longer(-age) %>% 
#   # mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>% 
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(age, value, color = name)) + 
#   geom_line() +
#   scico::scale_color_scico_d(palette = 'roma')
# 
# m25$rpt$slx_fish %>% 
#   as.data.frame() %>% 
#   # select(V36:63) %>%
#   mutate(age = 2:29) %>% 
#   pivot_longer(-age) %>% 
#   mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>%
#   filter(!year %in% c(1:35)) %>%
#   ggplot(aes(age, value, color = factor(year), group = year)) + 
#   geom_line() +
#   scico::scale_color_scico_d(palette = 'vik')
# 
# 
# 
# # research models 
# m25.2a = run_model(slx_scale_norm, data1, pars)
# model_test(m25.c, m25.2a)
# m25.2b = run_model(slx_scale_log, data1, pars)
# model_test(m25.c, m25.2b)
# m25.3 = run_model(slx_fish_rw, data2, pars2, map = list(sigma_rw = factor(NA)))
# m25.4 = run_model(slx_fish_rw_tvq, data3, pars3, map = list(sigma_q = factor(NA),
#                                                             sigma_rw = factor(NA)))
# m25.4$fit
# m25.5 = run_model(slx_srv, data1, pars4, map = list(sigma_slx_srv = factor(NA)))
# m25.5$fit
# 
# m25.6 = run_model(cohort, data1, pars)
# m25.6$fit
# 
# # free M
# m25.2 = run_model(slx_scale, data5, pars1)
# m25.2_rwt = run_model_reweight(model=slx_scale, data=data5, pars=pars1,  iters = 10)
# model_test(m25.2_rwt, m25.2)
# 
# 
# # tv q
# m25.9 = run_model(tvq, data1, pars,  random = c("log_q", "init_log_Rt", "log_Rt"))
# # data1$cv_M = 0.3
# 
# m25.10 = run_model(slx_scale, data6, pars)
# 
# 
# left_join(get_likes(m25$rpt, model = 'm25'),
#           get_likes(m25a$rpt, model = 'm25a')) %>% 
#   left_join(get_likes(m25.1$rpt, model = 'm25.1')) %>% 
#   left_join(
#     get_likes(m25.2$rpt, model = '25.2')) %>% 
#   left_join(get_likes(m25.2a$rpt, model = '25.2a')) %>%
#   left_join(get_likes(m25.2b$rpt, model = '25.2b')) %>%
#   left_join(get_likes(m25.3$rpt, model = '25.3')) %>%
#   left_join(get_likes(m25.4$rpt, model = '25.4')) %>%
#   left_join(get_likes(m25.5$rpt, model = '25.5')) %>%
#   left_join(get_likes(m25.6$rpt, model = '25.6')) %>%
#   left_join(get_likes(m25.7$rpt, model = '25.7')) %>%
#   mutate(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total")) %>% 
#   relocate(Likelihood) %>% 
#   select(-item) %>% 
#   # vroom::vroom_write(here::here(2025, "sep_pt", "tables", "like_tbl_slx.csv"), delim=",")
#   flextable::flextable() %>%
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 1.5) 
# 
# 
# left_join(get_pars(m25$rpt, m25$proj, model = 'm25'),
#           get_pars(m25.1$rpt, m25.1$proj, model = 'm25.1')) %>% 
#   left_join(get_pars(m25.2$rpt, m25.2$proj, model = 'm25.2')) %>% 
#   left_join(get_pars(m25.2a$rpt, m25.2a$proj, model = 'm25.2a')) %>%
#   left_join(get_pars(m25.2b$rpt, m25.2b$proj, model = 'm25.2b')) %>%
#   left_join(get_pars(m25.3$rpt, m25.3$proj, model = 'm25.3')) %>%
#   left_join(get_pars(m25.4$rpt, m25.4$proj, model = 'm25.4')) %>%
#   left_join(get_pars(m25.5$rpt, m25.5$proj, model = 'm25.5')) %>%
#   left_join(get_pars(m25.6$rpt, m25.6$proj, model = 'm25.6')) %>%
#   left_join(get_pars(m25.7$rpt, m25.7$proj, model = 'm25.7')) %>%
#   mutate(Item = c("M",  'a50-1', 'a50-2', 'a50-3', 'delta-1', 'delta-2', 
#                   'delta-3', 'a50 survey', 'delta survey', "q", "Log mean recruitment", 
#                   "Log mean F", "2024 Total biomass", "2024 Spawning biomass", "2024 OFL", 
#                   "2024 F OFL", " 2024 ABC", "2024 F ABC")) %>% 
#   relocate(Item) %>% 
#   select(-item) %>% 
#   # vroom::vroom_write(here::here(2025, "sep_pt", "tables", "par_tbl_slx.csv"), delim=",")
#   flextable::flextable() %>% 
#   flextable::autofit() %>%
#   flextable::width(j = 1, width = 2)  %>% 
#   flextable::colformat_double(
#     i = c(13:15,17),
#     big.mark = ",", 
#     digits = 0, 
#     na_str = "N/A"
#   ) %>% 
#   flextable::colformat_double(
#     i = c(1:12,16,18),
#     big.mark = ",", 
#     digits = 4, 
#     na_str = "N/A"
#   ) 
# 
# 
# 
# 
# 
# plot(data$srv_yrs, data$srv_obs, pch=19, col='gray')
# lines(data1$srv_yrs, data1$srv_obs, pch=19, type = "p")
# lines(data$srv_yrs, m25$rpt$srv_pred, col=1)
# lines(data$srv_yrs, m25.1$rpt$srv_pred, col=2)
# lines(data$srv_yrs, m25.2$rpt$srv_pred, col=3)
# lines(data$srv_yrs, m25.2a$rpt$srv_pred, col=4)
# # lines(data$srv_yrs, m25.2b$rpt$srv_pred, col=5)
# # lines(data$srv_yrs, m25.3$rpt$srv_pred, col=6)
# # lines(data$srv_yrs, m25.4$rpt$srv_pred, col=7)
# lines(data$srv_yrs, m25.5$rpt$srv_pred, col=7)
# # lines(data$srv_yrs, m25.6$rpt$srv_pred, col=7)
# # lines(data$srv_yrs, m25.7$rpt$srv_pred, col=7)
# 
# plot(data1$years, m25$rpt$spawn_bio, ylim = c(0, 400000), type="l")
# lines(data$years, m25a$rpt$spawn_bio, col=2)
# lines(data$years, m25b$rpt$spawn_bio, col=3)
# lines(data$years, m25b_rwt$rpt$spawn_bio, col=4)
# lines(data$years, m25c$rpt$spawn_bio, col=5)
# lines(data$years, m25c_rwt$rpt$spawn_bio, col=6)
# # lines(data$years, m25.3$rpt$spawn_bio, col=7)
# # lines(data$years, m25.4$rpt$spawn_bio, col=8)
# lines(data$years, m25.5$rpt$spawn_bio, col=9)
# # lines(data$years, m25.6$rpt$spawn_bio, col=10)
# # lines(data$years, m25.7$rpt$spawn_bio, col=11)
# 
# plot(data1$years, m25$rpt$Ft, type="l")
# lines(data$years, m25a$rpt$Ft, col=2)
# lines(data$years, m25b$rpt$Ft, col=3)
# lines(data$years, m25b_rwt$rpt$Ft, col=4)
# lines(data$years, m25c$rpt$Ft, col=5)
# lines(data$years, m25c_rwt$rpt$Ft, col=6)
# 
# plot(data1$years, m25$rpt$tot_bio, ylim = c(0, 1800000), type="l")
# lines(data$years, m25a$rpt$tot_bio, col=2)
# lines(data$years, m25.1$rpt$tot_bio, col=3)
# lines(data$years, m25.2$rpt$tot_bio, col=4)
# # lines(data$years, m25.2a$rpt$tot_bio, col=5)
# # lines(data$years, m25.2b$rpt$tot_bio, col=6)
# # lines(data$years, m25.3$rpt$tot_bio, col=7)
# # lines(data$years, m25.4$rpt$tot_bio, col=8)
# lines(data$years, m25.5$rpt$tot_bio, col=9)
# # lines(data$years, m25.6$rpt$tot_bio, col=10)
# # lines(data$years, m25.7$rpt$tot_bio, col=11)
# 
# 
# plot(data1$years, m25$rpt$recruits, ylim = c(0, 1000), type="l")
# lines(data$years, m25a$rpt$recruits, col=2)
# lines(data$years, m25.1$rpt$recruits, col=3)
# lines(data$years, m25.2$rpt$recruits, col=4)
# # lines(data$years, m25.2a$rpt$recruits, col=5)
# # lines(data$years, m25.2b$rpt$recruits, col=6)
# # lines(data$years, m25.3$rpt$recruits, col=7)
# # lines(data$years, m25.4$rpt$recruits, col=8)
# lines(data$years, m25.5$rpt$recruits, col=9)
# # lines(data$years, m25.6$rpt$recruits, col=10)
# # lines(data$years, m25.7$rpt$recruits, col=11)
# 
# 
# 
# m25$rpt$slx_block %>% 
#   as.data.frame() %>% 
#   mutate(age = 2:29,
#          id = "m25") %>% 
#   bind_rows(
#     m25a$rpt$slx_block %>% 
#       as.data.frame() %>% 
#       mutate(age = 2:29,
#              id = "m25a")) %>% 
#   bind_rows(
#     m25_rwt$rpt$slx_block %>% 
#       as.data.frame() %>% 
#       mutate(age = 2:29,
#              id = "m25-rwt")) %>% 
#   bind_rows(
#     m25a_rwt$rpt$slx_block %>% 
#       as.data.frame() %>% 
#       mutate(age = 2:29,
#              id = "m25a-rwt")) %>% 
#   bind_rows(
#     m25b_rwt$rpt$slx_block %>% 
#       as.data.frame() %>% 
#       mutate(age = 2:29,
#              id = "m25b-rwt")) %>% 
#   bind_rows(
#     m25b$rpt$slx_block %>% 
#       as.data.frame() %>% 
#       mutate(age = 2:29,
#              id = "m25b")
#   ) %>% 
#   pivot_longer(-c(age, id)) %>% 
#   mutate(block = as.numeric(gsub("V", "", name))) %>%
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(age, value, color = id, group = interaction(id, block))) + 
#   geom_line() +
#   scico::scale_color_scico_d("Model", palette = 'romaO') +
#   expand_limits(x = 0, y = 0) +
#   facet_wrap(~block) +
#   theme(legend.position = c(0.8, 0.22))
#   ggsave(here::here(2025, 'sep_pt', 'figs', 'fish_slx.png'), units = "in", width=6.5, height=6.5)
# 
# 
# m25b_rwt$rpt$slx_block %>% 
#   as.data.frame() %>% 
#   # select(V36:63) %>% 
#   mutate(age = 2:29) %>% 
#   pivot_longer(-age) %>% 
#   # mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>% 
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(age, value, color = name)) + 
#   geom_line() +
#   scico::scale_color_scico_d(palette = 'roma')
# 
# m25.4$rpt$slx_fish %>% 
#   as.data.frame() %>% 
#   select(V36:63) %>%
#   mutate(age = 2:29) %>% 
#   pivot_longer(-age) %>% 
#   mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>%
#   filter(!year %in% c(1:35)) %>%
#   ggplot(aes(age, value, color = factor(year), group = year)) + 
#   geom_line() +
#   scico::scale_color_scico_d(palette = 'vik')
# 
# 
# # figs ----
# 
# # spawn bio ----
# m25$rpt$spawn_bio %>% 
#   as.data.frame() %>% 
#   mutate(year = m25$rpt$years,
#          id = "m25") %>% 
#   bind_rows(
#     m25a$rpt$spawn_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25a")) %>% 
#   bind_rows(
#     m25_rwt$rpt$spawn_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25-rwt")) %>% 
#   bind_rows(
#     m25a_rwt$rpt$spawn_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25a-rwt")) %>% 
#   bind_rows(
#     m25b_rwt$rpt$spawn_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25b-rwt")) %>% 
#   bind_rows(
#     m25b$rpt$spawn_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25b")
#   ) %>% 
#   pivot_longer(-c(year, id)) %>% 
#   # mutate(block = as.numeric(gsub("V", "", name))) %>%
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(year, value, color = id)) + 
#   geom_line() +
#   scico::scale_color_scico_d("Model", palette = 'romaO') +
#   expand_limits(y = 0) +
#   theme(legend.position = c(0.8, 0.22))
# ggsave(here::here(2025, 'sep_pt', 'figs', 'spawn_bio.png'), units = "in", width=6.5, height=6.5)
# 
# # tot bio ----
# m25$rpt$tot_bio %>% 
#   as.data.frame() %>% 
#   mutate(year = m25$rpt$years,
#          id = "m25") %>% 
#   bind_rows(
#     m25a$rpt$tot_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25a")) %>% 
#   bind_rows(
#     m25_rwt$rpt$tot_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25-rwt")) %>% 
#   bind_rows(
#     m25a_rwt$rpt$tot_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25a-rwt")) %>% 
#   bind_rows(
#     m25b_rwt$rpt$tot_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25b-rwt")) %>% 
#   bind_rows(
#     m25b$rpt$tot_bio %>% 
#       as.data.frame() %>% 
#       mutate(year = m25$rpt$years,
#              id = "m25b")
#   ) %>% 
#   pivot_longer(-c(year, id)) %>% 
#   filter(id %in% c("m25", "m25a", "m25b")) %>% 
#   # mutate(block = as.numeric(gsub("V", "", name))) %>%
#   # filter(!year %in% c(1:35)) %>% 
#   ggplot(aes(year, value, color = id)) + 
#   geom_line() +
#   scico::scale_color_scico_d("Model", palette = 'romaO') +
#   expand_limits(y = 0) +
#   theme(legend.position = c(0.8, 0.22))
# ggsave(here::here(2025, 'sep_pt', 'figs', 'tot_bio.png'), units = "in", width=6.5, height=6.5)
# 
# 
