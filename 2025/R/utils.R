# have to update the data$catch_wt with the change in likelihood
run_retro_update <- function(rpt, data, pars, model, n_peels = 10, year, folder, subfolder = NULL,
                             quantities = c("spawn_bio", "tot_bio", "recruits"),
                             map = list(sigmaR = factor(NA)),
                             control = list(iter.max = 100000, eval.max = 20000),
                             save_outputs = TRUE) {
  
  # Setup and Load Base Model Components ---
  message("Starting retrospective analysis...")
  retro_path <- herein(year, folder, subfolder, "retro")
  if (save_outputs && !dir.exists(retro_path)) {
    dir.create(retro_path, recursive = TRUE)
  }
  
  d0 <- data
  p0 <- pars
  model <- model
  
  n_years <- length(d0$srv_ind)
  
  
  # Run Retrospective Peels ---
  message(paste("Running", n_peels, "retrospective peels..."))
  reps <- list()
  N_base_sum_catch <- sum(d0$catch_ind)
  
  for (i in 1:n_peels) {
    data <- d0
    pars <- p0
    
    peel_indices <- (n_years - i + 1):n_years
    
    # Peel data based on the logic in the original script
    data$years <- head(d0$years, -i)
    data$catch_ind <- head(d0$catch_ind, -i)
    data$catch_obs <- head(d0$catch_obs, -i)
    # data$catch_wt <- head(d0$catch_wt, -i)
    
    data$srv_ind[peel_indices] <- 0
    data$fish_age_ind[peel_indices] <- 0
    data$srv_age_ind[peel_indices] <- 0
    data$fish_size_ind[peel_indices] <- 0
    
    data <- lapply(data, unname)
    
    # Peel parameters
    pars$log_Ft <- head(p0$log_Ft, -i)
    pars$log_Rt <- head(p0$log_Rt, -i)
    
    pars <- lapply(pars, unname)
    
    # Refit model for the peel
    obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, silent = TRUE)
    fit <- nlminb(start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = control)
    reps[[paste0('sd', i)]] <- sdreport(obj)
  }
  
  if (save_outputs) {
    saveRDS(reps, file.path(retro_path, 'reps.RDS'))
    message(paste("Retrospective fits saved to:", file.path(retro_path, 'reps.RDS')))
  }
  
  # --- 3. Process and Tidy All Outputs ---
  message("Processing results and calculating Mohn's rho...")
  # Get sdreport summary from the base model fit
  base_sds <- as.data.frame(summary(sdreport(rpt$obj)))
  base_sds <- tibble::rownames_to_column(base_sds, "item")
  
  # Process the retrospective fits
  retro_sds_list <- purrr::map(reps, summary) %>%
    purrr::map(., as.data.frame) %>%
    purrr::map(., tibble::rownames_to_column, "item")
  
  # Tidy the retrospective data
  retro_df <- dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(
      peel = as.numeric(gsub("sd", "", id)),
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Tidy the base model data
  base_df <- base_sds %>%
    dplyr::mutate(
      peel = 0,
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Combine and add years
  base_years <- rpt$rpt$years
  all_data <- dplyr::bind_rows(retro_df, base_df) %>%
    dplyr::group_by(peel, item) %>%
    dplyr::mutate(year = base_years[1:dplyr::n()]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lci = Estimate - 1.96 * Std_Error,
      uci = Estimate + 1.96 * Std_Error,
      lci = ifelse(lci < 0, 0, lci),
      peel = factor(peel, levels = as.character(0:n_peels))
    )
  
  # --- 4. Calculate Mohn's Rho and MAE ---
  terminal_peels <- all_data %>%
    dplyr::filter(peel != 0) %>%
    dplyr::group_by(item, peel) %>%
    dplyr::summarise(year = max(year), peel_est = Estimate, .groups = 'drop')
  
  terminal_base <- all_data %>%
    dplyr::filter(peel == 0) %>%
    dplyr::select(item, year, base_est = Estimate)
  
  rho_metrics <- dplyr::left_join(terminal_peels, terminal_base, by = c("item", "year")) %>%
    dplyr::mutate(pdiff = (peel_est - base_est) / base_est) %>% # Using your formula
    tidyr::drop_na() %>%
    dplyr::group_by(item) %>%
    dplyr::summarise(
      rho = mean(pdiff, na.rm = TRUE),
      mae = median(abs(pdiff), na.rm = TRUE)
    )
  
  # --- 5. Generate and Save Plots ---
  message("Generating plots...")
  plot_list <- list()
  
  for (qty in quantities) {
    metrics <- rho_metrics %>% dplyr::filter(item == qty)
    l1 <- paste("Mohn's rho =", round(metrics$rho, 2))
    l2 <- paste("MAE =", round(metrics$mae, 2))
    plot_title <- gsub("_", " ", qty)
    plot_title <- paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))
    
    plot_data_qty <- all_data %>% dplyr::filter(item == qty)
    annot_x <- min(plot_data_qty$year) + 1 # Place annotation near the start
    annot_y <- max(plot_data_qty$uci, na.rm = TRUE) * 0.9 # Place it high up
    
    # Plot 1: Time series
    p1 <- ggplot2::ggplot(plot_data_qty, ggplot2::aes(x = year, y = Estimate, color = peel, group = peel, fill = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      scico::scale_fill_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::annotate(geom = 'text', x = annot_x, y = annot_y, label = paste(l1, l2, sep = "\n"), hjust = 0) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = paste(plot_title, "(t)"), x = "Year", title = paste("Retrospective Pattern:", plot_title)) +
      ggplot2::expand_limits(y = 0)
    
    # Plot 2: Relative difference
    pdiff_data <- all_data %>%
      dplyr::filter(item == qty, peel != 0) %>%
      dplyr::left_join(
        all_data %>% dplyr::filter(item == qty, peel == 0) %>% dplyr::select(year, base_est = Estimate),
        by = "year"
      ) %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est) %>% # Standard relative diff for plotting
      tidyr::drop_na()
    
    p2 <- ggplot2::ggplot(pdiff_data, ggplot2::aes(x = year, y = pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE, linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "Relative Difference from Base", x = "Year")
    
    # Combine and save
    combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1.5)) +
      patchwork::plot_layout(guides = "collect")
    
    plot_list[[qty]] <- combined_plot
    
    if (save_outputs) {
      file_name <- paste0(gsub(" ", "_", tolower(plot_title)), "_retro.png")
      ggplot2::ggsave(
        filename = file.path(retro_path, file_name),
        plot = combined_plot,
        width = 7, height = 7, units = "in", dpi = 300
      )
    }
  }
  
  message("Done.")
  # --- 6. Return Results ---
  return(list(
    retro_data = all_data,
    rho_metrics = rho_metrics,
    plots = plot_list
  ))
}

run_prospective_update <- function(rpt, data, pars, model, n_peels = 10, year, folder, subfolder = NULL,
                            quantities = c("spawn_bio", "tot_bio", "recruits"),
                            map = list(sigmaR = factor(NA)),
                            control = list(iter.max = 100000, eval.max = 20000),
                            save_outputs = TRUE) {
  
  # Setup and Load Base Model Components ---
  message("Starting prospective analysis...")
  retro_path <- herein(year, folder, subfolder, "prospective")
  if (save_outputs && !dir.exists(retro_path)) {
    dir.create(retro_path, recursive = TRUE)
  }
  
  d0 <- data
  p0 <- pars
  model <- model
  
  n_years <- length(data$srv_ind)
  
  
  # Run Retrospective Peels ---
  message(paste("Running", n_peels, "Prospective peels..."))
  reps <- list()
  N_base_sum_catch <- sum(d0$catch_ind)
  
  for (i in 1:n_peels) {
    data <- d0
    pars <- p0
    peel_indices <- (n_years - i + 1):n_years
    
    # Peel data based on the logic in the original script
    data$years <- tail(d0$years, -i)
    data$catch_ind <- tail(d0$catch_ind, -i)
    data$catch_obs <- tail(d0$catch_obs, -i)
    # data$catch_wt <- tail(d0$catch_wt, -i)
    
    data$srv_ind = tail(d0$srv_ind, -i)
    data$fish_age_ind = tail(d0$fish_age_ind, -i)
    data$srv_age_ind = tail(d0$srv_age_ind, -i)
    data$fish_size_ind = tail(d0$fish_size_ind, -i)
    
    data <- lapply(data, unname)
    
    # Peel parameters
    pars$log_Ft <- tail(p0$log_Ft, -i)
    pars$log_Rt <- tail(p0$log_Rt, -i)
    
    pars <- lapply(pars, unname)
    
    # Refit model for the peel
    obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, silent = TRUE)
    fit <- nlminb(start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = control)
    reps[[paste0('sd', i)]] <- sdreport(obj)
  }
  
  if (save_outputs) {
    saveRDS(reps, file.path(retro_path, 'reps.RDS'))
    message(paste("Retrospective fits saved to:", file.path(retro_path, 'reps.RDS')))
  }
  
  # --- 3. Process and Tidy All Outputs ---
  message("Processing results and calculating Mohn's rho...")
  # Get sdreport summary from the base model fit
  base_sds <- as.data.frame(summary(sdreport(rpt$obj)))
  base_sds <- tibble::rownames_to_column(base_sds, "item")
  
  # Process the retrospective fits
  retro_sds_list <- purrr::map(reps, summary) %>%
    purrr::map(., as.data.frame) %>%
    purrr::map(., tibble::rownames_to_column, "item")
  
  # Tidy the retrospective data
  retro_df <- dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(
      peel = as.numeric(gsub("sd", "", id)),
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Tidy the base model data
  base_df <- base_sds %>%
    dplyr::mutate(
      peel = 0,
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Combine and add years
  base_years <- rpt$rpt$years
  all_data <- dplyr::bind_rows(retro_df, base_df) %>%
    dplyr::group_by(peel, item) %>%
    dplyr::mutate(year = base_years[(peel + 1):(peel + dplyr::n())]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lci = Estimate - 1.96 * Std_Error,
      uci = Estimate + 1.96 * Std_Error,
      lci = ifelse(lci < 0, 0, lci),
      peel = factor(peel, levels = as.character(0:n_peels))
    )
  
  # Calculate Mohn's Rho and MAE ---
  # terminal_peels <- all_data %>%
  #   dplyr::filter(peel != 0) %>%
  #   dplyr::group_by(item, peel) %>%
  #   dplyr::summarise(year = max(year), peel_est = Estimate, .groups = 'drop')
  
  # terminal_base <- all_data %>%
  #   dplyr::filter(peel == 0) %>%
  #   dplyr::select(item, year, base_est = Estimate)
  
  # rho_metrics <- dplyr::left_join(terminal_peels, terminal_base, by = c("item", "year")) %>%
  #   dplyr::mutate(pdiff = (peel_est - base_est) / base_est) %>% # Using your formula
  #   tidyr::drop_na() %>%
  #   dplyr::group_by(item) %>%
  #   dplyr::summarise(
  #     rho = mean(pdiff, na.rm = TRUE),
  #     mae = median(abs(pdiff), na.rm = TRUE)
  #   )
  
  # --- 5. Generate and Save Plots ---
  message("Generating plots...")
  plot_list <- list()
  
  for (qty in quantities) {
    # metrics <- rho_metrics %>% dplyr::filter(item == qty)
    # l1 <- paste("Mohn's rho =", round(metrics$rho, 2))
    # l2 <- paste("MAE =", round(metrics$mae, 2))
    plot_title <- gsub("_", " ", qty)
    plot_title <- paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))
    
    plot_data_qty <- all_data %>% dplyr::filter(item == qty)
    annot_x <- min(plot_data_qty$year) + 1 # Place annotation near the start
    annot_y <- max(plot_data_qty$uci, na.rm = TRUE) * 0.9 # Place it high up
    
    # Plot 1: Time series
    p1 <- ggplot2::ggplot(plot_data_qty, ggplot2::aes(x = year, y = Estimate, color = peel, group = peel, fill = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      scico::scale_fill_scico_d("Peel", palette = 'roma', direction = -1) +
      # ggplot2::annotate(geom = 'text', x = annot_x, y = annot_y, label = paste(l1, l2, sep = "\n"), hjust = 0) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = paste(plot_title, "(t)"), x = "Year", title = paste("Prospective Pattern:", plot_title)) +
      ggplot2::expand_limits(y = 0)
    
    # Plot 2: Relative difference
    pdiff_data <- all_data %>%
      dplyr::filter(item == qty, peel != 0) %>%
      dplyr::left_join(
        all_data %>% dplyr::filter(item == qty, peel == 0) %>% dplyr::select(year, base_est = Estimate),
        by = "year"
      ) %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est) %>% # Standard relative diff for plotting
      tidyr::drop_na()
    
    p2 <- ggplot2::ggplot(pdiff_data, ggplot2::aes(x = year, y = pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE, linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "Relative Difference from Base", x = "Year")
    
    # Combine and save
    combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1.5)) +
      patchwork::plot_layout(guides = "collect")
    
    plot_list[[qty]] <- combined_plot
    
    if (save_outputs) {
      file_name <- paste0(gsub(" ", "_", tolower(plot_title)), "_pro.png")
      ggplot2::ggsave(
        filename = file.path(retro_path, file_name),
        plot = combined_plot,
        width = 7, height = 7, units = "in", dpi = 300
      )
    }
  }
  
  message("Done.")
  # --- 6. Return Results ---
  return(list(
    prospect_data = all_data,
    # rho_metrics = rho_metrics,
    plots = plot_list
  ))
}

#' Run a retrospective analysis (peeling) for an RTMB stock assessment model.
#'
#' This function automates the process of retrospective analysis. It iteratively
#' removes recent years of data, refits the model, calculates Mohn's rho,
#' and generates diagnostic plots for specified quantities of interest.
#'
#' @param rpt The output list from the `run_model()` function for the
#'   full dataset. Must contain `obj` (from RTMB::MakeADFun) and `rpt` (the report)
#' @param n_peels The number of retrospective peels to perform (integer)
#' @param year The assessment year, used for constructing file paths (e.g., 2024)
#' @param folder The folder name for the model or type (e.g., "mgmt", "research"), used for file paths (e.g., "mgmt")
#' @param subfolder used for file paths (e.g., "m24"), default: NULL
#' @param quantities A character vector of variable names from the model report
#'   to analyze (e.g., c("spawn_bio", "tot_bio", "recruits"))
#' @param map A named list specifying the mapping for `RTMB::MakeADFun`
#' @param control A list of control parameters for the `nlminb` optimizer
#' @param save_outputs A logical flag. If TRUE, saves the retrospective results
#'   (.RDS) and plots (.png) to a 'retro' subfolder
#'
#' @return A list containing three elements:
#'   1. `retro_data`: A tidy data frame with all estimates from all peels.
#'   2. `rho_metrics`: A data frame with Mohn's rho and MAE for each quantity.
#'   3. `plots`: A named list of the generated ggplot objects.
#'
#' @details
#' The function relies on the following packages: `RTMB`, `dplyr`, `purrr`,
#' `ggplot2`, `patchwork`, `scico`, `tibble`, and `here`. Ensure they are installed.
#' The Mohn's rho is calculated as the mean of `(Peel_terminal - Base_terminal) / Base_terminal`.
#'
#' @export
run_retro <- function(rpt, data, pars, model, n_peels = 10, year, folder, subfolder = NULL,
                      quantities = c("spawn_bio", "tot_bio", "recruits"),
                      map = list(sigmaR = factor(NA)),
                      control = list(iter.max = 100000, eval.max = 20000),
                      save_outputs = TRUE) {
  
  # Setup and Load Base Model Components ---
  message("Starting retrospective analysis...")
  retro_path <- herein(year, folder, subfolder, "retro")
  if (save_outputs && !dir.exists(retro_path)) {
    dir.create(retro_path, recursive = TRUE)
  }
  
  d0 <- data
  p0 <- pars
  model <- model
  
  n_years <- length(data$srv_ind)
  
  
  # Run Retrospective Peels ---
  message(paste("Running", n_peels, "retrospective peels..."))
  reps <- list()
  N_base_sum_catch <- sum(d0$catch_ind)
  
  for (i in 1:n_peels) {
    data <- d0
    pars <- p0
    
    peel_indices <- (n_years - i + 1):n_years
    
    # Peel data based on the logic in the original script
    data$years <- head(d0$years, -i)
    data$catch_ind <- head(d0$catch_ind, -i)
    data$catch_obs <- head(d0$catch_obs, -i)
    data$catch_wt <- head(d0$catch_wt, -i)
    
    data$srv_ind[peel_indices] <- 0
    data$fish_age_ind[peel_indices] <- 0
    data$srv_age_ind[peel_indices] <- 0
    data$fish_size_ind[peel_indices] <- 0
    
    data <- lapply(data, unname)
    
    # Peel parameters
    pars$log_Ft <- head(p0$log_Ft, -i)
    pars$log_Rt <- head(p0$log_Rt, -i)
    
    pars <- lapply(pars, unname)
    
    # Refit model for the peel
    obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, silent = TRUE)
    fit <- nlminb(start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = control)
    reps[[paste0('sd', i)]] <- sdreport(obj)
  }
  
  if (save_outputs) {
    saveRDS(reps, file.path(retro_path, 'reps.RDS'))
    message(paste("Retrospective fits saved to:", file.path(retro_path, 'reps.RDS')))
  }
  
  # --- 3. Process and Tidy All Outputs ---
  message("Processing results and calculating Mohn's rho...")
  # Get sdreport summary from the base model fit
  base_sds <- as.data.frame(summary(sdreport(rpt$obj)))
  base_sds <- tibble::rownames_to_column(base_sds, "item")
  
  # Process the retrospective fits
  retro_sds_list <- purrr::map(reps, summary) %>%
    purrr::map(., as.data.frame) %>%
    purrr::map(., tibble::rownames_to_column, "item")
  
  # Tidy the retrospective data
  retro_df <- dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(
      peel = as.numeric(gsub("sd", "", id)),
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Tidy the base model data
  base_df <- base_sds %>%
    dplyr::mutate(
      peel = 0,
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Combine and add years
  base_years <- rpt$rpt$years
  all_data <- dplyr::bind_rows(retro_df, base_df) %>%
    dplyr::group_by(peel, item) %>%
    dplyr::mutate(year = base_years[1:dplyr::n()]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lci = Estimate - 1.96 * Std_Error,
      uci = Estimate + 1.96 * Std_Error,
      lci = ifelse(lci < 0, 0, lci),
      peel = factor(peel, levels = as.character(0:n_peels))
    )
  
  # --- 4. Calculate Mohn's Rho and MAE ---
  terminal_peels <- all_data %>%
    dplyr::filter(peel != 0) %>%
    dplyr::group_by(item, peel) %>%
    dplyr::summarise(year = max(year), peel_est = Estimate, .groups = 'drop')
  
  terminal_base <- all_data %>%
    dplyr::filter(peel == 0) %>%
    dplyr::select(item, year, base_est = Estimate)
  
  rho_metrics <- dplyr::left_join(terminal_peels, terminal_base, by = c("item", "year")) %>%
    dplyr::mutate(pdiff = (peel_est - base_est) / base_est) %>% # Using your formula
    tidyr::drop_na() %>%
    dplyr::group_by(item) %>%
    dplyr::summarise(
      rho = mean(pdiff, na.rm = TRUE),
      mae = median(abs(pdiff), na.rm = TRUE)
    )
  
  # --- 5. Generate and Save Plots ---
  message("Generating plots...")
  plot_list <- list()
  
  for (qty in quantities) {
    metrics <- rho_metrics %>% dplyr::filter(item == qty)
    l1 <- paste("Mohn's rho =", round(metrics$rho, 2))
    l2 <- paste("MAE =", round(metrics$mae, 2))
    plot_title <- gsub("_", " ", qty)
    plot_title <- paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))
    
    plot_data_qty <- all_data %>% dplyr::filter(item == qty)
    annot_x <- min(plot_data_qty$year) + 1 # Place annotation near the start
    annot_y <- max(plot_data_qty$uci, na.rm = TRUE) * 0.9 # Place it high up
    
    # Plot 1: Time series
    p1 <- ggplot2::ggplot(plot_data_qty, ggplot2::aes(x = year, y = Estimate, color = peel, group = peel, fill = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      scico::scale_fill_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::annotate(geom = 'text', x = annot_x, y = annot_y, label = paste(l1, l2, sep = "\n"), hjust = 0) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = paste(plot_title, "(t)"), x = "Year", title = paste("Retrospective Pattern:", plot_title)) +
      ggplot2::expand_limits(y = 0)
    
    # Plot 2: Relative difference
    pdiff_data <- all_data %>%
      dplyr::filter(item == qty, peel != 0) %>%
      dplyr::left_join(
        all_data %>% dplyr::filter(item == qty, peel == 0) %>% dplyr::select(year, base_est = Estimate),
        by = "year"
      ) %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est) %>% # Standard relative diff for plotting
      tidyr::drop_na()
    
    p2 <- ggplot2::ggplot(pdiff_data, ggplot2::aes(x = year, y = pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE, linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "Relative Difference from Base", x = "Year")
    
    # Combine and save
    combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1.5)) +
      patchwork::plot_layout(guides = "collect")
    
    plot_list[[qty]] <- combined_plot
    
    if (save_outputs) {
      file_name <- paste0(gsub(" ", "_", tolower(plot_title)), "_retro.png")
      ggplot2::ggsave(
        filename = file.path(retro_path, file_name),
        plot = combined_plot,
        width = 7, height = 7, units = "in", dpi = 300
      )
    }
  }
  
  message("Done.")
  # --- 6. Return Results ---
  return(list(
    retro_data = all_data,
    rho_metrics = rho_metrics,
    plots = plot_list
  ))
}

run_prospective <- function(rpt, data, pars, model, n_peels = 10, year, folder, subfolder = NULL,
                            quantities = c("spawn_bio", "tot_bio", "recruits"),
                            map = list(sigmaR = factor(NA)),
                            control = list(iter.max = 100000, eval.max = 20000),
                            save_outputs = TRUE) {
  
  # Setup and Load Base Model Components ---
  message("Starting prospective analysis...")
  retro_path <- herein(year, folder, subfolder, "prospective")
  if (save_outputs && !dir.exists(retro_path)) {
    dir.create(retro_path, recursive = TRUE)
  }
  
  d0 <- data
  p0 <- pars
  model <- model
  
  n_years <- length(data$srv_ind)
  
  
  # Run Retrospective Peels ---
  message(paste("Running", n_peels, "prospective peels..."))
  reps <- list()
  N_base_sum_catch <- sum(d0$catch_ind)
  
  for (i in 1:n_peels) {
    data <- d0
    pars <- p0
    peel_indices <- (n_years - i + 1):n_years
    
    # Peel data based on the logic in the original script
    data$years <- tail(d0$years, -i)
    data$catch_ind <- tail(d0$catch_ind, -i)
    data$catch_obs <- tail(d0$catch_obs, -i)
    data$catch_wt <- tail(d0$catch_wt, -i)
    
    data$srv_ind = tail(d0$srv_ind, -i)
    data$fish_age_ind = tail(d0$fish_age_ind, -i)
    data$srv_age_ind = tail(d0$srv_age_ind, -i)
    data$fish_size_ind = tail(d0$fish_size_ind, -i)
    
    data <- lapply(data, unname)
    
    # Peel parameters
    pars$log_Ft <- tail(p0$log_Ft, -i)
    pars$log_Rt <- tail(p0$log_Rt, -i)
    
    pars <- lapply(pars, unname)
    
    # Refit model for the peel
    obj <- RTMB::MakeADFun(cmb(model, data), pars, map = map, silent = TRUE)
    fit <- nlminb(start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = control)
    reps[[paste0('sd', i)]] <- sdreport(obj)
  }
  
  if (save_outputs) {
    saveRDS(reps, file.path(retro_path, 'reps.RDS'))
    message(paste("Prrospective fits saved to:", file.path(retro_path, 'reps.RDS')))
  }
  
  # --- 3. Process and Tidy All Outputs ---
  message("Processing results and calculating Mohn's rho...")
  # Get sdreport summary from the base model fit
  base_sds <- as.data.frame(summary(sdreport(rpt$obj)))
  base_sds <- tibble::rownames_to_column(base_sds, "item")
  
  # Process the retrospective fits
  retro_sds_list <- purrr::map(reps, summary) %>%
    purrr::map(., as.data.frame) %>%
    purrr::map(., tibble::rownames_to_column, "item")
  
  # Tidy the retrospective data
  retro_df <- dplyr::bind_rows(retro_sds_list, .id = 'id') %>%
    dplyr::mutate(
      peel = as.numeric(gsub("sd", "", id)),
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Tidy the base model data
  base_df <- base_sds %>%
    dplyr::mutate(
      peel = 0,
      item = gsub("\\..*", "", item)
    ) %>%
    dplyr::filter(item %in% quantities) %>%
    dplyr::select(item, Estimate, Std_Error = `Std. Error`, peel)
  
  # Combine and add years
  base_years <- rpt$rpt$years
  all_data <- dplyr::bind_rows(retro_df, base_df) %>%
    dplyr::group_by(peel, item) %>%
    dplyr::mutate(year = base_years[(peel + 1):(peel + dplyr::n())]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lci = Estimate - 1.96 * Std_Error,
      uci = Estimate + 1.96 * Std_Error,
      lci = ifelse(lci < 0, 0, lci),
      peel = factor(peel, levels = as.character(0:n_peels))
    )
  
  # Calculate Mohn's Rho and MAE ---
  # terminal_peels <- all_data %>%
  #   dplyr::filter(peel != 0) %>%
  #   dplyr::group_by(item, peel) %>%
  #   dplyr::summarise(year = max(year), peel_est = Estimate, .groups = 'drop')
  
  # terminal_base <- all_data %>%
  #   dplyr::filter(peel == 0) %>%
  #   dplyr::select(item, year, base_est = Estimate)
  
  # rho_metrics <- dplyr::left_join(terminal_peels, terminal_base, by = c("item", "year")) %>%
  #   dplyr::mutate(pdiff = (peel_est - base_est) / base_est) %>% # Using your formula
  #   tidyr::drop_na() %>%
  #   dplyr::group_by(item) %>%
  #   dplyr::summarise(
  #     rho = mean(pdiff, na.rm = TRUE),
  #     mae = median(abs(pdiff), na.rm = TRUE)
  #   )
  
  # --- 5. Generate and Save Plots ---
  message("Generating plots...")
  plot_list <- list()
  
  for (qty in quantities) {
    # metrics <- rho_metrics %>% dplyr::filter(item == qty)
    # l1 <- paste("Mohn's rho =", round(metrics$rho, 2))
    # l2 <- paste("MAE =", round(metrics$mae, 2))
    plot_title <- gsub("_", " ", qty)
    plot_title <- paste0(toupper(substring(plot_title, 1, 1)), substring(plot_title, 2))
    
    plot_data_qty <- all_data %>% dplyr::filter(item == qty)
    annot_x <- min(plot_data_qty$year) + 1 # Place annotation near the start
    annot_y <- max(plot_data_qty$uci, na.rm = TRUE) * 0.9 # Place it high up
    
    # Plot 1: Time series
    p1 <- ggplot2::ggplot(plot_data_qty, ggplot2::aes(x = year, y = Estimate, color = peel, group = peel, fill = peel)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      scico::scale_fill_scico_d("Peel", palette = 'roma', direction = -1) +
      # ggplot2::annotate(geom = 'text', x = annot_x, y = annot_y, label = paste(l1, l2, sep = "\n"), hjust = 0) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(y = paste(plot_title, "(t)"), x = "Year", title = paste("Prospective Pattern:", plot_title)) +
      ggplot2::expand_limits(y = 0)
    
    # Plot 2: Relative difference
    pdiff_data <- all_data %>%
      dplyr::filter(item == qty, peel != 0) %>%
      dplyr::left_join(
        all_data %>% dplyr::filter(item == qty, peel == 0) %>% dplyr::select(year, base_est = Estimate),
        by = "year"
      ) %>%
      dplyr::mutate(pdiff = (Estimate - base_est) / base_est) %>% # Standard relative diff for plotting
      tidyr::drop_na()
    
    p2 <- ggplot2::ggplot(pdiff_data, ggplot2::aes(x = year, y = pdiff, color = peel, group = peel)) +
      ggplot2::geom_line(show.legend = FALSE, linewidth = 0.8) +
      scico::scale_color_scico_d("Peel", palette = 'roma', direction = -1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(y = "Relative Difference from Base", x = "Year")
    
    # Combine and save
    combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1.5)) +
      patchwork::plot_layout(guides = "collect")
    
    plot_list[[qty]] <- combined_plot
    
    if (save_outputs) {
      file_name <- paste0(gsub(" ", "_", tolower(plot_title)), "_retro.png")
      ggplot2::ggsave(
        filename = file.path(retro_path, file_name),
        plot = combined_plot,
        width = 7, height = 7, units = "in", dpi = 300
      )
    }
  }
  
  message("Done.")
  # --- 6. Return Results ---
  return(list(
    prospect_data = all_data,
    # rho_metrics = rho_metrics,
    plots = plot_list
  ))
}
