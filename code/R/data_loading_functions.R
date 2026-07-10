library(R.matlab)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)

library(ggplot2)

load_channel_labels <- function(chlabels_mat) {
  
  extract_one_group <- function(x) {
    as.character(
      unlist(
        x,
        recursive = TRUE,
        use.names = FALSE
      )
    )
  }
  
  chlabels_by_group <- lapply(
    chlabels_mat,
    extract_one_group
  )
  
  chlabels_by_group
}

load_trials <- function(x, chlabels_by_group, value_name = "rpeak") {
  
  make_participant_long_df <- function(arr, chlabels, subID) {
    
    stopifnot(length(dim(arr)) == 3)
    stopifnot(dim(arr)[2] == length(chlabels))
    
    idx <- arrayInd(seq_along(arr), .dim = dim(arr))
    
    df <- data.frame(
      subID = subID,
      clean = idx[, 1],
      channel = chlabels[idx[, 2]],
      trial = idx[, 3],
      value = as.vector(arr)
    )
    
    names(df)[names(df) == "value"] <- value_name
    
    df
  }
  
  bind_rows(lapply(seq_along(x), function(i) {
    
    arr <- x[[i]][[1]]
    chs <- chlabels_by_group[[i]]
    
    make_participant_long_df(
      arr = arr,
      chlabels = chs,
      subID = i
    )
  })) %>%
    filter(!is.na(.data[[value_name]])) %>%
    mutate(
      subID = factor(subID),
      clean = as.integer(clean),
      channel = factor(channel),
      trial = as.integer(trial)
    )
}

load_avg<- function(x_avg, chlabels_by_group, value_name = "rpeak") {
  
  make_avg_participant_df <- function(mat, chlabels, subID) {
    
    stopifnot(is.matrix(mat))
    stopifnot(ncol(mat) == length(chlabels))
    
    colnames(mat) <- chlabels
    
    df <- as.data.frame(mat, check.names = FALSE)
    
    df$subID <- subID
    df$clean <- seq_len(nrow(mat))
    
    df %>%
      pivot_longer(
        cols = all_of(chlabels),
        names_to = "channel",
        values_to = value_name
      )
  }
  
  bind_rows(lapply(seq_along(x_avg), function(i) {
    
    mat <- x_avg[[i]][[1]]
    chs <- chlabels_by_group[[i]]
    
    make_avg_participant_df(
      mat = mat,
      chlabels = chs,
      subID = i
    )
  })) %>%
    filter(!is.na(.data[[value_name]])) %>%
    mutate(
      subID = factor(subID),
      clean = as.integer(clean),
      channel = factor(channel)
    )
}

load_avg_from_trials <- function(x_trials, chlabels_by_group, value_name = "rpeak") {
  
  make_avg_from_trials_participant_df <- function(arr, chlabels, subID) {
    
    stopifnot(length(dim(arr)) == 3)
    stopifnot(dim(arr)[2] == length(chlabels))
    
    # clean × channel × trial -> clean × channel
    mat_avg <- apply(arr, c(1, 2), mean, na.rm = TRUE)
    
    stopifnot(ncol(mat_avg) == length(chlabels))
    
    colnames(mat_avg) <- chlabels
    
    df <- as.data.frame(mat_avg, check.names = FALSE)
    
    df$subID <- subID
    df$clean <- seq_len(nrow(mat_avg))
    
    df %>%
      pivot_longer(
        cols = all_of(chlabels),
        names_to = "channel",
        values_to = value_name
      )
  }
  
  bind_rows(lapply(seq_along(x_trials), function(i) {
    
    arr <- x_trials[[i]][[1]]
    chs <- chlabels_by_group[[i]]
    
    make_avg_from_trials_participant_df(
      arr = arr,
      chlabels = chs,
      subID = i
    )
  })) %>%
    filter(!is.na(.data[[value_name]])) %>%
    mutate(
      subID = factor(subID),
      clean = as.integer(clean),
      channel = factor(channel)
    )
}

load_imcoh_long <- function(
    imcoh_band,
    imcoh_ecg_band = NULL,
    band_names = c("Delta", "Theta", "Alpha")
) {
  
  stopifnot(length(dim(imcoh_band)) == 3)
  
  n_sub <- dim(imcoh_band)[1]
  n_clean <- dim(imcoh_band)[2]
  n_band <- dim(imcoh_band)[3]
  
  stopifnot(n_band == length(band_names))
  
  idx_band <- arrayInd(seq_along(imcoh_band), .dim = dim(imcoh_band))
  
  df_band <- data.frame(
    subID = factor(idx_band[, 1]),
    clean = idx_band[, 2],
    band = factor(band_names[idx_band[, 3]], levels = band_names),
    imcoh = as.vector(imcoh_band)
  )
  
  df_ecg <- NULL
  
  if (!is.null(imcoh_ecg_band)) {
    
    stopifnot(length(dim(imcoh_ecg_band)) == 3)
    stopifnot(dim(imcoh_ecg_band)[1] == n_sub)
    stopifnot(dim(imcoh_ecg_band)[3] == n_band)
    stopifnot(dim(imcoh_ecg_band)[2] %in% c(1, n_clean))
    
    idx_ecg <- arrayInd(seq_along(imcoh_ecg_band), .dim = dim(imcoh_ecg_band))
    
    df_ecg <- data.frame(
      subID = factor(idx_ecg[, 1]),
      clean = idx_ecg[, 2],
      band = factor(band_names[idx_ecg[, 3]], levels = band_names),
      imcoh_ecg = as.vector(imcoh_ecg_band)
    )
  }
  
  list(
    band = df_band,
    ecg = df_ecg
  )
}

make_comb_lookup <- function(n_total_ecg_leads = 15, n_leads, ecg_lead_names = NULL) {
  
  comb_mat <- t(combn(seq_len(n_total_ecg_leads), n_leads))
  
  comb_lookup <- as_tibble(comb_mat)
  names(comb_lookup) <- paste0("ecg_lead_", seq_len(n_leads))
  
  comb_lookup <- comb_lookup %>%
    mutate(combination_id = row_number()) %>%
    relocate(combination_id)
  
  lead_cols <- paste0("ecg_lead_", seq_len(n_leads))
  
  comb_lookup$ecg_leads <- apply(
    comb_lookup[, lead_cols, drop = FALSE],
    1,
    paste,
    collapse = ","
  )
  
  if (!is.null(ecg_lead_names)) {
    stopifnot(length(ecg_lead_names) == n_total_ecg_leads)
    
    for (k in seq_len(n_leads)) {
      comb_lookup[[paste0("ecg_lead_", k, "_name")]] <-
        ecg_lead_names[comb_lookup[[paste0("ecg_lead_", k)]]]
    }
    
    name_cols <- paste0("ecg_lead_", seq_len(n_leads), "_name")
    
    comb_lookup$ecg_leads_named <- apply(
      comb_lookup[, name_cols, drop = FALSE],
      1,
      paste,
      collapse = ","
    )
  }
  
  comb_lookup
}

make_channel_lookup <- function(chlabels_by_group) {
  
  bind_rows(lapply(seq_along(chlabels_by_group), function(i) {
    tibble(
      participant = factor(i),
      channel_id = seq_along(chlabels_by_group[[i]]),
      channel = chlabels_by_group[[i]]
    )
  })) %>%
    mutate(
      channel_id = as.integer(channel_id),
      channel = factor(channel)
    )
}

load_comb_rpeak <- function(
    comb,
    n_leads,
    method = c("ICA", "CCA"),
    chlabels_by_group,
    n_total_ecg_leads = 15,
    ecg_lead_names = NULL,
    value_name = "rpeak",
    warn_padded = TRUE
) {
  
  method <- match.arg(method)
  field_name <- paste0("rpeak.", method, ".maxabs")
  
  if (!field_name %in% names(comb)) {
    stop("Field not found: ", field_name)
  }
  
  x <- comb[[field_name]]
  n_participants <- length(x)
  
  if (length(chlabels_by_group) != n_participants) {
    stop(
      "chlabels_by_group has length ", length(chlabels_by_group),
      ", but outer object has ", n_participants, " participant cells."
    )
  }
  
  expected_n_combs <- choose(n_total_ecg_leads, n_leads)
  expected_n_clean <- n_leads + 1
  
  comb_lookup <- make_comb_lookup(
    n_total_ecg_leads = n_total_ecg_leads,
    n_leads = n_leads,
    ecg_lead_names = ecg_lead_names
  )
  
  channel_lookup <- make_channel_lookup(chlabels_by_group)
  
  df <- bind_rows(lapply(seq_along(x), function(sub_i) {
    
    cell <- x[[sub_i]]
    
    arr <- if (is.list(cell) && length(cell) == 1) {
      cell[[1]]
    } else {
      cell
    }
    
    if (length(dim(arr)) != 4) {
      stop(
        "Participant ", sub_i,
        ": expected a 4D array, got dim: ",
        paste(dim(arr), collapse = " x ")
      )
    }
    
    arr_dim <- dim(arr)
    
    if (arr_dim[1] < n_leads) {
      stop(
        "Participant ", sub_i,
        ": first dimension is ", arr_dim[1],
        ", but n_leads = ", n_leads,
        ". Cannot extract row ", n_leads, "."
      )
    }
    
    if (arr_dim[2] < expected_n_combs) {
      stop(
        "Participant ", sub_i,
        ": expected at least ", expected_n_combs,
        " combinations for ", n_leads,
        " ECG leads, but array has only ", arr_dim[2], "."
      )
    }
    
    if (arr_dim[3] < expected_n_clean) {
      stop(
        "Participant ", sub_i,
        ": expected at least ", expected_n_clean,
        " cleaning levels for ", n_leads,
        " ECG leads, but array has only ", arr_dim[3], "."
      )
    }
    
    if (arr_dim[4] != length(chlabels_by_group[[sub_i]])) {
      stop(
        "Participant ", sub_i,
        ": channel dimension is ", arr_dim[4],
        ", but chlabels_by_group[[", sub_i, "]] has length ",
        length(chlabels_by_group[[sub_i]]), "."
      )
    }
    
    if (warn_padded && sub_i == 1) {
      if (arr_dim[2] > expected_n_combs) {
        warning(
          "For n_leads = ", n_leads, ", array has ", arr_dim[2],
          " combination slots, but only first ", expected_n_combs,
          " are valid. Ignoring padded/unused combinations."
        )
      }
      
      if (arr_dim[3] > expected_n_clean) {
        warning(
          "For n_leads = ", n_leads, ", array has ", arr_dim[3],
          " cleaning-level slots, but only first ", expected_n_clean,
          " are valid. Ignoring padded/unused cleaning levels."
        )
      }
    }
    
    valid_combinations <- seq_len(expected_n_combs)
    valid_clean <- seq_len(expected_n_clean)
    
    # Extract valid block only:
    # comb1 from combined file: arr[1, 1:15, 1:2, ]
    # comb2 from combined file: arr[2, 1:105, 1:3, ]
    # comb3 from combined file: arr[3, 1:455, 1:4, ]
    # comb4: arr[4, 1:1365, 1:5, ]
    # comb5: arr[5, 1:3003, 1:6, ]
    arr_selected <- arr[
      n_leads,
      valid_combinations,
      valid_clean,
      ,
      drop = FALSE
    ]
    
    dim(arr_selected) <- c(
      expected_n_combs,
      expected_n_clean,
      arr_dim[4]
    )
    
    out <- expand.grid(
      combination_id = seq_len(expected_n_combs),
      clean = seq_len(expected_n_clean),
      channel_id = seq_len(arr_dim[4]),
      KEEP.OUT.ATTRS = FALSE
    ) %>%
      as_tibble()
    
    out[[value_name]] <- as.vector(arr_selected)
    
    out %>%
      mutate(
        participant = factor(sub_i),
        n_ecg_leads = n_leads,
        method = method,
        combination_id = as.integer(combination_id),
        clean = as.integer(clean),
        channel_id = as.integer(channel_id)
      )
  }))
  
  df %>%
    left_join(channel_lookup, by = c("participant", "channel_id")) %>%
    left_join(comb_lookup, by = "combination_id") %>%
    filter(!is.na(.data[[value_name]])) %>%
    relocate(
      n_ecg_leads,
      method,
      participant,
      clean,
      channel,
      channel_id,
      combination_id
    )
}

load_comb_rpeak_methods <- function(
    comb,
    n_leads,
    methods = c("ICA", "CCA"),
    chlabels_by_group,
    n_total_ecg_leads = 15,
    ecg_lead_names = NULL,
    value_name = "rpeak",
    warn_padded = TRUE
) {
  
  methods <- match.arg(methods, choices = c("ICA", "CCA"), several.ok = TRUE)
  
  purrr::map_dfr(methods, function(m) {
    load_comb_rpeak(
      comb = comb,
      n_leads = n_leads,
      method = m,
      chlabels_by_group = chlabels_by_group,
      n_total_ecg_leads = n_total_ecg_leads,
      ecg_lead_names = ecg_lead_names,
      value_name = value_name,
      warn_padded = warn_padded
    )
  })
}


format_pid <- function(local_id, offset = 0) {
  sprintf("P%02d", as.integer(local_id) + offset)
}

load_chlabels_dataset <- function(chlabel_file) {
  ch_data <- readMat(file.path(datapath, chlabel_file))
  load_channel_labels(ch_data$chlabels)
}

make_channel_count_df <- function(chlabels_by_group, dataset, offset = 0) {
  tibble(
    dataset = dataset,
    participant_local = seq_along(chlabels_by_group),
    participant = format_pid(participant_local, offset),
    n_channels = map_int(chlabels_by_group, length)
  )
}

load_full15_rpeak_dataset <- function(
    data_file,
    chlabels_by_group,
    dataset,
    offset = 0,
    methods = "CCA",
    measure = c("rpeak", "alpha_peak"),
    value_name = NULL
) {
  
  measure <- match.arg(measure)
  
  if (is.null(value_name)) {
    value_name <- measure
  }
  
  mat_data <- readMat(file.path(datapath, data_file))
  
  out <- map_dfr(methods, function(m) {
    
    field_name <- switch(
      measure,
      rpeak      = paste0("rpeak.", m, ".maxabs"),
      alpha_peak = paste0("alpha.peak.", m)
    )
    
    if (!field_name %in% names(mat_data)) {
      stop("Field not found in ", data_file, ": ", field_name)
    }
    
    load_avg(
      mat_data[[field_name]],
      chlabels_by_group,
      value_name = value_name
    ) %>%
      mutate(method = m)
  })
  
  mat_data <- NULL
  gc()
  
  out %>%
    mutate(
      dataset = dataset,
      participant_local = as.integer(as.character(subID)),
      participant = format_pid(participant_local, offset),
      subID = participant,
      channel = as.character(channel),
      method = factor(method, levels = c("ICA", "CCA")),
      clean = as.integer(clean)
    ) %>%
    relocate(
      dataset,
      participant,
      participant_local,
      subID,
      method,
      clean,
      channel
    )
}

load_comb_dataset <- function(
    file_plan,
    chlabels_by_group,
    dataset,
    offset = 0,
    methods = "CCA"
) {
  
  stopifnot(all(c("n_ecg_leads", "file") %in% names(file_plan)))
  
  out <- split(file_plan, file_plan$file) %>%
    imap_dfr(function(this_plan, file_name) {
      
      comb <- readMat(file.path(datapath, file_name))
      
      this_out <- map_dfr(this_plan$n_ecg_leads, function(k) {
        load_comb_rpeak_methods(
          comb = comb,
          n_leads = k,
          methods = methods,
          chlabels_by_group = chlabels_by_group
        )
      })
      
      comb <- NULL
      gc()
      
      this_out
    })
  
  out %>%
    mutate(
      dataset = dataset,
      participant_local = as.integer(as.character(participant)),
      participant = format_pid(participant_local, offset),
      participant = factor(participant),
      channel = as.character(channel),
      method = factor(method, levels = c("ICA", "CCA")),
      clean = as.integer(clean),
      n_ecg_leads = as.integer(n_ecg_leads),
      combination_id = as.integer(combination_id),
      ecg_leads = as.character(ecg_leads)
    ) %>%
    relocate(
      dataset,
      participant,
      participant_local,
      n_ecg_leads,
      method,
      clean,
      channel,
      channel_id,
      combination_id,
      ecg_leads
    )
}