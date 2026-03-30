# Estimates effects of Q and A on choice using an approximate data analysis method
est_bq_ba_s2 <- function(data){
  # Get Q values
  data_list <- by(data, data$sub_index, add_assoc, cue_cols = c("fA_ix","fB_ix"), out_cols = c("out_a","out_b"), 
                  num_assoc = 8, assoc_name = "q", lrn_rate = "alpha", for_rate = "forget")
  data <- do.call(rbind,data_list)
  
  for(r in 1:nrow(data)){
    # Get chosen Q value
    data$q_ch[r] <- data[r,paste0("q_",data$chosen_frac[r])]
  }
  
  # Add estimates of the affective impacts of outcomes
  data_list <- list()
  for(s in 1:max(data$sub_index)){
    sub_data <- filter(data,sub_index == s)
    
    val_fit <- lm(feed_rate_z ~ chosen_out + unchosen_out + q_ch + prat + block + trial_nl, sub_data)
    
    sub_data <- sub_data %>% mutate(mod_val = val_fit$coefficients[1] + val_fit$coefficients[2]*chosen_out +
                                    val_fit$coefficients[3]*unchosen_out + val_fit$coefficients[4]*q_ch)
    data_list[[s]] <- sub_data
  }
  data <- do.call(rbind,data_list)
  
  # Get remaining associations used to predict choice
  data <- data %>%
    mutate(fA_chosen=ifelse(choice==1,1,0)) %>%
    mutate(fB_chosen=ifelse(choice==2,1,0))
  
  data <- data %>%
    split(data$sub_index) %>%
    lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "mod_val",
           num_assoc = 8, assoc_name = "a", lrn_rate = "alpha", for_rate = "forget") %>%
    lapply(add_assoc, cue_cols = c("fA_ix","fB_ix"), out_cols = c("fA_chosen","fB_chosen"),
           num_assoc = 8, assoc_name = "c", lrn_rate = "tau") %>%
    bind_rows()
  
  # Get choice predictor differences
  for(r in 1:nrow(data)){
    data[r,"c_A"] <- data[r,paste0("c_",data[r,"fA_ix"])]
    data[r,"c_B"] <- data[r,paste0("c_",data[r,"fB_ix"])]
    data[r,"q_A"] <- data[r,paste0("q_",data[r,"fA_ix"])]
    data[r,"q_B"] <- data[r,paste0("q_",data[r,"fB_ix"])]
    data[r,"a_A"] <- data[r,paste0("a_",data[r,"fA_ix"])]
    data[r,"a_B"] <- data[r,paste0("a_",data[r,"fB_ix"])]
  }
  data <- data %>%
    mutate(c_diff = c_A - c_B) %>%
    mutate(q_diff = q_A - q_B) %>%
    mutate(a_diff = a_A - a_B)
  
  # Estimate standardized effects on choice
  choice_fit <- glm(fA_chosen ~ scale(q_diff) + scale(a_diff) + scale(c_diff),
                    data, family = "binomial")
  
  list("bQ"=choice_fit$coefficients[2],"bA"=choice_fit$coefficients[3])
}

#Adds columns to a trials df corresponding to the fractals at play in the current block. 
#df: a data frame with a specific block's trials
create_pair_cols_s2 <- function(df){
  unique_pairs <- df[!duplicated(df$fA_ix), ] #get the first 2 rows with unique pairs
  #add columns to the df with these unique pairs
  df_ret <- df %>% mutate(fA1 = unique_pairs$fA_ix[1],fB1 = unique_pairs$fB_ix[1],
                          fA2 = unique_pairs$fA_ix[2],fB2 = unique_pairs$fB_ix[2])
  return(df_ret)
}

#A wrapper functino that runs create_pair_cols on each of the blocks for a particular subject and stitches them together.
create_pair_cols_sub_s2 <- function(df){
  new_df_list <- by(df,df$block,create_pair_cols_s2)
  df_ret <- do.call(rbind,new_df_list)
  return(df_ret)
}

#Returns a list of data for input to stan, given trial-level data
#trials: trial-level data
#n_t: number of trials; must be set manually
stan_data_s22fu <- function(trials,n_t){

  #NEED TO CHECK THIS, ESP TOWARD THE END
  library(gdata)
  
  n_s <- length(unique(trials$id)) #get number of subjects
  
  frac_vec <- as.vector(as.matrix(trials[c("fA_ix","fB_ix")])) #get a big vector of fractals...
  n_f <- length(unique(frac_vec)) #...so we can get number of fractals
  
  fA <- sub_by_trial_matrix(trials,"fA_ix") #get a matrix with subjects on the rows and trials in the columns, containing indices for fractal A                                                    on that trial
  fB <- sub_by_trial_matrix(trials,"fB_ix") #ditto fractal B
  
  chosen_frac <- sub_by_trial_matrix(trials,"chosen_frac")
  unchosen_frac <- sub_by_trial_matrix(trials,"unchosen_frac")
  
  choice <- sub_by_trial_matrix(trials,"choice_numeric") #diddo whether fractal A (1) or B(2) was chosen
  
  out_a <- sub_by_trial_vec_list(trials,"out_a") #get a list of vectors, one per subject, each containing outcome a for each trial
  out_b <- sub_by_trial_vec_list(trials,"out_b") #diddo outcome b
  
  chosen_out <- sub_by_trial_vec_list(trials,"chosen_out") #get a list of vectors, one per subject, each containing outcome a for each trial
  unchosen_out <- sub_by_trial_vec_list(trials,"unchosen_out") #diddo outcome b
  
  number_blocks <- max(unique(trials$block))
  n_tb <- n_t/number_blocks #number of trials per block
  t_in_b <- sub_by_trial_vec_list(trials,"trial_nl") #trial-in-block; useful for V calculation
  
  decrate_first <- trials[!duplicated(trials$sub_index),"decrate_first"] #get a vector with one item per sub_index, 
                                                                         #indicating whether they completed the decision ratings in the first block
  
  
  dp_num <- sub_by_trial_matrix(trials,"dec_probe_number")
  fp_num <- sub_by_trial_matrix(trials,"feed_probe_number")
  
  #number of trials with each type of affect probe
  n_dp <- max(trials$dec_probe_number)
  n_fp <- max(trials$feed_probe_number)
  
  #get a vector of decision ratings, z-scored, one for each trial on which there was a decision rating
  dp_trials <- filter(trials,dec_probe_number != 0)
  dec_rate <- dp_trials$dec_rate_z
  
  #ditto feedback
  fp_trials <- filter(trials,feed_probe_number != 0)
  feed_rate <- fp_trials$feed_rate_z
  
  block_c <- sub_by_trial_vec_list(trials,"block_cent") #block predictor for each subject/trial
  tinb_c <- sub_by_trial_vec_list(trials,"trial_nl_cent") #trial predictor for each subject/trial
  
  prev_rate <- sub_by_trial_vec_list(trials,"prev_rate") #previous valence rating
  
  #delete the below if not doing ARL analyses
  trials <- trials %>% mutate(feed_rate_0 = ifelse(is.na(feed_rate_z),0,feed_rate_z)) #create column that is feed_rate_z with NA replaced by 0
  feed_rate_0 <- sub_by_trial_vec_list(trials,"feed_rate_0")

  pairs <- array(dim = c(n_s,n_t,2,2)) #declare the array of fractal pair indices
  for(s in 1:n_s){
    df <- trials %>% filter(sub_index == s) #grab the trials for this subject
    for(t in 1:n_t){
      #fill out the array for this subject and trial
      pairs[s,t,1,1] <- df$fA1[t]
      pairs[s,t,1,2] <- df$fB1[t]
      pairs[s,t,2,1] <- df$fA2[t]
      pairs[s,t,2,2] <- df$fB2[t]
    }
  }
  
  data <- list(
    n_t = n_t,
    n_s = n_s,
    n_f = n_f,
    fA = fA,
    fB = fB,
    chosen_frac = chosen_frac,
    unchosen_frac = unchosen_frac,
    choice = choice,
    out_a = out_a,
    out_b = out_b,
    chosen_out = chosen_out,
    unchosen_out = unchosen_out,
    n_tb = n_tb,
    t_in_b = t_in_b,
    decrate_first = decrate_first,
    dp_num = dp_num,
    fp_num = fp_num,
    n_dp = n_dp,
    n_fp = n_fp,
    dec_rate = dec_rate,
    feed_rate = feed_rate,
    block_c = block_c,
    tinb_c = tinb_c,
    prev_rate = prev_rate,
    feed_rate_0 = feed_rate_0,
    pairs = pairs
  )
  return(data)
}

#adds a column for previous valence rating to the trial-level dataset, returning that dataset
add_prevrate <- function(df,rat_col_name,pr_col_name="prev_rate"){
  df[pr_col_name] <- 0
  #for each trial in the block past the first...
  for(t in 2:nrow(df)){
    #if there was a rating on the previous trial...
    if(!is.na(df[t-1,rat_col_name])){
      df[t,pr_col_name] <- df[t-1,rat_col_name] #that's the prev_rate for this trial...
    } else{
      #otherwise, the prev_rate is the previous rating from the previous trial (the most recent rating)
      df[t,pr_col_name] <- df[t-1,pr_col_name]
    } 
  }
  return(df)
}

# Study 2 vector-length recovery simulation helpers

# This block implements the Study 2 vector-length recovery workflow.
# High-level flow:
# 1) Draw subject-level non-RL parameters from population distributions (posterior-median anchored).
# 2) Draw composite-variable effects for one generating condition (reward / PE / CC / blend).
# 3) Convert component draws to the 7 raw coefficients used in the affect model.
# 4) Simulate trial-level predictors and affect ratings.
# 5) Recover RL effects with OLS + subject fixed effects.
# 6) Convert recovered effects to SV/PE/CC/residual vector lengths.
# 7) Return tidy outputs and plots.

# Returns the median estimate for a named parameter from fsml$sum.
# sum_df: posterior summary table from read_fsml(... )$sum
# var_name: parameter name in the summary table
get_median <- function(sum_df, var_name) {
  row <- sum_df[sum_df$variable == var_name, , drop = FALSE]
  if (nrow(row) == 0) {
    stop(paste0("Could not find parameter in summary: ", var_name))
  }
  as.numeric(row$median[1])
}

# Draws one set of subject-level nuisance/choice parameters for one simulated dataset.
# Subject-level values are sampled from normal population distributions whose centers
# are set by posterior medians of group-level parameters.
s2_draw_subject_params <- function(sum_df, sub_df) {
  n_s <- nrow(sub_df)

  # Choice/RL parameter hyperparameters
  alpha_mu <- get_median(sum_df, "alpha_mu")
  alpha_sigma <- get_median(sum_df, "alpha_sigma")
  forget_mu <- get_median(sum_df, "forget_mu")
  forget_sigma <- get_median(sum_df, "forget_sigma")
  beta_mu <- get_median(sum_df, "beta_mu")
  beta_sigma <- get_median(sum_df, "beta_sigma")
  phi_mu <- get_median(sum_df, "phi_mu")
  phi_sigma <- get_median(sum_df, "phi_sigma")
  tau_mu <- get_median(sum_df, "tau_mu")
  tau_sigma <- get_median(sum_df, "tau_sigma")

  # Affect nuisance parameter hyperparameters
  base_mu <- get_median(sum_df, "base_mu")
  base_sigma <- get_median(sum_df, "base_sigma")
  blk_num_mu <- get_median(sum_df, "blk_num_mu")
  dec_mu <- get_median(sum_df, "dec_mu")
  bnum_dec_sigma <- get_median(sum_df, "bnum_dec_sigma")
  tr_mu <- get_median(sum_df, "tr_mu")
  tr_sigma <- get_median(sum_df, "tr_sigma")
  pr_mu <- get_median(sum_df, "pr_mu")
  pr_sigma <- get_median(sum_df, "pr_sigma")

  # Draw one row per subject.
  # Constrained parameters are sampled on latent scale then transformed.
  data.frame(
    sub_index = sub_df$sub_index,
    decrate_first = sub_df$decrate_first,
    alpha = plogis(rnorm(n_s, mean = alpha_mu, sd = alpha_sigma)),
    forget = plogis(rnorm(n_s, mean = forget_mu, sd = forget_sigma)),
    beta = rnorm(n_s, mean = beta_mu, sd = beta_sigma),
    phi = rnorm(n_s, mean = phi_mu, sd = phi_sigma),
    tau = plogis(rnorm(n_s, mean = tau_mu, sd = tau_sigma)),
    base = rnorm(n_s, mean = base_mu, sd = base_sigma),
    blk_z = rnorm(n_s, mean = 0, sd = 1),
    tr = rnorm(n_s, mean = tr_mu, sd = tr_sigma),
    pr = rnorm(n_s, mean = pr_mu, sd = pr_sigma),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      # Keep the same block-effect construction used in combined_shrink:
      # subjects who did decision ratings first vs second get opposite signs on dec_mu.
      blk = ifelse(
        decrate_first == 0,
        blk_num_mu + dec_mu + bnum_dec_sigma * blk_z,
        blk_num_mu - dec_mu + bnum_dec_sigma * blk_z
      )
    )
}

# Returns draws from a normal distribution truncated to [lower, upper].
draw_trunc_norm <- function(n, mean, sd, lower, upper) {
  p_low <- stats::pnorm(lower, mean = mean, sd = sd)
  p_high <- stats::pnorm(upper, mean = mean, sd = sd)
  u <- stats::runif(n, min = p_low, max = p_high)
  stats::qnorm(u, mean = mean, sd = sd)
}

# Returns the 7 prediction direction vectors used in decomposition/recovery.
# Rows align with RL coefficients in order:
# V_block, Q_ch_dec, Q_unch_dec, V_trial, Q_ch_feed, r_ch, r_unch.
prediction_direction_matrix <- function() {
  pe_ch <- c(-0.5, 0.5, 0, 0, 0, 0, 0)
  r_ch <- c(0, 1, 0, 0, 0, 0, 0)
  cc_ch <- c(0, 0.5, -0.5, 0, 0, 0, 0)
  pe_out_v <- c(0, 0, 0, -0.5, 0, 0.5, 0)
  pe_out_q <- c(0, 0, 0, 0, -0.5, 0.5, 0)
  r_out <- c(0, 0, 0, 0, 0, 1, 0)
  cc_out <- c(0, 0, 0, 0, 0, 0.5, -0.5)
  cbind(pe_ch, r_ch, cc_ch, pe_out_v, pe_out_q, r_out, cc_out)
}

# Draws subject-level true prediction-vector lengths for one generating condition.
# Active prediction vectors are drawn from N(0.5, 0.25), truncated to [0.05, 0.95].
# Inactive prediction vectors are fixed at 0.
draw_prediction_lengths <- function(n_s, condition) {
  pred_names <- c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out")
  out <- as.data.frame(matrix(0, nrow = n_s, ncol = length(pred_names)))
  names(out) <- pred_names

  active <- switch(
    condition,
    reward = c("r_choice", "r_out"),
    pe = c("pe_choice", "pe_out_q", "pe_out_v"),
    cc = c("cc_choice", "cc_out"),
    blend = pred_names,
    stop(paste0("Unknown condition: ", condition))
  )

  for (nm in active) {
    out[[nm]] <- draw_trunc_norm(
      n = n_s, mean = 0.5, sd = 0.25, lower = 0.05, upper = 0.95
    )
  }

  out
}

# Converts prediction-vector lengths to subject-level RL-variable coefficients.
# Returned columns correspond to the 7 RL predictors used in simulation:
# V_block, Q_ch_dec, Q_unch_dec, V_trial, Q_ch_feed, r_ch, r_unch.
prediction_lengths_to_rl_effects <- function(pred_length_df) {
  pred_order <- c("pe_choice", "r_choice", "cc_choice", "pe_out_v", "pe_out_q", "r_out", "cc_out")
  dvec_mat <- prediction_direction_matrix()
  length_mat <- as.matrix(pred_length_df[, c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out"), drop = FALSE])
  length_mat <- length_mat[, pred_order, drop = FALSE]
  eff_mat <- length_mat %*% t(dvec_mat)

  out <- as.data.frame(eff_mat, stringsAsFactors = FALSE)
  names(out) <- c("b_v_block", "b_q_ch_dec", "b_q_unch_dec", "b_v_trial", "b_q_ch_feed", "b_r_ch", "b_r_unch")
  out
}

# Simulates subject-level latent states and RL predictors used by the affect model.
# Choices are sampled from the model's softmax policy on each trial.
# This function mirrors the update structure in combined_shrink:
# - Q and C states evolve trial-by-trial
# - V resets at block starts and updates from chosen outcomes
# - choice probability uses beta*Q + phi*C
# The function outputs trial-level predictor columns that later feed simulated ratings.
s2_simulate_predictors_one_subject <- function(
  df_sub,
  alpha,
  forget,
  beta,
  phi,
  tau
) {
  # Ensure deterministic trial order within subject.
  df_sub <- df_sub[order(df_sub$overall_trial_nl), , drop = FALSE]
  n_t <- nrow(df_sub)
  n_f <- max(c(df_sub$fA_ix, df_sub$fB_ix), na.rm = TRUE)

  # Initialize latent states at 0 at start of subject.
  Q <- rep(0, n_f)
  C <- rep(0, n_f)
  V <- 0

  # Pre-allocate output predictor columns.
  ret <- df_sub
  ret$V_block_sim <- NA_real_
  ret$Q_ch_dec_sim <- NA_real_
  ret$Q_unch_dec_sim <- NA_real_
  ret$V_trial_sim <- NA_real_
  ret$Q_ch_feed_sim <- NA_real_
  ret$r_ch_sim <- NA_real_
  ret$r_unch_sim <- NA_real_
  ret$choice_numeric_sim <- NA_integer_
  ret$chosen_frac_sim <- NA_integer_
  ret$unchosen_frac_sim <- NA_integer_
  ret$chosen_out_sim <- NA_real_
  ret$unchosen_out_sim <- NA_real_

  for (t in seq_len(n_t)) {
    if (ret$trial_nl[t] == 1) {
      V <- 0
    }

    fA <- ret$fA_ix[t]
    fB <- ret$fB_ix[t]
    outA <- ret$out_a[t]
    outB <- ret$out_b[t]

    # Stable softmax for P(choose A), then sample a simulated choice.
    etaA <- beta * Q[fA] + phi * C[fA]
    etaB <- beta * Q[fB] + phi * C[fB]
    eta_max <- max(etaA, etaB)
    pA <- exp(etaA - eta_max) / (exp(etaA - eta_max) + exp(etaB - eta_max))
    chose_A <- stats::rbinom(1, size = 1, prob = pA)
    choice_num_sim <- ifelse(chose_A == 1, 1L, 2L)
    ch_sim <- ifelse(choice_num_sim == 1L, fA, fB)
    unch_sim <- ifelse(choice_num_sim == 1L, fB, fA)
    ch_out_sim <- ifelse(choice_num_sim == 1L, outA, outB)
    unch_out_sim <- ifelse(choice_num_sim == 1L, outB, outA)
    p_ch <- ifelse(choice_num_sim == 1L, pA, 1 - pA)

    ret$V_block_sim[t] <- V
    ret$Q_ch_dec_sim[t] <- Q[ch_sim]
    ret$Q_unch_dec_sim[t] <- Q[unch_sim]
    ret$V_trial_sim[t] <- ch_out_sim * p_ch + unch_out_sim * (1 - p_ch)
    ret$Q_ch_feed_sim[t] <- Q[ch_sim]
    ret$r_ch_sim[t] <- ch_out_sim
    ret$r_unch_sim[t] <- unch_out_sim
    ret$choice_numeric_sim[t] <- choice_num_sim
    ret$chosen_frac_sim[t] <- ch_sim
    ret$unchosen_frac_sim[t] <- unch_sim
    ret$chosen_out_sim[t] <- ch_out_sim
    ret$unchosen_out_sim[t] <- unch_out_sim

    # Q/C updates follow combined_shrink structure.
    old_Q <- Q
    old_C <- C
    Q <- (1 - forget) * Q
    Q[fA] <- old_Q[fA] + alpha * (outA - old_Q[fA])
    Q[fB] <- old_Q[fB] + alpha * (outB - old_Q[fB])

    # Update C based on simulated choice.
    if (choice_num_sim == 1L) {
      choice_a <- 1
      choice_b <- 0
    } else {
      choice_a <- 0
      choice_b <- 1
    }

    # Choice-autocorrelation state update for the two shown cues.
    C <- old_C
    C[fA] <- old_C[fA] + tau * (choice_a - old_C[fA])
    C[fB] <- old_C[fB] + tau * (choice_b - old_C[fB])

    # V update from simulated chosen outcome.
    V <- V + alpha * (ch_out_sim - V)
  }

  ret
}

# Simulates one full trial-level dataset (all subjects) for a given parameter draw set.
# For each subject:
# 1) simulate latent predictors
# 2) generate observed decision/outcome ratings with nuisance + RL components + residual noise
# 3) carry forward previous simulated rating for prev_rate_sim
# trials must contain scheduled probe columns:
# dec_probe_number_sched and feed_probe_number_sched
s2_simulate_affect_dataset <- function(
  trials,
  subj_params,
  rl_effects,
  d_resid,
  f_resid
) {
  # Pre-allocate one list element per subject for speed and clarity.
  out_list <- vector("list", length = nrow(subj_params))

  for (i in seq_len(nrow(subj_params))) {
    s <- subj_params$sub_index[i]
    df_sub <- trials[trials$sub_index == s, , drop = FALSE]
    b <- rl_effects[i, , drop = FALSE]

    pred_sub <- s2_simulate_predictors_one_subject(
      df_sub = df_sub,
      alpha = subj_params$alpha[i],
      forget = subj_params$forget[i],
      beta = subj_params$beta[i],
      phi = subj_params$phi[i],
      tau = subj_params$tau[i]
    )

    prev_sim <- 0

    pred_sub$prev_rate_sim <- NA_real_
    pred_sub$dec_rate_sim <- NA_real_
    pred_sub$feed_rate_sim <- NA_real_
    pred_sub$val_rate_sim <- NA_real_

    for (t in seq_len(nrow(pred_sub))) {
      # prev_rate_sim is the most recent simulated rating available before trial t.
      pred_sub$prev_rate_sim[t] <- prev_sim

      # Nuisance terms shared by decision and feedback probes.
      nuis <- subj_params$base[i] +
        subj_params$blk[i] * pred_sub$block_cent[t] +
        subj_params$tr[i] * pred_sub$trial_nl_cent[t] +
        subj_params$pr[i] * prev_sim

      # RL components from subject-specific coefficients implied by prediction-vector lengths.
      rl_dec <- b$b_v_block * pred_sub$V_block_sim[t] +
        b$b_q_ch_dec * pred_sub$Q_ch_dec_sim[t] +
        b$b_q_unch_dec * pred_sub$Q_unch_dec_sim[t]

      rl_feed <- b$b_v_trial * pred_sub$V_trial_sim[t] +
        b$b_q_ch_feed * pred_sub$Q_ch_feed_sim[t] +
        b$b_r_ch * pred_sub$r_ch_sim[t] +
        b$b_r_unch * pred_sub$r_unch_sim[t]

      # Generate whichever rating type exists on this trial.
      if (pred_sub$dec_probe_number_sched[t] != 0) {
        y <- nuis + rl_dec + rnorm(1, mean = 0, sd = d_resid)
        pred_sub$dec_rate_sim[t] <- y
        pred_sub$val_rate_sim[t] <- y
        prev_sim <- y
      } else if (pred_sub$feed_probe_number_sched[t] != 0) {
        y <- nuis + rl_feed + rnorm(1, mean = 0, sd = f_resid)
        pred_sub$feed_rate_sim[t] <- y
        pred_sub$val_rate_sim[t] <- y
        prev_sim <- y
      } else {
        pred_sub$val_rate_sim[t] <- prev_sim
      }
    }

    # Store subject dataset.
    out_list[[i]] <- pred_sub
  }

  # Return one long trial-level simulated dataset.
  dplyr::bind_rows(out_list)
}

# Fits recovery models to simulated ratings using OLS with subject fixed effects.
# This avoids heavier mixed-model fitting while still controlling baseline subject differences.
fit_recovery_models <- function(sim_df) {
  # Separate by probe type to match study design and original model structure.
  dec_df <- sim_df[!is.na(sim_df$dec_rate_sim), , drop = FALSE]
  feed_df <- sim_df[!is.na(sim_df$feed_rate_sim), , drop = FALSE]

  dec_fit <- stats::lm(
    dec_rate_sim ~ V_block_sim + Q_ch_dec_sim + Q_unch_dec_sim +
      block_cent + trial_nl_cent + prev_rate_sim + factor(sub_index),
    data = dec_df
  )

  feed_fit <- stats::lm(
    feed_rate_sim ~ V_trial_sim + Q_ch_feed_sim + r_ch_sim + r_unch_sim +
      block_cent + trial_nl_cent + prev_rate_sim + factor(sub_index),
    data = feed_df
  )

  list(dec_fit = dec_fit, feed_fit = feed_fit)
}

# Extracts the 7 RL coefficients in the exact order needed by prediction_lengths_from_effects.
extract_effect_vector <- function(fit_list) {
  dcoef <- stats::coef(fit_list$dec_fit)
  fcoef <- stats::coef(fit_list$feed_fit)

  c(
    dcoef["V_block_sim"],
    dcoef["Q_ch_dec_sim"],
    dcoef["Q_unch_dec_sim"],
    fcoef["V_trial_sim"],
    fcoef["Q_ch_feed_sim"],
    fcoef["r_ch_sim"],
    fcoef["r_unch_sim"]
  )
}

# Decomposes one 7D RL effect vector into raw prediction-vector lengths plus residual.
# Unlike get_vec_lens(), this works on one effect vector rather than posterior draws.
prediction_lengths_from_effects <- function(eff_vec) {
  dvec_mat <- prediction_direction_matrix()

  ws <- vec_optim(target = eff_vec, dvec = dvec_mat, init_pars = rep(0, 7))
  resid_raw <- sum(abs(eff_vec)) - sum(ws)
  raw_lengths <- c(ws, resid_raw)
  names(raw_lengths) <- c("pe_ch", "r_ch", "cc_ch", "pe_out_v", "pe_out_q", "r_out", "cc_out", "resid")

  # Raw prediction-vector lengths
  lens <- c(
    r_choice = as.numeric(raw_lengths["r_ch"]),
    pe_choice = as.numeric(raw_lengths["pe_ch"]),
    cc_choice = as.numeric(raw_lengths["cc_ch"]),
    r_out = as.numeric(raw_lengths["r_out"]),
    pe_out_q = as.numeric(raw_lengths["pe_out_q"]),
    pe_out_v = as.numeric(raw_lengths["pe_out_v"]),
    cc_out = as.numeric(raw_lengths["cc_out"]),
    resid = as.numeric(raw_lengths["resid"])
  )
  lens
}

# Returns nominal reference lengths used as target markers in prediction-recovery plots.
vector_recovery_reference_lengths <- function() {
  dplyr::bind_rows(
    data.frame(
      condition = "reward",
      prediction = c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out", "resid"),
      target = c(0.5, 0, 0, 0.5, 0, 0, 0, 0)
    ),
    data.frame(
      condition = "pe",
      prediction = c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out", "resid"),
      target = c(0, 0.5, 0, 0, 0.5, 0.5, 0, 0)
    ),
    data.frame(
      condition = "cc",
      prediction = c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out", "resid"),
      target = c(0, 0, 0.5, 0, 0, 0, 0.5, 0)
    ),
    data.frame(
      condition = "blend",
      prediction = c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out", "resid"),
      target = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0)
    )
  )
}

# Dot-cloud plot for recovered vector lengths.
# Uses vertical jitter (only) to reduce overplotting while preserving exact x values.
# dot_alpha/dot_size/dot_jitter_height control point rendering in dot mode.
plot_vector_recovery_dot <- function(
  length_df,
  ref_df = vector_recovery_reference_lengths(),
  x_limits = c(0, 1),
  hide_axis_text = TRUE,
  dot_alpha = 0.06,
  dot_size = 1.50,
  dot_jitter_height = 0.12
) {
  prediction_key <- data.frame(
    prediction = c("resid", "cc_out", "pe_out_v", "pe_out_q", "r_out", "cc_choice", "pe_choice", "r_choice"),
    y_id = c(1, 2, 3, 4, 5, 6, 7, 8),
    y_lab = c("Residual", "CC_out", "PE_out,V", "PE_out,Q", "R_out", "CC_choice", "PE_choice", "R_choice"),
    stringsAsFactors = FALSE
  )
  condition_levels <- c("reward", "pe", "cc", "blend")
  condition_labels <- c("Reward", "PE", "CC", "Blend")

  length_plot_df <- length_df |>
    dplyr::left_join(prediction_key, by = "prediction") |>
    dplyr::mutate(
      prediction = factor(prediction, levels = c("r_choice", "pe_choice", "cc_choice", "r_out", "pe_out_q", "pe_out_v", "cc_out", "resid")),
      condition = factor(condition, levels = condition_levels, labels = condition_labels)
    )

  ref_plot_df <- ref_df |>
    dplyr::left_join(prediction_key, by = "prediction") |>
    dplyr::mutate(
      condition = factor(condition, levels = condition_levels, labels = condition_labels)
    )

  # Median x-position per condition x prediction row (drawn slightly below cloud).
  median_df <- length_plot_df |>
    dplyr::group_by(condition, prediction, y_id) |>
    dplyr::summarise(length_med = stats::median(length), .groups = "drop")

  x_text <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_text()
  y_text <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_text()
  x_ticks <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_line(color = "black", linewidth = 0.47)
  y_ticks <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_blank()

  ggplot2::ggplot(
    length_plot_df,
    ggplot2::aes(x = length, y = y_id)
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.47) +
    # Draw an explicit x-axis baseline in every facet (including top row).
    ggplot2::geom_hline(yintercept = 0.5, color = "black", linewidth = 0.47) +
    ggplot2::geom_point(
      ggplot2::aes(color = prediction),
      alpha = dot_alpha,
      size = dot_size,
      shape = 16,
      position = ggplot2::position_jitter(height = dot_jitter_height, width = 0)
    ) +
    ggplot2::geom_point(
      data = median_df,
      ggplot2::aes(x = length_med, y = y_id - 0.19),
      inherit.aes = FALSE,
      shape = 23,
      size = 1.7,
      stroke = 0.45,
      alpha = 1,
      fill = "white",
      color = "black"
    ) +
    # Draw target ticks above points.
    ggplot2::geom_segment(
      data = ref_plot_df,
      ggplot2::aes(x = target, xend = target, y = y_id - 0.24, yend = y_id + 0.24),
      inherit.aes = FALSE,
      color = "#22A652",
      linewidth = 0.85,
      lineend = "round",
      alpha = 1
    ) +
    ggplot2::facet_wrap(~condition, ncol = 2) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = c(0.5, 8.5), clip = "off") +
    ggplot2::scale_y_continuous(
      breaks = prediction_key$y_id,
      labels = prediction_key$y_lab
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::scale_color_manual(
      values = c(
        r_choice = "#5C9ED6", r_out = "#5C9ED6",
        pe_choice = "#8B5FBF", pe_out_q = "#8B5FBF", pe_out_v = "#8B5FBF",
        cc_choice = "#D9534F", cc_out = "#D9534F",
        resid = "#414141"
      ),
      guide = "none"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#DCDCDC", linewidth = 0.47),
      axis.text.x = x_text,
      axis.text.y = y_text,
      axis.ticks.x = x_ticks,
      axis.ticks.y = y_ticks,
      axis.line.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

# Density/violin plot for recovered vector lengths.
plot_vector_recovery_density <- function(
  length_df,
  ref_df = vector_recovery_reference_lengths(),
  x_limits = c(0, 1),
  hide_axis_text = TRUE
) {
  prediction_levels <- c("resid", "cc_out", "pe_out_v", "pe_out_q", "r_out", "cc_choice", "pe_choice", "r_choice")
  prediction_labels <- c("Residual", "CC_out", "PE_out,V", "PE_out,Q", "R_out", "CC_choice", "PE_choice", "R_choice")
  condition_levels <- c("reward", "pe", "cc", "blend")
  condition_labels <- c("Reward", "PE", "CC", "Blend")

  length_plot_df <- length_df |>
    dplyr::mutate(
      prediction_raw = prediction,
      prediction = factor(prediction, levels = prediction_levels, labels = prediction_labels),
      condition = factor(condition, levels = condition_levels, labels = condition_labels)
    )

  ref_plot_df <- ref_df |>
    dplyr::mutate(
      prediction = factor(prediction, levels = prediction_levels, labels = prediction_labels),
      condition = factor(condition, levels = condition_levels, labels = condition_labels)
    )

  x_text <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_text()
  y_text <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_text()
  x_ticks <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_line(color = "black", linewidth = 0.47)
  y_ticks <- if (hide_axis_text) ggplot2::element_blank() else ggplot2::element_line(color = "black", linewidth = 0.47)

  ggplot2::ggplot(
    length_plot_df,
    ggplot2::aes(x = length, y = prediction)
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.47) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = prediction_raw),
      orientation = "y",
      trim = TRUE,
      scale = "width",
      adjust = 1,
      width = 0.82,
      alpha = 0.35,
      color = NA
    ) +
    ggplot2::geom_point(
      data = ref_plot_df,
      ggplot2::aes(x = target, y = prediction),
      inherit.aes = FALSE,
      shape = 124,
      size = 4.2,
      stroke = 1.3,
      color = "#22A652"
    ) +
    ggplot2::facet_wrap(~condition, ncol = 2) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::scale_fill_manual(
      values = c(
        r_choice = "#5C9ED6", r_out = "#5C9ED6",
        pe_choice = "#8B5FBF", pe_out_q = "#8B5FBF", pe_out_v = "#8B5FBF",
        cc_choice = "#D9534F", cc_out = "#D9534F",
        resid = "#414141"
      ),
      guide = "none"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#DCDCDC", linewidth = 0.47),
      axis.text.x = x_text,
      axis.text.y = y_text,
      axis.ticks.x = x_ticks,
      axis.ticks.y = y_ticks,
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.47),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

# Dispatcher: choose dot-cloud or density rendering.
# dot_* arguments apply only when plot_type = "dot".
plot_vector_recovery <- function(
  length_df,
  ref_df = vector_recovery_reference_lengths(),
  x_limits = c(0, 1),
  hide_axis_text = TRUE,
  plot_type = "dot",
  dot_alpha = 0.06,
  dot_size = 1.50,
  dot_jitter_height = 0.12
) {
  # Validate/normalize plot_type against the allowed options in c("dot", "density").
  # This keeps user choices ("dot" or "density") and errors on invalid values.
  plot_type <- match.arg(plot_type, c("dot", "density"))
  if (plot_type == "dot") {
    plot_vector_recovery_dot(
      length_df = length_df,
      ref_df = ref_df,
      x_limits = x_limits,
      hide_axis_text = hide_axis_text,
      dot_alpha = dot_alpha,
      dot_size = dot_size,
      dot_jitter_height = dot_jitter_height
    )
  } else {
    plot_vector_recovery_density(
      length_df = length_df,
      ref_df = ref_df,
      x_limits = x_limits,
      hide_axis_text = hide_axis_text
    )
  }
}

# Builds both panel and single-condition plots from a precomputed estimates data frame.
build_vector_recovery_plots <- function(
  length_df,
  ref_df = vector_recovery_reference_lengths(),
  x_limits = c(0, 1),
  hide_axis_text = TRUE,
  plot_type = "dot",
  dot_alpha = 0.06,
  dot_size = 1.50,
  dot_jitter_height = 0.12
) {
  # Validate/normalize plot_type against the allowed options in c("dot", "density").
  plot_type <- match.arg(plot_type, c("dot", "density"))
  panel_plot <- plot_vector_recovery(
    length_df = length_df,
    ref_df = ref_df,
    x_limits = x_limits,
    hide_axis_text = hide_axis_text,
    plot_type = plot_type,
    dot_alpha = dot_alpha,
    dot_size = dot_size,
    dot_jitter_height = dot_jitter_height
  )

  condition_plots <- lapply(
    split(length_df, length_df$condition),
    function(d) {
      plot_vector_recovery(
        length_df = d,
        ref_df = ref_df[ref_df$condition == unique(d$condition), , drop = FALSE],
        x_limits = x_limits,
        hide_axis_text = hide_axis_text,
        plot_type = plot_type,
        dot_alpha = dot_alpha,
        dot_size = dot_size,
        dot_jitter_height = dot_jitter_height
      ) + ggplot2::theme(legend.position = "none")
    }
  )

  list(panel_plot = panel_plot, condition_plots = condition_plots)
}

# Convenience wrapper: read an estimates CSV and return panel/condition plots.
build_vector_recovery_plots_from_file <- function(
  estimates_file,
  ref_df = vector_recovery_reference_lengths(),
  x_limits = c(0, 1),
  hide_axis_text = TRUE,
  plot_type = "dot",
  dot_alpha = 0.06,
  dot_size = 1.50,
  dot_jitter_height = 0.12
) {
  # Validate/normalize plot_type against the allowed options in c("dot", "density").
  plot_type <- match.arg(plot_type, c("dot", "density"))
  length_df <- utils::read.csv(estimates_file, stringsAsFactors = FALSE)
  build_vector_recovery_plots(
    length_df = length_df,
    ref_df = ref_df,
    x_limits = x_limits,
    hide_axis_text = hide_axis_text,
    plot_type = plot_type,
    dot_alpha = dot_alpha,
    dot_size = dot_size,
    dot_jitter_height = dot_jitter_height
  )
}

# Summarizes how often recovered lengths are exactly zero
# for predictions whose true target length is zero.
#
# length_df:
#   Long-format output with columns condition, prediction, and length.
#   This comes from run_s2_vector_recovery()$estimates or a saved estimates CSV.
# ref_df:
#   Reference target table from vector_recovery_reference_lengths().
# include_resid:
#   If FALSE, the residual row is excluded from zero-recovery summaries.
summarize_zero_recovery <- function(
  length_df,
  ref_df = vector_recovery_reference_lengths(),
  include_resid = FALSE
) {
  # Keep only prediction rows whose true target length is zero.
  zero_targets <- ref_df |>
    dplyr::filter(target == 0)

  # Optionally remove residual from evaluation.
  if (!include_resid) {
    zero_targets <- zero_targets |>
      dplyr::filter(prediction != "resid")
  }

  # Restrict recovered draws to zero-target prediction rows
  # and create an indicator for "recovered as exactly zero".
  check_df <- length_df |>
    dplyr::inner_join(
      zero_targets[, c("condition", "prediction"), drop = FALSE],
      by = c("condition", "prediction")
    ) |>
    dplyr::mutate(is_zero = length == 0)

  # Summary 1: percentage recovered as zero within each condition x prediction row.
  by_prediction <- check_df |>
    dplyr::group_by(condition, prediction) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_zero = sum(is_zero),
      pct_zero = 100 * mean(is_zero),
      .groups = "drop"
    )

  # Summary 2: percentage recovered as zero pooled across predictions within each condition.
  by_condition <- check_df |>
    dplyr::group_by(condition) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_zero = sum(is_zero),
      pct_zero = 100 * mean(is_zero),
      .groups = "drop"
    )

  # Summary 3: percentage recovered as zero pooled across all evaluated rows.
  overall <- check_df |>
    dplyr::summarise(
      n = dplyr::n(),
      n_zero = sum(is_zero),
      pct_zero = 100 * mean(is_zero)
    )

  list(
    by_prediction = by_prediction,
    by_condition = by_condition,
    overall = overall
  )
}

# Main wrapper for Study 2 vector-length recovery.
# n_sims is the number of simulations per generating condition.
# Total simulated datasets = 4 * n_sims (reward, pe, cc, blend).
# plot_type controls plotting style for returned figures: "dot" (default) or "density".
# dot_* arguments pass through to the dot-plot renderer when plot_type = "dot".
run_s2_vector_recovery <- function(
  trials,
  model_out_dir,
  n_sims = 250,
  seed = 19,
  plot_type = "dot",
  dot_alpha = 0.06,
  dot_size = 1.50,
  dot_jitter_height = 0.12
) {
  # Validate plot_type against the allowed options
  plot_type <- match.arg(plot_type, c("dot", "density"))
  set.seed(seed)

  # Ensure required columns exist before expensive looping begins.
  needed_cols <- c(
    "sub_index", "decrate_first", "overall_trial_nl", "trial_nl", "block_cent", "trial_nl_cent",
    "fA_ix", "fB_ix", "choice_numeric", "chosen_frac", "unchosen_frac",
    "out_a", "out_b", "chosen_out", "unchosen_out"
  )
  missing_cols <- setdiff(needed_cols, names(trials))
  if (length(missing_cols) > 0) {
    stop(paste0("Missing required trial columns: ", paste(missing_cols, collapse = ", ")))
  }

  fsml <- read_fsml("combined_shrink", model_out_dir = model_out_dir)
  sum_df <- fsml$sum

  # Residual SDs used when generating simulated ratings.
  d_resid <- get_median(sum_df, "d_resid")
  f_resid <- get_median(sum_df, "f_resid")

  # Build a subject-level scaffold from trial data.
  # decrate_first is retained because it determines block-effect construction.
  sub_df <- trials |>
    dplyr::distinct(sub_index, decrate_first) |>
    dplyr::arrange(sub_index)

  # Build scheduled (non-skip-dependent) probe numbers for simulation.
  # Decision ratings are scheduled in one block per subject (determined by decrate_first);
  # feedback ratings are scheduled in the other block.
  trials_sim <- trials |>
    dplyr::arrange(sub_index, overall_trial_nl) |>
    dplyr::mutate(
      dec_sched = (decrate_first == 1 & block == 1) | (decrate_first == 0 & block == 2),
      feed_sched = !dec_sched
    ) |>
    dplyr::group_by(sub_index) |>
    dplyr::mutate(
      dec_probe_number_sched = ifelse(dec_sched, cumsum(dec_sched), 0L),
      feed_probe_number_sched = ifelse(feed_sched, cumsum(feed_sched), 0L)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-dec_sched, -feed_sched)

  conditions <- c("reward", "pe", "cc", "blend")
  # Pre-allocate list containers for speed; one entry per simulated dataset.
  results <- vector("list", length = length(conditions) * n_sims)
  idx <- 1

  for (cond in conditions) {
    for (sim_i in seq_len(n_sims)) {
      # 1) Draw subject parameters and condition-specific true prediction lengths.
      subj_params <- s2_draw_subject_params(sum_df = sum_df, sub_df = sub_df)
      pred_len_draws <- draw_prediction_lengths(n_s = nrow(sub_df), condition = cond)
      rl_effect_draws <- prediction_lengths_to_rl_effects(pred_len_draws)

      # 2) Simulate ratings and fit recovery regressions.
      sim_df <- s2_simulate_affect_dataset(
        trials = trials_sim,
        subj_params = subj_params,
        rl_effects = rl_effect_draws,
        d_resid = d_resid,
        f_resid = f_resid
      )

      fit_list <- fit_recovery_models(sim_df)
      est_eff <- extract_effect_vector(fit_list)
      est_len <- prediction_lengths_from_effects(est_eff)

      # 3) Store long-format outputs for easy summarization/plotting.
      results[[idx]] <- data.frame(
        condition = cond,
        sim = sim_i,
        prediction = names(est_len),
        length = as.numeric(est_len),
        stringsAsFactors = FALSE
      )

      idx <- idx + 1
    }
  }

  # Stitch together outputs and construct plots.
  est_df <- dplyr::bind_rows(results)
  ref_df <- vector_recovery_reference_lengths()
  plot_bundle <- build_vector_recovery_plots(
    length_df = est_df,
    ref_df = ref_df,
    x_limits = c(0, 1),
    hide_axis_text = TRUE,
    plot_type = plot_type,
    dot_alpha = dot_alpha,
    dot_size = dot_size,
    dot_jitter_height = dot_jitter_height
  )
  zero_summary <- summarize_zero_recovery(
    length_df = est_df,
    ref_df = ref_df,
    include_resid = FALSE
  )

  list(
    estimates = est_df,
    zero_off_summary = zero_summary,
    panel_plot = plot_bundle$panel_plot,
    condition_plots = plot_bundle$condition_plots
  )
}
