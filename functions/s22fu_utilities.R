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