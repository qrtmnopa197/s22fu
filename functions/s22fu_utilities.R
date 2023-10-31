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
    feed_rate_0 = feed_rate_0
  )
  return(data)
}

#adds a column for previous valence rating to the trial-level dataset, returning that dataset
add_prevrate <- function(df,rat_col_name){
  df$prev_rate <- 0
  #for each trial in the block past the first...
  for(t in 2:nrow(df)){
    #if there was a rating on the previous trial...
    if(!is.na(df[t-1,rat_col_name])){
      df$prev_rate[t] <- df[t-1,rat_col_name] #that's the prev_rate for this trial...
    } else{
      #otherwise, the prev_rate is the previous rating from the previous trial (the most recent rating)
      df$prev_rate[t] <- df$prev_rate[t-1]
    } 
  }
  return(df)
}