#This function takes a single subject's raw psychopy data and reformats it into a long dataset usable for analysis. 
#It creates all variables of interest that can be created from the raw data alone.
s22fu_wrangle_psychopy_data <- function(csv_path){
  library(ddpcr)
  library(tidyverse)
  library(sigmoid)
  if(!is.function(get_csvs_to_analyze)){
    stop("get_csvs_to_analyze() is not loaded into the environment of arlbv_wrangle_psychopy_data(), suggesting that the functions in s22_utilities.R are not present in the function's environment.")
  }
  print(csv_path) #print the csv path for debugging purposes - so you know what subject you're on

  df_full <- read_csv(csv_path) #read in subject's data
  #if the subject never pressed the right arrow key on affect probes, the RT columns won't appear. Thus, you need to create fully empty RT columns for them.

  
  #make one df with all the subject-level info, and a second with all the trial-level info
  sub_info <- select(df_full, any_of(c(id = "id", "age_gender_formresponse", "date", "earnings_before_ec", "earnings",total_experiment_time = "total_experiment_time", 
                                       comp_qs = "mouse_3.clicked_name",instruct_keypress = "instruct_key.keys",
                                       feed_check = "pressq_feed.keys", decrate_check = "pressq_decrate.keys", feedrate_check = "pressq_feedrate.keys"))) 
  trials <- select(df_full, any_of(c("fA_img","fB_img", "fA_win_prob","fB_win_prob",choice_rt = "mouse.rt", too_slow = "too_slow",choice = "click",
                                     "out_a","out_b",block_decrate = "dec_rate",block_feedrate = "feed_rate", decrate_raw = "val_slider_7.response", feedrate_raw = "val_slider_8.response", 
                                     decrate_rt = "dec_rate.rt", feedrate_rt = "feed_rate.rt", 
                                     trial_raw = "trials.thisN", block_raw = "blocks.thisN", makeup_repetition = "trials_and_makeups.thisN",
                                     pt_choice_n = "pt_choice_outer.thisN",pt_fdaff_n = "pt_fdaff_outer.thisN",pt_dcaff_n = "pt_dcaff_outer.thisN",
                                     choice_instruct_n = "choice_instruct.thisN", affect_instruct_n = "affect_instruct.thisN",instruct_keys = "instruct_key.keys")),
                   ends_with(".ran"))
  
  #if they never clicked the continue button on a valence rating, fill in that column with NAs
  if(!("feedrate_rt" %in% names(trials))){
    trials$feedrate_rt <- NA
  }
  if(!("decrate_rt" %in% names(trials))){
    trials$decrate_rt <- NA
  }
  
  #Count number of practice trials completed for the sub-level df. This must be done before the trials df is reformatted
  #get the last row on which each set of practice trials was comleted
  pt_choice_last <- which(trials$pt_choice_n == 2)
  pt_fdaff_last <- which(trials$pt_fdaff_n == 2)
  pt_dcaff_last <- which(trials$pt_dcaff_n == 2)
  
  trials_ptchoice <- trials[1:pt_choice_last,] #get all rows before the end of the choice PTs
  choice_pt_completed <- length(which(!is.na(trials_ptchoice$choice)))#count how many many choices were made in these PTs
  #get all rows between the end of the choice PTs and the fdaff PTs. Repeat the process...
  trials_ptfdaff <- trials[(pt_choice_last+1):pt_fdaff_last,]
  feedrate_pt_completed <- length(which(!is.na(trials_ptfdaff$choice)))
  #repeat the process for dcaff PTs...
  trials_ptdcaff <- trials[(pt_fdaff_last+1):pt_dcaff_last,]
  decrate_pt_completed <- length(which(!is.na(trials_ptdcaff$choice)))
  
  #Similarly, calculate the longest string of instructions trial on which the subject did not press the right arrow key.
  trials_choice_inst <- trials[which(trials$choice_instruct_n==0):which(trials$choice_instruct_n == max(trials$choice_instruct_n,na.rm=TRUE)),] #grab the rows corresponding to choice instruction trials
  inst_key_runs <- rle(is.na(trials_choice_inst$instruct_keys)) #get the runs of NAs and non-NAs (NA runs are TRUE)
  #if some key presses were NA
  if(any(inst_key_runs$values)){
    na_indices <- which(inst_key_runs$values) #get the indices that correspond to runs of NA
    consecutive_autoprocess_choice <- max(inst_key_runs$lengths[na_indices]) #get the longest NA run
  } else{
    consecutive_autoprocess_choice <- 0 #if none of the runs are NA, then the subject had no late choices
  }
  #do the same thing for the affect instructions
  trials_affect_inst <- trials[which(trials$affect_instruct_n==0):which(trials$affect_instruct_n == max(trials$affect_instruct_n,na.rm=TRUE)),] #grab the rows corresponding to choice instruction trials
  inst_key_runs <- rle(is.na(trials_affect_inst$instruct_keys)) #get the runs of NAs and non-NAs (NA runs are TRUE)
  #if some key presses were NA
  if(any(inst_key_runs$values)){
    na_indices <- which(inst_key_runs$values) #get the indices that correspond to runs of NA
    consecutive_autoprocess_affect <- max(inst_key_runs$lengths[na_indices]) #get the longest NA run
  } else{
    consecutive_autoprocess_affect <- 0 #if none of the runs are NA, then the subject had no late choices
  }
  consecutive_autoprocess <- max(c(consecutive_autoprocess_affect,consecutive_autoprocess_choice))
  
  trials <- trials %>% select(-c(pt_choice_n,pt_fdaff_n,pt_dcaff_n,choice_instruct_n,affect_instruct_n,instruct_keys)) #don't need these anymore
  
  #REFORMAT THE TRIALS DF...
  #Find the rows containing the first and last trials, getting rid of everything above and below those...
  aff_qs_run <- which(trials$aff_qs2.ran == 1) #identify the rows on which aff_qs2 - the last loop before the main task - ran
  main_task_start <- max(aff_qs_run) + 1 #the row after the last run of this loop should be the first trial
  block_runs <- which(trials$blocks.ran == 1) #identify all rows on which the blocks loop ran. 
  main_task_end <- max(block_runs) #the last run of this loop signifies the end of the main task
  trials <- trials[main_task_start:main_task_end,] #keep only rows within the main task
  
  #Create columns for block number and makeup/repetition a trial represents, as well as what the fA_img and fB_img were
  #These things are all are output at the end of each block/repetition loop,
  #so you need to fill out all the rows above each loop end
  #with that value. You can use the custom "fill_vec" function to do this (see s22_functions.R in this folder)
  trials$fA_img <- fill_vec(trials$fA_img, bottom_up = TRUE)
  trials$fB_img <- fill_vec(trials$fB_img, bottom_up = TRUE)
  trials$block_raw <- fill_vec(trials$block_raw, bottom_up = TRUE)
  trials$makeup_repetition <- fill_vec(trials$makeup_repetition, bottom_up = TRUE)
  trials$block_decrate <- fill_vec(trials$block_decrate, bottom_up = TRUE)
  trials$block_feedrate <- fill_vec(trials$block_feedrate, bottom_up = TRUE)
  
  # You want one row per trial, meaning you want to ratchet down one row at the end of every trial. 
  # Fortunately, the end of each trial is marked by a loop end, so psychopy does indeed ratchet down one row
  # at the end of every trial. However, there are a few loops that sometimes end in between trials - 
  #  i.e., before the trial starts, with the first trial always starting on the right row - 
  # leading to an unnecessary ratchet-down before the start of the next trial, before any trial data has been collected.
  # To resolve this issue, you can simply delete all rows on which these loops have run/the pre-trial ratchet-down 
  # has occurred (since no trial data is collected before the ratchet-down).
  # A second issue is that there’s a loop which sometimes ends in the middle of a trial (remaining_trial), 
  # which results in a single trial’s data being spread across two lines. 
  # To address this, you should identify trials in which there was a ratchet-down mid-trial. 
  # In these cases, you know that the current row and the row below it represent a single trial’s data, 
  # and that you need to combine them into one row. The simplest way to do this is to copy the trial data
  # on the second row to the first row, thus ensuring that the first row contains the full trial’s data. 
  # Then, delete the second row, which contains only redundant information.
  trials$delete_row <- 0
  for(row in 1:nrow(trials)){
    #On rows where a pre-trial ratchet-down has occurred...
    if (!is.na(trials$trials_and_makeups.ran[row]) | !(is.na(trials$blocks.ran[row]))){
      trials$delete_row[row] <- 1 #mark row for deletion
    } else if (!(is.na(trials$remaining_trial.ran[row]))){
      #Mid-trial ratchet-down. First, copy all the valuable data from the row below to the row on which the trial started.
      row_below <- row + 1
      trials$trial_raw[row] <- trials$trial_raw[row_below]
      trials$fB_img[row] <- trials$fB_img[row_below]
      trials$fA_img[row] <- trials$fA_img[row_below]
      if(any(names(trials) == "too_slow")){
        trials$too_slow[row] <- trials$too_slow[row_below]
      }
      #Now, mark the row below (containing only redundant data) for deletion
      trials$delete_row[row_below] <- 1 
    }
  }
  trials <- filter(trials, delete_row == 0) #actually delete the marked rows
  trials <- select(trials, -delete_row,-ends_with(".ran")) #get rid of .ran and delete_row columns
  
  
  #Do some variable recoding...
  
  #if the data has a "too slow" column in it...
  if(any(names(trials) == "too_slow")){
    #if the subject was "too slow" on a trial, mark the trial as late
    for(row in 1:nrow(trials)){
      if(trials$too_slow[row] == "TRUE"){
        trials$choice[row] <- "late"
      }
    }
    trials <- select(trials,-too_slow) #don't need this anymore
  }
  
  #python indexing starts at 0, so this rectifies that
  trials$trial <- trials$trial_raw + 1
  trials$block <- trials$block_raw + 1
  trials <- select(trials,-c(trial_raw,block_raw)) #don't need raw rows anymore
  
  #Trials actually give you the trial within the trials_and_makeups loop, not the trial within the block.
  #You can rectify this by going through the df top to bottom, making each repetition trial equal to the preceding trial plus one
  #NB:This only works assuming the trials are in ascending order and there are no gaps between trials. It's very difficult to do this
  #in a more principled way, so you'll just have to go with this for now.
  for(row in 1:nrow(trials)){
    #if you're on a repetition
    if(trials$makeup_repetition[row] != 0){
      trials$trial[row] <- trials$trial[row-1] + 1 #set the trial value to one more than the previous trial
    }
  }
  
  #Reduce the fA/B_img columns to just the fractal number - not the whole path to the JPEG
  trials <- mutate(trials, fA_img = str_extract(fA_img,"\\d+"), fB_img = str_extract(fB_img,"\\d+"))
  
  trials$row_index <- c(1:nrow(trials)) #get the row numbers for this subject, which will be useful for a number of things to follow
  
  #create choice numeric column, assigning a value of 1 for choosing fractal A and 2 for choosing fractal B
  trials$choice_numeric <- NA
  for(row in 1:nrow(trials)){
    if(trials$choice[row] == "fA"){
      trials$choice_numeric[row] <- 1
    } else if(trials$choice[row] == "fB"){
      trials$choice_numeric[row] <- 2
    }
  }
  
  #create a "stay" column, indicating whether the participant made the same choice the next time the pair was presented (if so, 1) or not (if so, 0)
  stay_data <- filter(trials,choice != "late") #just look at non-late choices (coding this is tricky if you look at everything)
  stay_data$stay <- NA #initialize stay column
  #go trial-by-trial...
  for(row in 1:nrow(stay_data)){
    pair_rows <- which(stay_data$fA_img == stay_data$fA_img[row]) #identify the rows on which the participant played with the same pair
    #if there are future rows on which they played with the same pair...
    if(any(pair_rows > row)){
      future_pres <- pair_rows[pair_rows > row] #identify all future trials/rows in which this pair was presented
      next_pres <- min(future_pres) #identify the next choice with this pair by taking the min of this vector
      if(stay_data$choice[next_pres] == stay_data$choice[row]){
        stay_data$stay[row] <- 1 #if the choice at the next presentation is the same as the current choice, the participant "stayed"
      } else if(stay_data$choice[next_pres] != stay_data$choice[row]){
        stay_data$stay[row] <- 0 #otherwise they switched
      }
    } 
  }
  #join this stay column back to the main trials df
  stay_join <- select(stay_data,row_index,stay) 
  trials <- left_join(trials,stay_join,by="row_index")
  
  trials <- trials %>% mutate(dec_rate = decrate_raw*-1 + 2,feed_rate = feedrate_raw*-1 + 2,) %>% select(-feedrate_raw,-decrate_raw) #recode valence to be on a 0-1 scale
  trials <- trials %>% mutate(dec_rate_z = dec_rate,feed_rate_z = feed_rate) #because col_zscore replaces the original columns, but I want to keep them, I'm duplicating them
  trials <- trials %>% col_zscore(c("dec_rate_z","feed_rate_z"))
  
  #create columns for chosen and unchosen outcomes
  trials <- trials %>% mutate(chosen_out = ifelse(choice == "fA",out_a,out_b),
                              unchosen_out = ifelse(choice == "fA",out_b,out_a))
  
  #get pair presentation numbers (for trials on which a choice was made)
  trials_nl <- trials %>% filter(choice != "late") #filter out late trials
  pair_pres_list <- by(trials_nl,trials_nl$fA_img,assign_pp_nums) #for each pair (identified by fA), assign a sequential vector of numbers to pair_pres, 
                                                            #returning the pair pres numbers and the row index
  pair_pres_df <- do.call(rbind,pair_pres_list)  #combine these into one big df
  trials <- left_join(trials,pair_pres_df,by="row_index") #join to the trials df
  
  
  trials <- trials %>% mutate(all_rate_z = ifelse(is.na(dec_rate_z),feed_rate_z,dec_rate_z)) #create single column of affect ratings
  block_list <- by(trials,trials$block,add_prevrate,rat_col_name="all_rate_z") #add prev_rate column to the df for each block
  trials <- do.call(rbind,block_list)
  trials$prev_rate <- as.numeric(trials$prev_rate)
    
  #NOW, THE SUBJECT INFO DF...
  age <- sub_info$age_gender_formresponse[grepl("\\d+",sub_info$age_gender_formresponse)] #identify row with the age by searching for a digit value
  gender <- sub_info$age_gender_formresponse[grepl("[a-zA-Z]+",sub_info$age_gender_formresponse)] #identify row with gender by searching for letter value
  if(length(age)==0){
    age <- 999
  }
  if(length(gender)==0){
    gender <- 999
  }
  earnings_before_ec <- sub_info$earnings_before_ec[grepl("\\d+",sub_info$earnings_before_ec)] #identify row with earnings before ec by searching for a digit value
  earnings <- sub_info$earnings[grepl("\\d+",sub_info$earnings)] #ditto for total earnings
  total_experiment_time <- sub_info$total_experiment_time[grepl("\\d+",sub_info$total_experiment_time)] #same strategy for total experiment time
  
  #count the number of each type of attention check passed
  feed_check_passed <- length(which(sub_info$feed_check == "q"))
  decrate_check_passed <- length(which(sub_info$decrate_check == "q"))
  feedrate_check_passed <- length(which(sub_info$feedrate_check == "q"))
  
  att_checks_passed <- feed_check_passed + decrate_check_passed + feedrate_check_passed #get the total
  
  id <- sub_info$id[1] #get id from one of the rows
  date <- sub_info$date[1] #get date from one of the rows
  
  #check the longest string of late choices the subject had. Long strings suggest that the participant left the task
  value_runs <- rle(trials$choice) #get the runs of each value
  #if late is one of the choice values...
  if(any(value_runs$values == "late")){
    late_indices <- which(value_runs$values == "late") #get the indices that correspond to runs of "late"
    consecutive_late_choices <- max(value_runs$lengths[late_indices]) #get the longest "late" run
  } else{
    consecutive_late_choices <- 0 #if none of the runs are "late"s, then the subject had no late choices
  }
  #get the standard deviations of valence ratings
  decrate_sd <- sd(trials$dec_rate,na.rm=TRUE) 
  feedrate_sd <- sd(trials$feed_rate,na.rm=TRUE) 
  valence_sd <- sd(c(trials$dec_rate,trials$feed_rate),na.rm=TRUE)
  
  instruct_keypress <- length(which(sub_info$instruct_keypress == "right")) #number of manual processions on instructions slides; 
                                                                            #a small number is a red flag they weren't paying attention
  
  answers_correct <- length(which(sub_info$comp_qs=="poly_true"))
  answers_incorrect <- length(which(sub_info$comp_qs=="poly_false"))
  
  trials_completed <- sum(trials$choice != "late") #get the total number of trials completed
  percent_left <- length(which(trials$choice == "fA"))/trials_completed  #times fA was chosen divided by total choices

  b1 <- filter(trials,choice != "late" & block == 1)
  percent_left_b1 <- length(which(b1$choice == "fA"))/nrow(b1)
  
  b2 <- filter(trials,choice != "late" & block == 2)
  percent_left_b2 <- length(which(b2$choice == "fA"))/nrow(b2)
    
  
  #get the percentage of trials on which the subject didn't make a choice in time
  late_numeric <- ifelse(trials$choice == "late",1,0)
  late_percent <- mean(late_numeric)
  
  #get the percentage of trials on which subjects made an affect rating
  trials_nl <- trials %>% filter(choice != "late") #provided affect rating on all non-late trials
  trials_decrate <- trials_nl %>% filter(block_decrate == 1)
  decrate_skipped_percent <- length(which(is.na(trials_decrate$dec_rate)))/nrow(trials_decrate)
  trials_feedrate <- trials_nl %>% filter(block_feedrate == 1)
  feedrate_skipped_percent <- length(which(is.na(trials_feedrate$feed_rate)))/nrow(trials_feedrate)
  valrate_skipped_percent <- (decrate_skipped_percent + feedrate_skipped_percent)/2
  
  ### Calculate some of the above variables for each block separately
  
  #late percent
  block1 <- filter(trials,block == 1)
  late_numeric_b1 <- ifelse(block1$choice == "late",1,0)
  late_percent_b1 <- mean(late_numeric_b1)
  
  block2 <- filter(trials,block == 2)
  late_numeric_b2 <- ifelse(block2$choice == "late",1,0)
  late_percent_b2 <- mean(late_numeric_b2)
  
  
  if(trials$block_decrate[1] == 1){
    decrate_first <- 1
  }else if(trials$block_feedrate[1] == 1){
    decrate_first <- 0
  }

  
  #CHECK YOU DIDN'T MISS ANYTHING
  sub_info_final <- data.frame(id,age,gender,date,earnings_before_ec,earnings,total_experiment_time,answers_correct,answers_incorrect,
                               feed_check_passed,feedrate_check_passed,decrate_check_passed,att_checks_passed,consecutive_late_choices,decrate_sd,feedrate_sd,
                               valence_sd,instruct_keypress,consecutive_autoprocess,percent_left,trials_completed,percent_left_b1,percent_left_b2,late_percent,decrate_skipped_percent,
                               feedrate_skipped_percent,valrate_skipped_percent,late_percent_b1,late_percent_b2,row.names=NULL,
                               choice_pt_completed,decrate_pt_completed,feedrate_pt_completed,decrate_first)
  
  base_data <- cbind(sub_info_final,trials) #merge the two dfs
  
  return(list(base_data,sub_info_final)) #return the trial-level data and subject-level-data
}
