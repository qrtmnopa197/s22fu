#This is the master script for initial data analysis steps
#It identifies data to use, wrangles it into a usable form, creates additional variables based on this data, and creates plots of certain variables for quality-checking.

##SET MANUALLY
path_to_project_directory <- "~/projects/s22_follow_up/"
path_to_s22 <- "~/projects/spring_2022_study/"
path_to_arlbv <- "~/projects/ARL_BV/"
ids_to_exclude <- c(1:39,522,142) #Ps whose data you don't want to analyze even if it looks good. Exclude 86041 - the test code - and 1 - obviously a pilot ID code. 
##############

#clear out results from old analyses
system("mv /Users/dp/projects/s22_follow_up/analysis_data/*.csv /Users/dp/projects/s22_follow_up/analysis_data/old_analysis_data") #move any CSV files in this folder to old_analysis_data folder
system("rm -rf /Users/dp/projects/s22_follow_up/analysis_data/qc_plots/trial_level_vars/*") #clear the folder containing trial-level plots
system("rm -f /Users/dp/projects/s22_follow_up/analysis_data/qc_plots/sub_level_variables.pdf") #remove old subject level variable plot

source(paste0(path_to_s22,"code/functions/s22_utilities.R")) #get s22 functions
source(paste0(path_to_project_directory,"code/functions/s22fu_wrangle_psychopy_data.R")) 
source(paste0(path_to_project_directory,"code/functions/s22fu_utilities.R")) 
source(paste0(path_to_arlbv,"code/functions/arlbv_utilities.R")) #get ARL-BV functions

csvs_to_analyze <- get_csvs_to_analyze(path_to_project_directory,ids_to_exclude) #get the usable CSVs to analyze
all_data <- lapply(csvs_to_analyze, s22fu_wrangle_psychopy_data) #reformats each CSV, turning it into a long dataset usable for analysis, and adds all variables of interest that can be created from the raw data alone.
                                                               #Returns a list - one element for each subject - where each element is itself a list containing dfs with the trial-level data and subject-level data
trials_list <- lapply(all_data, function(l) l[[1]]) #get a list of the trial-level dfs only
trial_level_data <- do.call(rbind,trials_list) #stack trial dfs into one big df
sub_list <- lapply(all_data, function(l) l[[2]]) #get a list of the subject-level dfs only
sub_level_data <- do.call(rbind,sub_list) #stack them into one big data frame
#write both to CSVs
date_time <- Sys.time() %>% chartr(" ","_",.) %>% chartr(":","_",.) #grab the date and time, reformatting ':' and '-' to  '_' so you can label the files with it
write.csv(trial_level_data,paste0(path_to_project_directory,"analysis_data/trial_level_data_all_subs_",date_time,".csv"),row.names = FALSE)
#write subject-level data
write.csv(sub_level_data,paste0(path_to_project_directory,"analysis_data/sub_level_data_all_subs_",date_time,".csv"),row.names = FALSE)
#HERE
#Create plot grids of select variables for quality checking. These are saved to the qc_plots folder
tl_hists <- c("decrate_rt","feedrate_rt", "choice_rt","dec_rate","feed_rate") #trial level variables to plot
plot_trial_level_vars(trial_level_data,tl_hists,path_to_project_directory) #create and save plot grids
sl_hists <- c("consecutive_late_choices","instruct_keypress","consecutive_autoprocess",
              "choice_pt_completed","decrate_pt_completed","feedrate_pt_completed",
              "answers_correct","answers_incorrect","earnings_before_ec",
              "feed_check_passed","feedrate_check_passed","decrate_check_passed","att_checks_passed",
              "percent_left","percent_left_b1","percent_left_b2","late_percent","trials_completed","decrate_skipped_percent","feedrate_skipped_percent","valrate_skipped_percent",
              "decrate_sd","feedrate_sd","valence_sd","total_experiment_time") #subject-level variables to plot
plot_sub_level_vars(sub_level_data,sl_hists,path_to_project_directory) #create and save plot grids




