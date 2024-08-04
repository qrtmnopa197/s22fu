library(dplyr)
df_full <- read.csv("/Users/dp/projects/s22_follow_up/analysis_data/trial_level_data_all_subs_2023-05-30_20_32_41.csv",na.strings = "")
df_full$fA_win_prob <- as.numeric(df_full$fA_win_prob)
df_full$fB_win_prob <- as.numeric(df_full$fB_win_prob)
11
12
13
14
17
21
3
4
5
7
9

df <- filter(df_full,id == 9)

get_trans <- function(df){
  vec <- c()
  for(row in 2:nrow(df)){
    delta <- df$fA_win_prob[row] - df$fA_win_prob[row-1]
    vec <- c(vec,delta)
  }
  return(vec)
}

fA_df <- select(df,fA_win_prob,fA_img)
fa_list <- by(fA_df,fA_df$fA_img,get_trans) 
fa_vec <- unlist(fa_list)

hist(fa_vec)
sd(fa_vec,na.rm=TRUE)


get_trans_b <- function(df){
  vec <- c()
  for(row in 2:nrow(df)){
    delta <- df$fB_win_prob[row] - df$fB_win_prob[row-1]
    vec <- c(vec,delta)
  }
  return(vec)
}

fB_df <- select(df,fB_win_prob,fB_img)
fb_list <- by(fB_df,fB_df$fB_img,get_trans_b) 
fb_vec <- unlist(fb_list)

hist(fb_vec)
sd(fb_vec,na.rm=TRUE)

