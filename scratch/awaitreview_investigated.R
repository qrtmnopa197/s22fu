qual <- read.csv("/Users/dp/Downloads/survey.csv")
prolif <- read.csv("/Users/dp/Downloads/prolific.csv")
prolif_nocodes <- prolif %>% filter(Status=="AWAITING REVIEW") %>% filter(Participant.id != "61006b0f1f3bfddb6b4c5d55")
qual_nocodes <- qual %>% filter(id %in% prolif_nocodes$Participant.id)
subs <- read.csv("/Users/dp/projects/s22_follow_up/analysis_data/sub_level_data_all_subs_2023-06-11_14_53_12.csv")
subs_nocodes <- subs %>% filter(id %in% prolif_nocodes$Participant.id)



#get ps and Ids
paths <- list.files("/Users/dp/projects/s22_follow_up/input_data",full.names = TRUE) %>% grep("/[0-9][^_/]*_(s22|arl|arlbv|s22fu)_task.*.csv",.,value = T)
ps_ids <- data.frame(ids=rep(NA,length(paths)),ps=rep(NA,length(paths)))
row <- 0

for(p in paths){
  df <- read_csv(p) #read in data
  if(nrow(df) > 0){
    #if this is a dataset
    row <- row + 1
    if(is.na(df$id[1])){
      ps_ids$id[row] <- NA
    } else {
      ps_ids$id[row] <- df$id[1]
    }
    
    if(is.na(df$participant[1])){
      ps_ids$ps[row] <- NA
    } else{
      ps_ids$ps[row] <- df$participant[1]
    }
  }
}

view(ps_ids %>% filter(id %in% prolif_nocodes$Participant.id))
write.csv(ps_ids,"/Users/dp/projects/s22_follow_up/analysis_data/ps_ids.csv")
