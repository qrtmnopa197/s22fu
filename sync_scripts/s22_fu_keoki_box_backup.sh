#Backs up the entire spring_2022_study folder to keoki, and then Box.

set -e #exit the script if there are errors; this guarantees that both syncs were successful if the script runs through

rclone copy /Users/dp/projects/s22_follow_up Duke_box:/PROJECT\ 3373\:\ Projects/s22_follow_up -v #copies to Box