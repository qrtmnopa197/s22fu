#This script pushes all changes in the repository to the master branch, not pulling first on the assumption that no one else is making changes to the repository.
#To work on server's that are not DP's Mac, change the path in the first line.

cd /Users/dp/projects/s22_follow_up/code
git add .
git commit -m "$1"
git push origin main

