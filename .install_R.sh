#!/usr/bin/env bash
while IFS=" " read -r package; 
do 
  RScript -e "if('"$package"' %in% rownames(installed.packages()) == FALSE) {install.packages('"$package"')}";
done < ".requirements_R.txt"