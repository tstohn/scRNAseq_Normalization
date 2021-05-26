#!/usr/bin/env bash
while IFS=" " read -r package; 
do	
  Rscript -e "if("'"$package"'" %in% rownames(installed.packages()) == FALSE) {install.packages(\"$package\", repo="'"http://cran.rstudio.com/"'")}";
done < ".requirements_R.txt"
