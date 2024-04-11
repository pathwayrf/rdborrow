# Change this to your local directory containing the myproj1 folder
mydir="/home/zhoux104/R/xinercode/"
#mydir="/Users/xiner/Library/CloudStorage/Dropbox/rdborrow_pre/"




# Function to check, install, and load packages
check_and_install_load_packages <- function(package_list) {
  # Check installed packages
  installed_packages <- rownames(installed.packages())
  
  # Identify missing packages
  missing_packages <- package_list[!(package_list %in% installed_packages)]
  
  # Install missing packages
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  } else {
    message("All packages are already installed.")
  }
  
  # Load the packages
  for (pkg in package_list) {
    library(pkg, character.only = TRUE)
  }
}

# List of packages to check, install, and load
packages_to_install_load <- c('devtools',"ggplot2","ggpubr", "dplyr", "tidyr","nnet","tableone","kableExtra","boot","batchtools","utils",
                              "reshape2","ggsci","knitr","xtable","latex2exp","patchwork","ggridges","ggh4x","CVXR","mmrm","cobalt","WeightIt",
                              "parallel","mice","rbmi")



# Call the function with the list of packages
check_and_install_load_packages(packages_to_install_load)


### RBMI

# if using RefBasedMI
#devtools::install_github("UCL/RefBasedMI")
#library(RefBasedMI)
 
# if using rbmi and need to update to latest version
# remove previous install
# library_path <- find.package("rbmi")
# remove.packages("rbmi", lib = dirname(library_path))
# # install the new version
# install.packages(paste(mydir,"myproj2/rbmi_1.2.3.tar.gz",sep=''), repos = NULL, type="source",
#                  lib=.libPaths()[2])
# library(rbmi)
 

