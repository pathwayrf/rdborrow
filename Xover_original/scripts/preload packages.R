
# Change this to your local directory containing the myproj1 folder
mydir="/home/zhoux104/R/xinercode/"





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
packages_to_install_load <- c("ggplot2", "dplyr", "tidyr","nnet","tableone","kableExtra","mmrm","boot","batchtools","utils",
                              "reshape2","ggsci","knitr","xtable","latex2exp","patchwork","ggridges","ggh4x","CVXR","rbmi","cobalt","WeightIt")



# Call the function with the list of packages
check_and_install_load_packages(packages_to_install_load)

