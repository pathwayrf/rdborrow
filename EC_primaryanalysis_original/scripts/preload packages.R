# Change this to your local directory containing the myproj1 folder
# mydir="/home/zhoux104/R/xinercode/"





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
                              "reshape2","ggsci","knitr","xtable","latex2exp","patchwork","ggridges","ggh4x")



# Call the function with the list of packages
check_and_install_load_packages(packages_to_install_load)



# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
# check_cmdstan_toolchain()
# install_cmdstan(release_url="https://github.com/stan-dev/cmdstan/releases/download/v2.32.1/cmdstan-2.32.1.tar.gz",
#                 dir ="/home/zhoux104/R/xinercode/myproj1")
# Set path to your CmdStan installation
# cmdstanr::set_cmdstan_path("/home/zhoux104/R/xinercode/myproj1/cmdstan-2.32.1")
#set_cmdstan_path(path="/home/zhoux104/.cmdstan")
# cmdstan_path()
    
