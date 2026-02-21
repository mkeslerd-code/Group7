# ---- Git quick push script ----

system("git add .")
system('git commit -m "Add project files"')
system("git push origin main")


# To use ADNI MERGE Library

#Install Hmisc package: The ADNIMERGE package depends on the Hmisc package, which can be installed from within R by running
install.packages(c("assertr", "bookdown", "Hmisc", "tidyverse"))

#Replace "/path/to/ADNIMERGE_0.0.1.tar.gz" with the actual location and filename of the file you downloaded. ( I have uploaded the file into our Git Repository to save us steps)
install.packages("/Users/makayla/Captsone project/Group7/ADNIMERGE2.tar.gz", repos = NULL, type = "source")

#Load the package

library(ADNIMERGE)

#View package documentation and vignettes

help(package = "ADNIMERGE")

#Load main data frame - This function loads the consolidated ADNI data set into an analysis-ready R data frame.
adnimerge_data <- adnimerge()

#View documentation for a specific data table (e.g., adas)

?adas
