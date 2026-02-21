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

library(ADNIMERGE2)

#View package documentation and vignettes

help(package = "ADNIMERGE2")

