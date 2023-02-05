# Dalhousie Fisheries Module


### To do before first class

1.  Install the program [[R](https://cran.r-project.org/)], there are options for Windows, Mac, and Linux users.  
2.  Install the program [R-studio](https://www.rstudio.com/products/rstudio/download/#download).  Again there are options for Windows, Mac, and Linux users.  
3.  Install the program [JAGS](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/), select version 4.3.  Again there are options for Windows, Mac, and Linux users.
4. *If using Windows* install the program [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) which is needed to install some R packages. This should not be necessary if using a Mac or Linux OS.

### Once you have these programs installed

1.  Open Rstudio  
    - RStudio is a Graphical User Interface (GUI) that makes using R much more intuitive.  
2.  In RStudio you will need to install a number of *packages*.  
    - *Packages* are additional tools that are developed by R users and they broaden the scope of what you can do within R.  
3.  The code needed to install these packages is given below
    - Simply copy this code and paste it into the *Console* window in R, (hit enter if it doesn't run immediately).
    - Lots of things you don't need to worry about will start flashing across your screen as the installation moves forward.
4.  If you already use R, still run the below code, it will skip over any packages you already have installed and simply install the packages you need for the course

```r
req.packages <- c("tidyverse","dplyr","sf","units","cowplot","knitr",'ggthemes',"marmap","RandomFields","RCurl",'readr','s2','plotly',"ggplot2","stars","tmaptools","maptools","rnaturalearth","rnaturalearthdata","raster","rgdal","RStoolbox","pals","ggnewscale","ggspatial",'devtools','raster','bookdown','gridExtra','scales','matlab','R2jags','zipfR','RandomFields')
# If you don't have the packages install them + give a heads up that you are
new.packages <- req.packages[!(req.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) 
{
  cat(paste0("Heads up, I have to install these packages for this to work:", new.packages ))
  #wanna.install <- readline(prompt = "If you want to install these package(s) enter 'y': ")
  #if(tolower(wanna.install) == 'y') 
  install.packages(new.packages,repos = "http://cran.us.r-project.org") #else { stop("You didn't want to install the packages so this script does not work.")}
}
# This is needed for our maps.
devtools::install_github('ropensci/rnaturalearthhires')

# You may/maynot need tinytex(), but once the above installations finish, I'm suggesting that you install it.
tinytex::install_tinytex()
# to uninstall TinyTeX, run tinytex::uninstall_tinytex() 

```
