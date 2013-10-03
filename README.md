###Overview
---------

This project contains all the scripts require to recreate the data and analyses present in the paper Isaac *et al* (2014). Range change simulations. *This Journal* xxx-xxx. In this paper we use simulated recording data to test methods of estimating range change under a range of recording scenarios.

###How do I run this myself?
---------------------------

We recommend you use RStudio since this has a good interface with GitHub. Once you have [installed RStudio](http://www.rstudio.com/ide/download/) you will also need to [install git](http://git-scm.com/downloads). Git allows RStudio to talk to GitHub.

With R studio and Git installed, open RStudio and go to Project> Create project> Version Control> Git. Under 'Repository URL' paste the address for this page 'https://github.com/BiologicalRecordsCentre/RangeChangeSims'.

Once you have created this project run the script 'Run_full_analysis.r' to recreate the majority of results presented in the paper.

The only results that are not produced using this script are those for the occupancy models. However we do provide a script that describes the model used ('Occupancy_model.r').

###How do I use your methods with my own data?
----------------------------------------------

Excluding the occupancy models the various methods used can be found in the R package 'sparta'. This package is included in this repository. In addition sparta has its own repository where you can download the development version, red tutorials and ask questions as well as report problems. 
