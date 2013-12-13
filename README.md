###Overview
---------

This project contains all the scripts require to recreate the data and analyses present in the paper Isaac *et al* (2014). Range change simulations. *This Journal* xxx-xxx. In this paper we use simulated recording data to test methods of estimating range change under a range of recording scenarios.

###How do I run this myself?
---------------------------

We recommend you use RStudio since this has a good interface with GitHub. You also need Git, this allows RStudio to talk to GitHub. [Install git](http://git-scm.com/downloads), first, and oncew you have this installed [install RStudio](http://www.rstudio.com/ide/download/). If you don't have R you will need to [download and install](http://cran.r-project.org/) that too.

With R studio and Git installed, open RStudio and go to Project> Create project> Version Control> Git. Under 'Repository URL' paste the address for this page 'https://github.com/BiologicalRecordsCentre/RangeChangeSims'. 

If you have just installed git for the first time you might get an error from RStudio saying it can't find git. To fix this go to Tools> Options> Git/SVN> and browse for the git.exe (usually in a folder such as 'program files/GIT/bin/git.exe').

Once you have created this project run the script 'Run_full_analysis.r' to recreate the majority of results presented in the paper.

The only results that are not produced using this script are those for the occupancy models. However, we do provide a script that describes the model used ('Occupancy_model.r').

###How do I use these methods with my own data?
-----------------------------------------------

Excluding the occupancy models the various methods used can be found in the R package 'sparta'. A version of this package is included in this repository, however we recommend you install the most recent version of this package from [sparta's repository on GitHub](https://github.com/biologicalrecordscentre/sparta). There you will also be able to read tutorials and ask questions as well as report problems. 
