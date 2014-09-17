###Overview
---------

This project contains all the scripts required to recreate the data and analyses present in the paper Isaac *et al* (2014). Range change simulations. *This Journal* xxx-xxx. In this paper we use simulated recording data to test methods of estimating range change under a range of recording scenarios.

###How do I run this myself?
---------------------------

Full details can be found in the vignette, which you can [download via this link](https://github.com/BiologicalRecordsCentre/RangeChangeSims/blob/master/vignette.pdf?raw=true).

We recommend you use RStudio since this has a good interface with GitHub. You also need Git, this allows RStudio to talk to GitHub. [Install git](http://git-scm.com/downloads), first, and once you have this installed [install RStudio](http://www.rstudio.com/ide/download/). If you don't have R you will need to [download and install](http://cran.r-project.org/) that too. Finally, to run occupancy models we need a programme called JAGS installed. This can be downloaded from [sourceforge](http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/)

With R studio and Git installed, open RStudio and go to File> New project> Version Control> Git. Under 'Repository URL' paste the address for this page 'https://github.com/BiologicalRecordsCentre/RangeChangeSims'. 

If you have just installed git for the first time you might get an error from RStudio saying it can't find git. To fix this go to Tools> Options> Git/SVN> and browse for the git.exe (usually in a folder such as 'program files/GIT/bin/git.exe').

Once you have created this project run the script 'Run_full_analysis.r' to recreate the majority of results presented in the paper.


###How do I use these methods with my own data?
-----------------------------------------------

The various modelling methods used can be found in the R package 'sparta'. A version of this package is included in this repository, however we recommend you install the most recent version of this package from [sparta's repository on GitHub](https://github.com/biologicalrecordscentre/sparta). There you will also be able to read tutorials and ask questions as well as report problems. 
