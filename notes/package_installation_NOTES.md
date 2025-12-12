# R/package notes

## To do sometime:

R packages to update on the cluster (not yet, as of Dec 2025. wait for R 4.5 to be available via Rstudio)
- first: S4Vectors  0.47.6 
Then:
- IRanges 2.43.8 
- XVector 0.49.3

Updated both laptops after the major BioC release on 10/30, and now I have those new versions. However on the cluster, as of 12/11/2025, the most recent version of Rstudio is 4.4.0. They need to update the OS of the G-class nodes before they can let us have R 4.5 with R studio. That means that I'm restricted to an older Bioconductor release (<3.22), and that means I don't have access to the newer versions of the above packages.

## General notes on package installation

For some projects, I'll use `renv` to manage package versions. See [`renv` notes](https://github.com/jayoung/Rtest_and_Rnotes/blob/main/notes/renv_notes.md).

To install a CRAN package:
```
install.packages("ape", dependencies = TRUE)
```

To install a Bioconductor package:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

# or perhaps:
BiocManager::install("Biostrings", lib="/home/jayoung/malik_lab_shared/linux_gizmo/R_packages")
```

To install a package from a tar.gz file:
```
install.packages("GEOquery_2.1.8.tar.gz")
# or from outside R:
R CMD INSTALL GEOquery_2.1.8.tar.gz
```


In rare cases, a package is hosted elsewhere, not CRAN or Bioconductor:
```
install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
## needed to do this on Mac
install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
```

Packages can be installed from github too, using `devtools::install_github()`. Examples:

```
library(devtools)
install_github("easyGgplot2", "kassambara")
install_github("js229/Vennerable")
```

It used to be useful to run `chooseCRANmirror()` but I don't bother now.


It may be helpful to specify the library location when we install/update packages, when we're running on the cluster. I am using some centrally installed packages in a lib location I don't have write permissions for. E.g.:
```
update.packages(repos=repos, 
                ask=FALSE, 
                lib.loc="/home/jayoung/R/x86_64-unknown-linux-gnu-library/2.13", 
                instlib="/home/jayoung/R/x86_64-unknown-linux-gnu-library/2.13")
```

Can also load packages from a specific location:
```
library(ShortRead,lib.loc="/home/jayoung/R/x86_64-unknown-linux-gnu-library/2.13")

```


To REMOVE a package:
```
remove.packages("subplex")
remove.packages("checkmate", "/fh/fast/malik_h/grp/malik_lab_shared/linux_gizmo/R_devel_packages")
```

## General notes on where packages get installed 

`.libPaths()` shows the locations where R looks for packages

To add another location from within R:  `.libPaths("/home/btrask/traskdata/lib_linux/R")`

And/or we can use the linux `R_LIBS_USER` environmental variable. As of Dec 11 2025, I do that in my bash login file (`~/gizmo_login_items.bash`) using this line:
```
export R_LIBS_USER=$HOME/malik_lab_shared/linux_gizmo/R_packages/v_fh_4.4.0-foss-2023b
```
(I also used to have aliases in there called `Rnew` and `Rscript` but I got rid of those)

When I am running gizmo-Rstudio-R-4.4.0, this is the output of `.libPaths()`:
```
[1] "/fh/fast/malik_h/grp/malik_lab_shared/linux_gizmo/R_packages/v_fh_4.4.0-foss-2023b"
[2] "/app/software/fhR/4.4.0-foss-2023b"                                                
[3] "/home/jayoung/R/x86_64-pc-linux-gnu-library/4.4"                                   
[4] "/app/software/R-Tidyverse/4.4.0-gfbf-2023b"                                        
[5] "/app/software/R/4.4.0-gfbf-2023b/lib/R/library"  
```


To find out a package version, including details on where it was loaded from:
```
installed.packages()["R.utils","Version"]
installed.packages()["seqLogo",]
```


# Packages that need additional modules

There are some packages that won't install within an rhino/gizmo Rstudio-server session, because we need to load additional modules first, e.g. flextable and gdtools require the cairo module.

Here's what I'm trying, from a gizmo command line session
```
module purge
module load cairo/1.18.0-GCCcore-13.3.0
module load fhR/4.4.0-foss-2023b

R

install.package("gdtools")
install.package("flextable")
q(save="no")

module purge
```
Seems like that worked, and I can load those libraries within Rstudio-server without needing to load the cairo module.


# Notes on installation of some specific packages

## ggmsa install (Feb 2025)

when trying to install ggmsa in rhino/gizmo Rstudio, I got an error

the error was with one of the dependencies = proj4 package, needed by a different dependency called ggalt. Here's the proj4 error:

```
configure: error: Cannot find working proj.h headers and library.
*** You may need to install libproj-dev or similar! ***
```

instead I did a command-line install of proj4 - seems like it worked

```
module purge
module load PROJ/9.3.1-GCCcore-13.2.0
module load fhR/4.4.0-foss-2023b

which Rscript
    /app/software/R/4.4.0-gfbf-2023b/bin/Rscript

printenv R_LIBS_USER
    /home/jayoung/malik_lab_shared/linux_gizmo/R_packages/v_fh_4.4.0-foss-2023b

Rscript -e 'install.packages("proj4", repos="https://cran.r-project.org")'
module purge
```

after that, from with Rstudio:
```
install.packages("ggalt")
BiocManager::install("ggmsa")
```

this worked.


(xxx havent' yet tested the package!)



## Rmpi install (very old notes!)

wants sprng.h and mpi.h (sprng is for the dependency package rsprng)

```
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib64
install.packages("Rmpi", dependencies=TRUE)

install.packages("Rmpi", dependencies=TRUE, configure.args="--with-mpi=/usr/lib64/mvapich2/1.4-gcc/include") 

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib64:/usr/lib64/mvapich2/1.4-gcc/lib
install.packages("Rmpi", configure.args="--with-mpi=/usr/lib64/mvapich2/1.4-gcc") 
```

