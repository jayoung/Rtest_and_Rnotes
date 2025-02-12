# ggmsa install (Feb 2025)

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