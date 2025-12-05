# Using renv for reproducibility

We use the `renv` package within each Rproject to manage package versions and allow reproducibility

When you use Rstudio to initialize a new project you can tell it to (a) use `renv` and (b) start an associated git repo.

When you start a fresh project using `renv` it doesn't see ANY of the packages we installed system-wide.

## TO DO list

Read the resources

Questions:
- should I be recording something like a snapshot ID or git commit ID, when I have working code, before I change package versions?   Right now I'm treating the snapshots anonymously

## renv resources

[renv documentation site](https://rstudio.github.io/renv/)

[Introduction to renv](https://rstudio.github.io/renv/articles/renv.html)

Using `renv` with [Bioconductor](https://rstudio.github.io/renv/articles/bioconductor.html)

[Posit guide](https://docs.posit.co/ide/user/ide/guide/environments/r/renv.html) to renv

https://erikgahner.dk/2025/using-renv-in-r/

https://carlosivanrodriguez.com/how-to/setup_renv.html

https://lmu-osc.github.io/introduction-to-renv/


Example repo where I used `renv`: [test_Rpackage_versions](https://github.com/jayoung/test_Rpackage_versions)

## How to use renv

### Initialize a new project

Don't use the Rstudio way, because I want to specify that we are using Bioconductor, and to specify a particular version ([version info and release dates](https://bioconductor.org/about/release-announcements/)).

```
library("renv")
renv::init(bioconductor = "3.22")
```

Status is useful:
```
renv::status()
```

Initialization does a lot of things:
- creates `renv.lock` - a json style file that records the versions of all packages used. Updated every time we run `renv::snapshot()`
- creates `renv` directory - package management will occur in here. It contains its own `.gitignore` file that means the packages themselves don't get synced to git, only the metadata about versions.
- creates `.Rprofile` file, which will source a long script called `renv/activate.R` every time R is restarted for this project


### If I cloned a project from git into a new place

It won't have the packages installed, but it will know what it needs.  

I can run `renv::restore()` and it will install the necessary packages

### Installing packages using renv

See [install docs](
https://rstudio.github.io/renv/reference/install.html).

To install a CRAN package I do this:
```
renv::install("rmarkdown")
```

To install a Bioconductor package I do this:
```
renv::install("bioc::ape")
```

In a DIFFERENT Rproject that's not controlled by `renv` the libPaths are unchanged, and it's still using the older package versions I had previously installed. That's good.
```
.libPaths()
[1] "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
```

### Updating packages

```
renv::update()
```


### Lock the environment


When I have a version of the code+packages that works, I "lock" that setup: 

```
renv::snapshot()
```

If things are messed up and I need to go to the locked setup I restore it:
```
renv::restore()
```

Restoring to older versions can also be done - see `renv::history()` and `renv::revert()`.

## Specific notes on use


### Bioconductor

Bioconductor advice [here](https://rstudio.github.io/renv/articles/bioconductor.html).

We initialize `renv` in a different way.

### Recommended habits

When I start a new project, I'll use renv, although I'll initialize from command-line NOT from Rstudio (so I can specify NBioconductor use).  

I'll start with the most recent version of R and packages that I can. Perhaps I also run `update.packages()` ASAP.

Each time we have a version of the code that works well, we should "lock" that setup:
- call `renv::snapshot()` to record the latest package versions in your lockfile.
- git add/commit/push to commit the lockfile



### Mac paths

System-wide cache of packages goes here:
```
~/Library/Caches/org.R-project.R/R/renv
```
That location includes a dir called `cache/v5/macos/R-4.5/aarch64-apple-darwin20/`, which contains a folder for each package. Within each package dir there's a subdir with the version ID.

Project-specific packages: in the project dir, there's a folder called `renv/library/macos/R-4.5/aarch64-apple-darwin20`, within which is a link for each package, linking to the correct version of the package in `/Users/jayoung/Library/Caches/org.R-project.R/R/renv/cache/v5/macos/R-4.5/aarch64-apple-darwin20`. 


In a fresh project called test_Rpackage_versions, these are the .libPaths():
```
.libPaths()
[1] "/Users/jayoung/Documents/local_files/Rprojects/test_Rpackage_versions/renv/library/macos/R-4.5/aarch64-apple-darwin20"
[2] "/Users/jayoung/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74"   
```

i.e. a named folder within the system-wide `renv/sandbox/` dir - this is where the standard base-R packages are, and the project's `renv/library/` dir, which is where additional packages are.



## Some tests

### Mac - test_Rpackage_versions 

On Mac, in `/Users/jayoung/Documents/local_files/Rprojects/`.

I initialized this using Rstudio in a new folder, checking the `use renv` and `start git repo` options.


### Mac - test_Rpackage_versions_v2

On Mac, in `/Users/jayoung/Documents/local_files/Rprojects/`.

I initialized this using Rstudio in a new folder, checking the `start git repo` option but NOT the `use renv` button, because I want to tell it to use Bioconductor when I initialize.

In this project, I explicitly initialized `renv` like this, telling it a specific version:

```
library("renv")
# use a specific version of Bioconductor
renv::init(bioconductor = "3.22")

update.packages()
```

Installed some packages:
```
renv::install("bioc::ape")
```

Locked the setup:
```
renv::snapshot()
```

Installed more packages

```
renv::install("ggplot2")
renv::install("tidyverse")
```

gdtools is needed for ggtree. It was a big pain to install. It was trying to compile from scratch rather than using binary, and that was failing due to dependencies. 

Fixed one dependency problem by using command-line `brew` to install `cairo` (see [notes](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/janets_NOTES_forMyself/local_infrastructure/mac_laptop_2024_setup_NOTES.md)) - I restart R before doing more, hoping it sees where cairo is. 

Then I got stuck for a while on another error regarding freetype/freetype2/ fontconfig/fontconfig.h

Here's my pak attempt for gdtools (didn't work)
```
# install pak
renv::install("pak")

# Configure `renv` to use `pak` 
options(renv.config.pak.enabled = TRUE)

# try installing gdtools  - it still fails
renv::install("gdtools")

# I got errors like this:
/opt/homebrew/Cellar/cairo/1.18.4/include/cairo/
    In file included from src/tests/sysdeps.c:1:
    /opt/homebrew/include/cairo/cairo-ft.h:50:10: fatal error: 'fontconfig/fontconfig.h' file not found
# I can see that file here:
/opt/homebrew/Cellar/fontconfig/2.17.1/include/fontconfig/fontconfig.h 
```

To keep things clean I'll separately install two other packages it tries to install alongside gdtools - this worked fine:
```
pak::pkg_install("digest")
pak::pkg_install("systemfonts")
```

`renv::install("gdtools")` still failed (same error)

It DOES work when I force it to install from the binary by supplying the URL (only when `renv.config.pak.enabled = FALSE`)

```
options(renv.config.pak.enabled = FALSE)
renv::install("https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.5/gdtools_0.4.4.tgz")
```

Now install ggtree - seems like it worked this time
```
renv::install("ggtree")
```

Locked the setup:
```
renv::snapshot()
```


The first time I do a git push, I need to create a new remote repo and tell it where to push to.
```
git remote add origin git@github.com:jayoung/test_Rpackage_versions_v2.git
     # seems like that worked
```
I initially couldn't do a git push (`ERROR: Repository not found. fatal: Could not read from remote repository.`). I used github website to create a new blank repo called test_Rpackage_versions_v2, and THEN the git push worked.
```
git push --set-upstream origin main
```



## Other maybe useful tools

`BiocManager::valid()` gives a report on versions of all installed Bioc packages
