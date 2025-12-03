# Using renv for reproducibility

We use the `renv` package within each Rproject to manage package versions and allow reproducibility

When you use Rstudio to initialize a new project you can tell it to (a) use `renv` and (b) start an associated git repo.

When you start a fresh project using `renv` it doesn't see ANY of the packages we installed system-wide.

## Specific notes on use

### Mac paths

In a fresh project called test_Rpackage_versions, these are the .libPaths():
```
.libPaths()
[1] "/Users/jayoung/Documents/local_files/Rprojects/test_Rpackage_versions/renv/library/macos/R-4.5/aarch64-apple-darwin20"
[2] "/Users/jayoung/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.5/aarch64-apple-darwin20/4cd76b74"   
```

## renv resources


[Introduction to renv](https://rstudio.github.io/renv/articles/renv.html)

Using `renv` with [Bioconductor](https://rstudio.github.io/renv/articles/bioconductor.html)

[Posit guide](https://docs.posit.co/ide/user/ide/guide/environments/r/renv.html) to renv

https://erikgahner.dk/2025/using-renv-in-r/

https://carlosivanrodriguez.com/how-to/setup_renv.html

https://lmu-osc.github.io/introduction-to-renv/





those can go here or perhaps I move that note to the R repo
https://github.com/jayoung/thoughts/learning_to_do.md
