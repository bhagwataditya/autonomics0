[![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/0.1.0/wip.svg)](http://www.repostatus.org/#wip)

# autonomics.support

Miscellaneous utility functions.  Part of [autonomics](https://bitbucket.org/account/user/graumannlabtools/projects/autonomics), the R suite for automated omics data analysis.

## Installation

To install the package, you first need the
[*devtools*](https://github.com/hadley/devtools) package.

```{r}
install.packages("devtools")
```

Then you can install the *autonomics.support* package using

```{r}
library(devtools)
install_bitbucket("graumannlabtools/autonomics.support")
```

## Functionality

`factorify` creates a factor where the levels are in the order that they appeared in the original vector.

`uniquify` makes a vector unique by appending spaces to the elements.

`cfread` wraps `data.table::fread` with `suppressWarnings` and some convenient defaults.

`print2pdf` wraps `grDevices::pdf`, closing the device on completion.

`print2txt` wraps `utils::write.table` with some convenient defaults.