# ChAMP
Clone of the Bioconductor package ChAMP, required to make bugfixes

# Installation

Install this package from this Github repo. This requires `devtools`:
```
install.packages('devtools')
library(devtools)
```

By default (on my Ubuntu setup), any warnings raised during installation of a remote package are converted to errors. This prevents installation for really innocuous reasons. To fix this, edit the `.Renviron` file (located in `$HOME`) to add/modify the following:
```
R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
```

Installation should then be as simple as
```
install_github("gaberosser/ChAMP")
```
