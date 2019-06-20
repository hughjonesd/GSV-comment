
This is source code and data for a comment on
Garbarino, Ellen, Robert Slonim, and Marie Claire Villeval (2018):
[“A Method to Estimate Mean Lying Rates and Their Full
Distribution.”](https://doi.org/10.1007/s40881-018-0055-4) Journal of the Economic Science Association. 

For the PDF, see https://github.com/hughjonesd/GSV-comment/raw/master/GSV-comment-brief.pdf.

The published version is available [here](https://link.springer.com/article/10.1007/s40881-019-00069-x), if you like 
bad formatting and supporting parasitic academic publishers.

To reproduce the paper, clone the repository from your command line:

```bash
git clone https://github.com/hughjonesd/GSV-comment.git
```

You'll also need GSV's original [lying calculator](http://lyingcalculator.gate.cnrs.fr/), and
the `checkpoint` R package.

When you're ready, fire up your favourite R IDE and:

1. On the command line, run:

```r
library(checkpoint)
checkpoint("2018-11-03")
```

This will install the versions of the packages used to create the comment.

2. Set `rerun_java <- TRUE` at the top of `GSV-comment-brief.Rmd`, and knit the document.
3. R will stop and ask you to run the GSV lying calculator on the file `"GSV-sims.csv"`. Do so.
4. Set `rerun_java <- FALSE` and knit the document again.

## Other files

`bayesian-heads-cts.R` contains functions to estimate the distribution of liars in a sample
using Bayesian methods. Here's an example:

```r
source("bayesian-heads-cts.R")
# heads is number of reported heads, N is sample size, P is probability of *bad* outcome
updated <- update_prior(heads = 33, N = 50, P = 0.5, prior = dunif)
dist_mean(updated)
dist_quantile(updated, c(.025, .975))
dist_hdr(updated)
curve(updated)
```

*Update*: code is now available as a [standalone R package](https://github.com/hughjonesd/truelies).

`GSV-heads.R` contains functions to estimate the distribution of *lies told* in a sample,
using GSV's methods.

`bayesian-heads.R` contains some obsolete code.

`GSV-sims.csv` contains input for the GSV Java program. `GSV-comment-brief.Rmd` will automatically
create this, but you can use it if you want to skip that step.
