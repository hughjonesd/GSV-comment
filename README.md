
This is source code and data for a comment on Garbarino, Slonim and Villeval (JESA 2018).

To reproduce the paper, clone the repository from your command line:

```bash
git clone https://github.com/hughjonesd/GSV-comment.git
```

You'll also need GSV's original [lying calculator](http://lyingcalculator.gate.cnrs.fr/), and
various R libraries, especially the `checkpoint` library to reproduce the state of CRAN in 
November 2018.

When you're ready, fire up your favourite R IDE and:

1. Set `rerun_java <- TRUE` at the top of `GSV-comment-brief.Rmd`, and knit the document. 
2. R will stop and ask you to run the lying calculator on the file `"GSV-sims.csv"`. Do so. 
3. Set `rerun_java <- FALSE` and knit the document again.  

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

`GSV-heads.R` contains functions to estimate the distribution of *lies told* in a sample,
using GSV's methods. 

`bayesian-heads.R` contains some obsolete code.

`GSV-sims.csv` contains input for the GSV Java program. `GSV-comment-brief.Rmd` will automatically
create this, but you can use it if you want to skip that step.
