
This is source code and data for a comment on Garbarino, Slonim and Villeval (JESA 2018).

To reproduce the paper, clone the repository and run `GSV-comment-brief.Rmd` from within RStudio. 
To reproduce starting from scratch, you'll need GSV's original [lying calculator](http://lyingcalculator.gate.cnrs.fr/).

The file `bayesian-heads.R` contains functions to estimate the distribution of liars in a sample
using Bayesian methods. Here's an example:

```r
# default is a uniform prior
# heads is number of reported heads, N is sample size, P is probability of *bad* outcome
updated <- update_prior(heads = 33, N = 50, P = 0.5)
posterior_mean(updated)
posterior_quantile(c(.025, .975), updated)
barplot(updated)

```