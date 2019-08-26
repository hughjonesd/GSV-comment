
  <!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hughjonesd/GSV-comment/binder?urlpath=rstudio)
  <!-- badges: end -->
  
This is source code and data for a comment on
Garbarino, Ellen, Robert Slonim, and Marie Claire Villeval (2018):
[“A Method to Estimate Mean Lying Rates and Their Full
Distribution.”](https://doi.org/10.1007/s40881-018-0055-4) Journal of the Economic Science Association. 

For the PDF, see https://github.com/hughjonesd/GSV-comment/raw/master/GSV-comment-brief.pdf.

The published version is available [here](https://link.springer.com/article/10.1007/s40881-019-00069-x), if you like 
bad formatting and supporting parasitic academic publishers.

To reproduce the paper, click the "build binder" badge above. 

Or, to reproduce manually, clone the repository from your command line:

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

Code for the Bayesian methods is now available as a [standalone R package](https://github.com/hughjonesd/truelies).


`GSV-sims.csv` contains input for the GSV Java program. `GSV-comment-brief.Rmd` will automatically
create this, but you can use it if you want to skip that step.
