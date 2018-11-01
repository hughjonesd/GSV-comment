
#' Calculate the probability of observing `heads`
#' good outcomes out of `N` total outcomes, conditional
#' on `L` subjects being liars.
#'
#' @param L Number of Liars. Can be a vector.
#' @param heads Number of heads reported
#' @param N Total number in sample
#' @param P Probability of *bad* outcome
#'
#' @return
#'
#' @examples
#' prob_report_given_L(0:3, 2, 3, 0.5)
prob_report_given_L <- function (L, heads, N, P) {
  stopifnot(all(L <= N), heads <= N)
  # can't be true if L > heads: liars always report heads
  result <- rep(NA_real_, length(L))
  result[L > heads] <- 0
  # if heads reported out of N, and L are liars
  # then N-L people flipped a coin; of them heads - L got heads.
  # coin flip is binomially distributed.
  # the chance of this is the probability of getting heads - L
  
  low_L <- L[L <= heads]
  result[L <= heads] <- sapply(low_L, function (l) {
    dbinom(heads - l, size = N - l, prob = 1 - P)
  })

  return(result)
}


#' Calculate posterior of Ls given prior, and heads reports of good
#' outcome
#'
#' @param L Number of liars. Can be a vector
#' @param heads Number of heads reported
#' @param N Total in sample.
#' @param P Probability of *bad* outcome
#' @param prior Prior over L. Must be a N+1 length vector that sums to 1. 
#'   `prior[1]` is the prior probability that L == 0
#'   `prior[2]` is the prob that L == 1, etc.
#'
#' @return A vector of posterior probabilities
#'
#' @examples
#' probs_0_to_10 <- prob_L_given_report(0:10, 13, 20, 0.5, uniform_prior(20))
prob_L_given_report <- function(L, heads, N, P, prior) {
  check_posterior(prior)
  stopifnot(length(prior) == N + 1)
  stopifnot(L <= N, L >= 0, heads <= N, heads >= 0, P <= 1, P >= 0)

  evidence <- prob_report_given_L(0:N, heads, N, P)
  numerator <- prior[L + 1] * evidence[L + 1]
  denominator <- sum(prior*evidence)
  
  return(numerator/denominator)
}

#' Uniform prior on 0 to N
#' 
#' Creates a prior with equal probability of every number from 0 to N inclusive.
#'
#' @param N sample size
#'
#' @return The prior
#'
#' @examples
#' uniform_prior(5)
uniform_prior <- function (N) {
  rep(1/(N + 1), N + 1)
}

#' Update prior distribution given heads/N heads reported
#' 
#'
#' @inherit prob_L_given_report params
#' @param prior By default, a uniform prior on there being between 
#'   0 and N liars inclusive.
#'
#' @return The updated posterior
#'
#' @examples
#' updated <- update_prior(33, 50, 0.5)
#' posterior_mean(updated)
#' posterior_quantile(c(.025, .975), updated)
#' barplot(updated)
#' 
#' # different prior:
#' updated <- update_prior(33, 50, 0.5, prior = dbinom(0:50, 50, 0.25))
#' posterior_mean(updated)
#' posterior_quantile(c(.025, .975), updated)
update_prior <- function(heads, N, P, prior = uniform_prior(N)) {
  pstr <- prob_L_given_report(0:N, heads = heads, N = N, P = P,
        prior = prior)
  
  check_posterior(pstr)
  
  return(pstr)
}


naive_prob_L <- function (heads, N, P) {
  stopifnot(0 <= heads, heads <= N, 0 <= P, P <= 1)
  
  naive_prob <- (heads/N - (1- P))/P
  if (naive_prob <= 0) {
    stop("Naive estimate gives non-positive proportion of liars.", 
          " Can't use this for an empirical prior.")
  }
  if (naive_prob >= 1) {
    stop("Naive estimate gives proportion of liars >= 1.",
          "Can't use this for an empirical prior.")
  }
  
  return(naive_prob)
}

#' Calculate full distribution of 
#'
#' @param heads1 Number of reported heads in group 1
#' @param N1 Sample size in group 1
#' @param heads2 Number of reported heads in group 2
#' @param N2 Sample size in group 2
#' @param P  Probability of a bad outcome
#' @param prior_prob_L1 Prior belief that a person in group 1 is a liar.
#'   The default pools the groups and takes the naive estimate
#'   (heads/N - (1-P))/P.
#' @param prior_prob_L2 Prior belief that a person in group 2 is a liar.
#'   Shared with group 1 by default
#'
#' @return 
#' A data frame with two columns, "diff" and "prob". "diff"
#' gives a possible difference between the proportion of liars
#' in group 1 and group 2 (positive values mean that group 1 has more).
#' "prob" gives the posterior probability of that outcome.
#' 
#' @examples
#' GB_RU <- compare_groups(46, 89, 71, 100, 0.5)
#'
#' cdf <- cumsum(GB_RU$prob)
#' plot(GB_RU$diff, cdf, type = "s", xlim = c(-.25, .25))
#' 
#' # 95% confidence intervals:
#' ci_low <- max(which(cdf < 0.025))
#' ci_hi  <- min(which(cdf > 0.975))
#' GB_RU$diff[c(ci_low, ci_hi)]
#' 
#' # mean difference:
#' sum(GB_RU$diff * GB_RU$prob)
#' 
#' # compare:
#' prop_L <- (46 + 71)/(89 + 100) * 2 - 1
#' pstr_GB <- update_prior(46, 89, 0.5, dbinom(0:89, 89, prop_L))
#' pstr_RU <- update_prior(71, 100, 0.5, dbinom(0:100, 100, prop_L))
#' posterior_mean(pstr_GB)/89 - posterior_mean(pstr_RU)/100
compare_groups <- function (heads1, N1, heads2, N2, P, 
     prior_prob_L1 = naive_prob_L(heads1 + heads2, N1 + N2, P),
     prior_prob_L2 = prior_prob_L1) {
  stopifnot(0 <= heads1, heads1 <= N1, 0 <= heads2, heads2 <= N2, 
      0 <= P, P <= 1, 0 <= prior_prob_L1, prior_prob_L2 <= 1,
      0 <= prior_prob_L2, prior_prob_L2 <= 1)
  
  prior1 <- dbinom(0:N1, N1, prior_prob_L1)
  prior2 <- dbinom(0:N2, N2, prior_prob_L2)
  
  pstr1 <- update_prior(heads1, N1, P, prior1)
  pstr2 <- update_prior(heads2, N2, P, prior2)
  
  probs_diffs <- pstr1 %*% t(pstr2)
    # for L1 = 0:N1 and for L2= 0:N2, work out difference in prob of liars
  diffs <- outer(0:N1/N1, 0:N2/N2, "-")
  
  diffs <- c(diffs)  
  probs_diffs <- c(probs_diffs)
  check_posterior(probs_diffs)
  
  probs_diffs <- probs_diffs[order(diffs)]
  diffs <- sort(diffs)
  
  # sum probs where the value of diffs are equal:
  # Use rle to define "different values" (tapply would turn it into
  # a factor, rounding after 3 dec places)
  rled <- rle(diffs)
  unique_diffs <- rled$values
  rled$values <- seq_along(rled$values) # ensure all values have own label
  same_diff_values <- inverse.rle(rled)
  probs_unique_diffs <- tapply(probs_diffs, same_diff_values, sum)
  
  dist <- data.frame(
          diff = unique_diffs, 
          prob  = probs_unique_diffs
        )
  return(dist)
}


check_posterior <- function(pstr) {
  stopifnot(isTRUE(all.equal(sum(pstr), 1)), all(pstr >= 0))
}

#' Returns the mean of a discrete distribution
#'
#' @param pstr A vector of positive numbers summing to 1
#'
#' @return
#'
#' @examples
posterior_mean <- function (pstr) {
  check_posterior(pstr)
  
  sum(pstr * seq(0, along.with = pstr))
}


#' Returns quantiles of a discrete distribution
#'
#' @param q A vector of quantiles
#' @param pstr A vector of positive numbers summing to 1
#'
#' @return A vector of values at the quantiles
#'
#' @examples
posterior_quantile <- function (q, pstr) {
  stopifnot(all(q >= 0), all(q <= 1))
  check_posterior(pstr)
  cdf <- cumsum(pstr) # prob L <= l, for l = 0 to N

  # so e.g. 30% quantile is x of 0:N s.t. cdf[x] < q < cdf[x+1]
  # the -1 returns the value of Lstarting from 0.
  L_quantiles <- findInterval(q, c(0, unique(cdf))) - 1
 
  return(L_quantiles) 
}

#' Checks coverage of confidence intervals by simulations
#'
#' @param L True number of liars
#' @param N Sample size
#' @param P Probability of bad outcome
#' @param CI Confidence interval between 0 and 1
#' @param nreps Number of simulations to run
#' @param prior_fn A one-argument function which returns a N+1 length vector of probabilities
#'
#' @return
#'
#' @examples
#' check_ci_coverage(15, 30, P = 0.5)
#' check_ci_coverage(15, 30, P = 0.5, CI = 0.99)
check_ci_coverage <- function (L, N, P, CI = 0.95, nreps = 1000, prior_fn = uniform_prior) {
  tail <- (1 - CI)/2
  stopifnot(L <= N, P <= 1, P >= 0, CI <= 1, CI >= 0, is.function(prior_fn))
  
  within_ci <- replicate(nreps, {
    heads <- L + rbinom(1, N - L, prob = 1 - P)
    pstr <- update_prior(heads, N, P, prior = prior_fn(N))
    bounds <- posterior_quantile(c(tail, 1 - tail), pstr)
    bounds[1] <= L && L <= bounds[2] 
  })
  
  return(mean(within_ci))
}
