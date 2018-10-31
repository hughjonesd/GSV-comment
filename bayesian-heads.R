
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

check_posterior <- function(pstr) {
  stopifnot(isTRUE(all.equal(sum(pstr), 1)), all(pstr >= 0))
}

#' Returns the mean of a discrete distribution
#'
#' @param pstr A vector of positive numbers summing to 1
#'
#' @return
#' @export
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
#' @export
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
#' @export
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
