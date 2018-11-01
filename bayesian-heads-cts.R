
# Bayesian ways to calculate the probability of people in a sample reporting heads 
# when they observe tails

is_prob <- function (p) all(p >= 0 & p <= 1)

#' Calculate the probability of observing `heads`
#' good outcomes out of `N` total outcomes, conditional
#' on subjects lying with prob `lambda` if they see heads.
#'
#' @param lambda Prob of lying. 
#' @param heads Number of heads reported
#' @param N Total number in sample
#' @param P Probability of *bad* outcome
#'
#' @return A vector of probabilities
#'
#' @examples
prob_report_given_lambda <- function (lambda, heads, N, P) {
  stopifnot(is_prob(lambda), 0 <= heads, heads <= N, is_prob(P))

  # on average, heads are reported 1-P + P*lambda of the time
  dbinom(heads, N, 1 - P + lambda * P)
}


try_integral <- function(f, a, b, npoints = 100) {
  integ <- stats::integrate(f, a, b, subdivisions = npoints) # numerical integration
  if (integ$message != "OK") stop("stats::integrate failed with message: ", integ$message)
  
  return(integ$value)
}

#' Calculate posterior of lambdas given prior, and heads reports of good
#' outcome
#'
#' @inherit prob_report_given_lambda params
#' @param prior Prior over lambda. A function which takes a vector of values between 0 and 1, 
#'   and returns the probability density
#' @param npoints How many points to integrate on?
#'
#' @return A vector of posterior probabilities
#'
#' @examples
update_prior <- function(heads, N, P, prior, npoints = 1e3) {
  stopifnot(heads <= N, heads >= 0, is_prob(P), is.function(prior))

  f <- function (lprime) prior(lprime) * prob_report_given_lambda(lprime, heads, N, P)
  integ <- try_integral(f, 0, 1)
  denominator <- integ
  
  posterior <- function(lambda) {
    prob_report_given_lambda(lambda, heads, N, P) * prior(lambda)/denominator
  }
  
  return(posterior)
}

mean_dist <- function (dist) {
  stopifnot(is.function(dist))
  
  f <- function (x) x * dist(x)
  EV <- try_integral(f, 0, 1)

  return(EV)
}

interval_prob <- function (dist, l, r) {
  stopifnot(l <= r, is.function(dist))
  try_integral(dist, l, r)
}
