
# Bayesian ways to calculate the probability of people in a sample reporting heads 
# when they observe tails


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



#' Estimate probability that two samples s1, s2 from 2 independent distributions have s1 > s2
#'
#' @param dist1 
#' @param dist2 
#'
#' @return A probability scalar
#' @export
#'
#' @examples
compare_dists <- function (dist1, dist2) {
  outer_f <- function (l1) {
    dist1(l1) * dist_prob_within(dist2, 0, l1)
  }
  outer_f <- Vectorize(outer_f)
  prob_1_more <- try_integral(outer_f, 0, 1)
  
  return(prob_1_more)
}


#' Given dist1 and dist2, find the pdf of dist1 - dist2
#' 
#' At the moment this only works when dist1 and dist2 are defined on [0, 1].
#' 
#' @param dist1,dist2 Probability distribution functions 
#'
#' @return A probability distribution function defined on [-1, 1].
#' 
#' @examples
difference_dist <- function (dist1, dist2) {
  
  # prob l1 - l2 = d is integral from 0 to 1-d of dist1(z+d)*dist2(z)
  # l1 - l2 = -.5; l1 = l2 + .5; integral from .5 to 1 of dist1(z + d)*dist2(z)
  dd_fun <- function (d) {
    stopifnot(d >= -1 && d <= 1)
    dist_prob <- function(x) dist1(x + d) * dist2(x)
    dist_prob <- Vectorize(dist_prob)
    ends <- if (d > 0) c(0, 1 - d) else c(-d, 1)
    try_integral(dist_prob, ends[1], ends[2])
  }
  
  return(Vectorize(dd_fun))
}


#' Find mean of a probability distribution function
#'
#' @param dist A one-argument function
#' @param l Minimum value
#' @param r Maximum value
#'
#' @return A scalar
#'
#' @examples
dist_mean <- function (dist, l = 0, r = 1) {
  stopifnot(is.function(dist))
  
  f <- function (x) x * dist(x)
  EV <- try_integral(f, l, r)

  return(EV)
}

dist_prob_within <- function (dist, l, r) {
  stopifnot(l <= r, is.function(dist))
  try_integral(dist, l, r)
}


#' Find quantiles given a probability distribution function
#'
#' @param dist A one argument function
#' @param probs A vector of probabilities
#'
#' @return A vector of quantiles
#'
#' @examples
dist_quantile <- function(dist, probs) {
  stopifnot(is_prob(probs), is.function(dist))
  
  qs <- sapply(probs, function (prob) {
    f <- function (x) (dist_prob_within(dist, 0, x) - prob)^2
    # want x s.t. interval_prob(dist, 0, x) == prob
    res <- optimize(f, c(0, 1))
    res$minimum
  })
  
  return(qs)
}

is_prob <- function (p) all(p >= 0 & p <= 1)


if (FALSE) {
  N <- 100
  P <- 0.5
  prior <- dunif
  CI <- 0.95
  cov_check <- replicate(1000, {
    lambda <- runif(1)
    heads <- rbinom(1, N, p = lambda * P + 1 - P)
    posterior <- update_prior(heads, N, P, prior)
    cis <- dist_quantiles(posterior, c(1/2 - CI/2, 1/2 + CI/2))
    
    cis[1] <= lambda && lambda <= cis[2]
  })
  table(cov_check)
}
