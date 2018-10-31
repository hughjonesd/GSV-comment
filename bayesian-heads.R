
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
#' @return
#' @export
#'
#' @examples
prob_L_given_report <- function(L, heads, N, P, prior) {
  check_posterior(prior)
  stopifnot(length(prior) == N + 1)
  stopifnot(L <= N, L >= 0, heads <= N, heads >= 0, P <= 1, P >= 0)

  evidence <- prob_report_given_L(0:N, heads, N, P)
  numerator <- prior[L + 1] * evidence[L + 1]
  denominator <- sum(prior*evidence)
  
  return(numerator/denominator)
}

uniform_prior <- function (N) rep(1/(N + 1), N + 1)


#' Update prior distribution given heads/N heads reported
#' 
#'
#' @inherit prob_L_given_report params
#' @param prior By default, a uniform prior on there being between 
#'   0 and N liars inclusive.
#'
#' @return
#' @export
#'
#' @examples
update_prior <- function(heads, N, P, prior = uniform_prior(N)) {
  pstr <- prob_L_given_report(0:N, heads = heads, N = N, P = P,
        prior = prior)
  
  check_posterior(pstr)
  
  return(pstr)
}

check_posterior <- function(pstr) {
  stopifnot(isTRUE(all.equal(sum(pstr), 1)), all(pstr >= 0))
}

posterior_mean <- function (pstr) {
  check_posterior(pstr)
  
  sum(pstr * seq(0, along.with = pstr))
}

posterior_quantile <- function (q, pstr) {
  check_posterior(pstr)
  cdf <- cumsum(pstr) # prob L <= l, for l = 0 to N

  # so e.g. 30% quantile is x of 0:N s.t. cdf[x] < q < cdf[x+1]
  # the -1 returns the value of Lstarting from 0.
  L_quantiles <- findInterval(q, c(0, unique(cdf))) - 1
 
  return(L_quantiles) 
}

