
# functions to compute distributions of lies told in coin flip experiments, according to GSV (2018).
# WARNING: IMHO these methods are flawed, use at your own risk!


dist_L_gsv <- function(heads, N, P) {
  # L liars means that T = heads - L
  # T is distributed binomially conditional on being no more than `heads`
  L = 0:heads
  result <- rep(NA_real_, length(L))
  result <- dbinom(heads - L, N, 1 - P) / sum(dbinom(L, N, 1 - P))
  
  result
}

prop_liars <- function(L, heads, N) {
  stopifnot(all(L <= heads), heads >= 0, heads <= N)
  # for any number of liars L, reported heads, and N:
  # true heads T = heads - L
  # number observing low outcome = N - T
  # proportion lying is L/(N-T)
  T <- heads - L
  prop_lied <- L/(N - T)
  prop_lied[N == T] <- 1 # if everyone saw true heads, then   
  # we can't estimate how many might have
  # GSV code solves this by assuming 1
  
  prop_lied
}

prop_liars_ci_gsv <- function (ci, heads, N, P) {
  stopifnot(0 <= ci, ci <= 1, heads <= N, 0 <= heads, 0 <= P, P <= 1)
  
  dist_L <- dist_L_gsv(heads, N, P)
  stopifnot(length(dist_L) == heads + 1)
  
  pl <- prop_liars(0:heads, heads, N)
  stopifnot(length(pl) == heads + 1)
  
  pl <- pl[dist_L > 0] # get rid of 0-probability events
  dist_L <- dist_L[dist_L > 0]
  dist_L <- dist_L[order(pl)]
  pl <- sort(pl)
  
  cdf <- cumsum(dist_L)
  
  # we know pl is always positive, so represent that:
  cdf <- c(0, cdf)
  pl  <- c(0, pl)
  
  tail <- (1 - ci)/2
  # lower bound: no more than tail% of probability has less than this.
  # this ought to be max(which(tail <= cdf)) but the paper says this instead - which
  # will lead to excessively narrow confidence intervals when N is low
  lb_idx <- min(which(cdf > tail)) # must exist since cdf starts at 0.
  lb <- pl[lb_idx]
  
  # upper bound: no more than 1-tail% of probability has more. This is what the paper says
  # but it is far from what they do:
  ub_idx <- max(which(cdf < 1 - tail))
  ub <- pl[ub_idx]
  
  return(c(lb, ub))
}

prop_liars_ev_gsv <- function(heads, N, P) {
  dist_L <- dist_L_gsv(heads, N, P)
  prop_lied <- prop_liars(0:heads, heads, N)
  sum(dist_L * prop_lied)
}
